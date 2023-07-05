import asyncio
import numpy as np
from Bio import SeqIO
import gzip
from collections import Counter
from typing import Callable
import os
import pickle
import taichi as ti
from itertools import chain
from typing import List

"""
Author: Lu Cheng, @chengl7
Date: July 5, 2023
Description: kmer counting script
"""


MISSING_VAL = 255

ti.init(arch=ti.cuda, default_ip=ti.i64)
if ti.cfg.arch == ti.cuda:
    print("GPU is available")
else:
    print("GPU is not available")


# create a directory if not exist
def mk_dir(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    else:
        print(f"Folder already exists: {folder_path}")


# remove all files in chunk_dir
def rm_files(chunk_dir="./chunks"):
    if os.path.exists(chunk_dir):
        files = os.listdir(chunk_dir)
        # Iterate over the files and remove them one by one
        for file in files:
            file_path = os.path.join(chunk_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)


def get_chunk_file_paths(chunk_dir="./chunks"):
    assert os.path.exists(chunk_dir)
    files = [file for file in os.listdir(chunk_dir) if file.endswith(".pkl")]
    files = sorted(files, key=lambda x: int(x.split("_")[1][:-4]))
    files = [os.path.join(chunk_dir, file) for file in files]
    return files


def get_taichi_dtype(np_dtype: np.uint32):
    dict = {np.uint8: ti.uint8, np.uint16: ti.uint16, np.uint32: ti.uint32, np.uint64: ti.uint64}
    return dict[np_dtype]


class Buffer:
    def __init__(self, buffer_size: int = 2 ** 26, dtype=np.uint8, id=None):
        # 2**26 bytes are 64MB, 2**30 bytes are 1GB
        self.buffer_size = buffer_size
        self.dtype = dtype
        self.buffer = np.empty(buffer_size, dtype=dtype)
        self.pointer = 0
        self.data_to_write = None
        self.is_full = False
        self.id = None

    def append(self, data: np.uint8):
        # append input data into the buffer, return flag of success
        assert data.dtype == self.buffer.dtype
        data_size = len(data)
        buffer_space = self.buffer_size - self.pointer

        # Check if the buffer has enough space to accommodate the new data
        if data_size > buffer_space:
            self.data_to_write = data
            self.is_full = True
            return False
        else:
            # Append the data to the buffer
            self.buffer[self.pointer: (self.pointer + data_size)] = data
            self.pointer += data_size
            return True

    def flush(self):
        # flush the buffer such that it can be used again, data_to_write needs to be handled first
        assert self.data_to_write is None
        self.pointer = 0
        self.data_to_write = None
        self.is_full = False
        self.id = None


class MaxSizeQueue:
    def __init__(self, maxsize=0):
        self.maxsize = maxsize
        self.queue = asyncio.Queue(maxsize=maxsize)
        self.op_lock = asyncio.Condition()  # lock for operations on the queue

    async def put(self, item):
        async with self.op_lock:
            while self.queue.full():
                await self.op_lock.wait()
            await self.queue.put(item)
            self.op_lock.notify_all()

    async def get(self):
        async with self.op_lock:
            while self.queue.empty():
                await self.op_lock.wait()
            item = await self.queue.get()
            self.op_lock.notify_all()
            return item

    def full(self):
        return self.queue.full()

    def empty(self):
        return self.queue.empty()


def dna2arr(dna_str, dtype=np.uint8) -> np.ndarray:
    """
    convert an input DNA string to numpy uint8 array, with a missing value appended
    Args:
        dna_str: a DNA sequence, all letters should be upper case
        dtype: data type for storing DNA string
    Returns:
        a numpy array
    """
    res = np.empty(len(dna_str) + 1, dtype=dtype)
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, b in enumerate(dna_str):
        res[i] = base_map.get(b, MISSING_VAL)
    res[-1] = MISSING_VAL  # add a separator to the end of the string
    return res


def read_dnaseq_file(file_name, file_type="fasta") -> np.ndarray:
    """
    file_name: input DNA sequence file name
    file_type: fasta, fastq,
    """

    def read_stream(fh):
        for rec in SeqIO.parse(fh, file_type):
            yield dna2arr(str(rec.seq).upper())

    if file_name.endswith(".gz"):
        with gzip.open(file_name, "rt") as fh:
            yield from read_stream(fh)
    else:
        with open(file_name, "r") as fh:
            yield from read_stream(fh)


# get the hash dtype for given kmer length
def get_hash_dtype(kmer_len):
    if 0 < kmer_len < 16:
        return np.uint32
    elif kmer_len < 32:
        return np.uint64
    else:
        raise Exception(f"max_kmer_len=31, kmer_len={kmer_len} is greater the maximum value.")


# generate a hash mask for kmers such that bits out of scope can be masked to 0
def gen_hash_mask(k: int, dtype: Callable[[int], np.dtype]):
    mask = dtype((1 << 2 * k) - 1)
    return mask


@ti.func
def revcom_hash_uint32(in_hash: ti.u32,
                       mask: ti.u32,
                       twobit_mask: ti.u32,
                       k: int):
    com_hash = mask - in_hash  # complement hash
    ret_hash = twobit_mask & com_hash
    for i in range(k - 1):
        ret_hash = ret_hash << 2
        com_hash = com_hash >> 2
        ret_hash += twobit_mask & com_hash
    return ret_hash


@ti.func
def revcom_hash_uint64(in_hash: ti.u64,
                       mask: ti.u64,
                       twobit_mask: ti.u64,
                       k: int):
    com_hash = mask - in_hash  # complement hash
    ret_hash = twobit_mask & com_hash
    for i in range(k - 1):
        ret_hash = ret_hash << 2
        com_hash = com_hash >> 2
        ret_hash += twobit_mask & com_hash
    return ret_hash


@ti.kernel
def revcom_hash_kernel_uint32(in_hash_arr: ti.types.ndarray(dtype=ti.u32),
                              out_hash_arr: ti.types.ndarray(dtype=ti.u32),
                              mask_arr: ti.types.ndarray(dtype=ti.u32),
                              kmer_len: int, in_hash_arr_size: int):
    for i in range(in_hash_arr_size):
        out_hash_arr[i] = revcom_hash_uint32(in_hash_arr[i], mask_arr[0], mask_arr[1], kmer_len)


@ti.kernel
def revcom_hash_kernel_uint64(in_hash_arr: ti.types.ndarray(dtype=ti.u64),
                              out_hash_arr: ti.types.ndarray(dtype=ti.u64),
                              mask_arr: ti.types.ndarray(dtype=ti.u64),
                              kmer_len: int, in_hash_arr_size: int):
    for i in range(in_hash_arr_size):
        out_hash_arr[i] = revcom_hash_uint64(in_hash_arr[i], mask_arr[0], mask_arr[1], kmer_len)


def get_revcom_hash_arr(in_hash_arr: np.ndarray, kmer_len: int):
    hash_dtype = get_hash_dtype(kmer_len)
    mask_arr = np.array([(1 << 2 * kmer_len) - 1, 3], dtype=hash_dtype) # mask and twobit_mask

    out_hash_arr = np.empty_like(in_hash_arr)
    hash_arr_size = len(in_hash_arr)
    if hash_dtype == np.uint32:
        revcom_hash_kernel_uint32(in_hash_arr, out_hash_arr, mask_arr, kmer_len, hash_arr_size)
    elif hash_dtype == np.uint64:
        revcom_hash_kernel_uint64(in_hash_arr, out_hash_arr, mask_arr, kmer_len, hash_arr_size)
    return out_hash_arr


def revcom_hash(in_hash: np.uint64, kmer_len: int):
    hash_dtype = get_hash_dtype(kmer_len)
    in_hash = hash_dtype(in_hash)
    mask_arr = np.array([(1 << 2 * kmer_len) - 1, 3], dtype=hash_dtype) # mask and twobit_mask
    mask = mask_arr[0]
    twobit_mask = mask_arr[1]

    com_hash = mask - in_hash  # complement hash
    # reverse
    ret_hash = twobit_mask & com_hash
    for i in range(kmer_len - 1):
        ret_hash = ret_hash << hash_dtype(2)
        com_hash = com_hash >> hash_dtype(2)
        ret_hash += twobit_mask & com_hash
    return ret_hash


# get the hash value for invalid kmers
def get_invalid_hash(dtype: Callable[[int], np.dtype]):
    return dtype(np.iinfo(dtype).max)


# def kmer2hash(arr: np.ndarray, arr_size: int, st_pos: int, k: int,
#               hash_arr: np.ndarray, hash_dtype: Callable[[int], np.dtype], invalid_hash: np.dtype):
#     """
#     compute the hash key for kmer starting at st_pos
#     arr: input numpy array, np.uint8, elements are 0,1,2,3, MISSING_VAL=255
#     arr_size: maximum size of input arr
#     st_pos: start position of kmer
#     k: kmer length
#     hash_arr: numpy arr to store the hash
#     hash_dtype: numpy dtype of hash key
#     invalid_hash: invalid hash value, reserved for invalid kmers
#     return: a hash
#     """
#     hash_arr[st_pos] = invalid_hash
#     if st_pos + k >= arr_size:
#         # note that the last letter at arr_size-1 is MISSING_VAL
#         return invalid_hash
#
#     kh = hash_dtype(0)
#     for i in range(k):
#         if arr[st_pos + i] == MISSING_VAL:
#             # if any letter in the kmer is MISSING_VAL
#             return invalid_hash
#         kh = kh << hash_dtype(2)
#         kh += hash_dtype(arr[st_pos + i])
#     hash_arr[st_pos] = kh
#     return kh

@ti.func
def cal_ham_dist_uint32(hash1: ti.u32, hash2: ti.u32, kmer_len: int):
    xor_result = hash1 ^ hash2
    twobit_mask = ti.cast(3, ti.u32)
    hamming_dist = 0
    for _ in range(kmer_len):
        cmp_res = xor_result & twobit_mask
        hamming_dist += cmp_res != 0
        xor_result >>= 2
    return hamming_dist


@ti.kernel
def cal_ham_dist_kernel_uint32(hash_arr: ti.types.ndarray(dtype=ti.u32),
                               target_hash: ti.types.ndarray(dtype=ti.u32),
                               ham_dist_arr: ti.types.ndarray(dtype=ti.u8),
                               hash_arr_size: int,
                               kmer_len: int):
    for i in range(hash_arr_size):
        ham_dist_arr[i] = ti.cast(cal_ham_dist_uint32(hash_arr[i], target_hash[0], kmer_len), ti.u8)


@ti.func
def cal_ham_dist_uint64(hash1: ti.u64, hash2: ti.u64, kmer_len: int):
    xor_result = hash1 ^ hash2
    twobit_mask = ti.cast(3, ti.u64)
    hamming_dist = 0
    for _ in range(kmer_len):
        cmp_res = xor_result & twobit_mask
        hamming_dist += cmp_res != 0
        xor_result >>= 2
    return hamming_dist


@ti.kernel
def cal_ham_dist_kernel_uint64(hash_arr: ti.types.ndarray(dtype=ti.u64),
                               target_hash: ti.types.ndarray(dtype=ti.u64),
                               ham_dist_arr: ti.types.ndarray(dtype=ti.u8),
                               hash_arr_size: int,
                               kmer_len: int):
    for i in range(hash_arr_size):
        ham_dist_arr[i] = ti.cast(cal_ham_dist_uint64(hash_arr[i], target_hash[0], kmer_len), ti.u8)


@ti.func
def kmer2hash_taichi_uint32(arr: ti.types.ndarray(dtype=ti.u8), arr_size: int, st_pos: int, k: int,
                            hash_arr: ti.types.ndarray(dtype=ti.u32),
                            invalid_hash: ti.types.u32,
                            missing_val: ti.types.u32):
    # hash_arr[st_pos] = invalid_hash
    invalid_hash_flag = 0
    if st_pos + k >= arr_size:
        invalid_hash_flag = 1

    kh = ti.u32(0)
    for i in range(k):
        if arr[st_pos + i] == missing_val:
            invalid_hash_flag = 1
        kh = kh << 2
        kh += arr[st_pos + i]
    hash_arr[st_pos] = kh

    if invalid_hash_flag > 0:
        hash_arr[st_pos] = invalid_hash


@ti.kernel
def kmer2hash_kernel_uint32(arr: ti.types.ndarray(dtype=ti.u8), arr_size: int, k: int,
                            hash_arr: ti.types.ndarray(dtype=ti.u32),
                            invalid_hash_arr: ti.types.ndarray(dtype=ti.u32),
                            missing_val_arr: ti.types.ndarray(dtype=ti.u8)):
    for i in range(arr_size):
        kmer2hash_taichi_uint32(arr, arr_size, i, k, hash_arr, invalid_hash_arr[0], missing_val_arr[0])


@ti.func
def kmer2hash_taichi_uint64(arr: ti.types.ndarray(dtype=ti.u8), arr_size: int, st_pos: int, k: int,
                            hash_arr: ti.types.ndarray(dtype=ti.u64), invalid_hash: ti.u64,
                            missing_val: ti.u8):
    # hash_arr[st_pos] = invalid_hash
    invalid_hash_flag = 0
    if st_pos + k >= arr_size:
        invalid_hash_flag = 1

    kh = ti.u64(0)
    for i in range(k):
        if arr[st_pos + i] == missing_val:
            invalid_hash_flag = 1
        kh = kh << 2
        kh += arr[st_pos + i]
    hash_arr[st_pos] = kh

    if invalid_hash_flag > 0:
        hash_arr[st_pos] = invalid_hash


@ti.kernel
def kmer2hash_kernel_uint64(arr: ti.types.ndarray(dtype=ti.u8), arr_size: int, k: int,
                            hash_arr: ti.types.ndarray(dtype=ti.u64),
                            invalid_hash_arr: ti.types.ndarray(dtype=ti.u64),
                            missing_val_arr: ti.types.ndarray(dtype=ti.u8)):
    for i in range(arr_size):
        kmer2hash_taichi_uint64(arr, arr_size, i, k, hash_arr, invalid_hash_arr[0], missing_val_arr[0])


async def comp_kmer_hash_taichi(buffer: Buffer, kmer_len: int) -> Counter:
    """
    Compute kmer hash for each kmer from the input buffer, get the
    Args:
        buffer: a Buffer object that contain DNA sequences, A-0, C-1,
        kmer_len: length of kmer
    Returns: a numpy array
    """

    # await asyncio.sleep(np.random.rand()) # simulation of a time-consuming job

    missing_val = buffer.dtype(MISSING_VAL)
    missing_val_arr = np.array([missing_val])

    hash_dtype = get_hash_dtype(kmer_len)
    invalid_hash = get_invalid_hash(hash_dtype)
    invalid_hash_arr = np.array([invalid_hash])

    hash_arr = np.empty(buffer.buffer_size, dtype=hash_dtype)
    if hash_dtype == np.uint32:
        kmer2hash_kernel_uint32(buffer.buffer, buffer.pointer, kmer_len, hash_arr, invalid_hash_arr, missing_val_arr)
    elif hash_dtype == np.uint64:
        kmer2hash_kernel_uint64(buffer.buffer, buffer.pointer, kmer_len, hash_arr, invalid_hash_arr, missing_val_arr)
    else:
        raise Exception(f"Unknown kmer hash type hash_dtype={hash_dtype}")

    unique_hash, counts = np.unique(hash_arr[0:buffer.pointer], return_counts=True)
    inds = unique_hash != invalid_hash
    unique_hash = unique_hash[inds]
    counts = counts[inds]
    hash_counts_dict = Counter(dict(zip(unique_hash, counts)))

    return hash_counts_dict


# producer
async def chunk_reader(chunk_dir: str, task_queue: MaxSizeQueue):
    file_paths = get_chunk_file_paths(chunk_dir)
    for i_chunk, file_path in enumerate(file_paths):
        with open(file_path, "rb") as fh:
            buffer = pickle.load(fh)
            await task_queue.put(buffer)

    # Signal the consumer that no more items will be produced
    await task_queue.put(None)


def convert_input_chunks(fasta_file: str, buffer, out_dir="./chunks"):
    def write_chunk(i_chunk, chunk):
        filename = f"{out_dir}/chunk_{i_chunk}.pkl"
        with open(filename, "wb") as fh:
            pickle.dump(chunk, fh)

    out_dir = os.path.normpath(out_dir)
    out_dir = out_dir.rstrip(os.path.sep)
    mk_dir(out_dir)
    i_chunk = 0
    for arr in read_dnaseq_file(fasta_file):
        flag = buffer.append(arr)
        if not flag:
            data_to_write = buffer.data_to_write
            buffer.data_to_write = None
            buffer.id = i_chunk
            write_chunk(i_chunk, buffer)

            i_chunk += 1
            buffer.flush()
            buffer.append(data_to_write)

    if buffer.pointer > 0:
        buffer.id = i_chunk
        write_chunk(i_chunk, buffer)


# consumer of task queue
async def kmer_counter_chunk(task_queue: MaxSizeQueue, res_queue: MaxSizeQueue, kmer_len: int):
    while True:
        buffer = await task_queue.get()
        if buffer is None:
            await res_queue.put(None)
            print("All tasks have finished.")
            return None

        # process a task
        print(f"Processing chunk id={buffer.id}")
        # tmp_counter = await comp_kmer_hash(buffer, kmer_len)
        tmp_counter = await comp_kmer_hash_taichi(buffer, kmer_len)
        await res_queue.put(tmp_counter)


# consumer of result queue
async def sum_counting_res(res_queue: MaxSizeQueue):
    res_counter = Counter()
    while True:
        counter = await res_queue.get()
        # await asyncio.sleep(np.random.rand())  # simulation of a time-consuming job
        if counter is None:
            return res_counter
        else:
            res_counter += counter


# async def comp_kmer_hash(buffer: Buffer, kmer_len: int) -> Counter:
#     """
#     Compute kmer hash for each kmer from the input buffer, get the
#     Args:
#         buffer: a Buffer object that contain DNA sequences, A-0, C-1,
#         kmer_len: length of kmer
#     Returns: a numpy array
#     """
#
#     # await asyncio.sleep(np.random.rand()) # simulation of a time-consuming job
#
#     hash_dtype = get_hash_dtype(kmer_len)
#     # hash_mask = gen_hash_mask(kmer_len, hash_dtype)
#     invalid_hash = get_invalid_hash(hash_dtype)
#     hash_arr = np.empty(buffer.buffer_size, dtype=hash_dtype)
#
#     for i in range(buffer.pointer):
#         hash_arr[i] = kmer2hash(buffer.buffer, buffer.pointer, i, kmer_len, hash_arr, hash_dtype, invalid_hash)
#
#     unique_hash, counts = np.unique(hash_arr[0:buffer.pointer], return_counts=True)
#     inds = unique_hash != invalid_hash
#     unique_hash = unique_hash[inds]
#     counts = counts[inds]
#     hash_counts_dict = Counter(dict(zip(unique_hash, counts)))
#
#     return hash_counts_dict


async def count_kmer_producer_consumer_chunk(chunk_dir: str, kmer_len: int, q_size: int = 10):
    task_queue = MaxSizeQueue(maxsize=q_size)
    res_queue = MaxSizeQueue(maxsize=q_size)

    assert os.path.exists(chunk_dir) and len(get_chunk_file_paths(chunk_dir)) > 0

    producer_task = asyncio.create_task(chunk_reader(chunk_dir, task_queue))
    consumer_task = asyncio.create_task(kmer_counter_chunk(task_queue, res_queue, kmer_len))
    sum_task = asyncio.create_task(sum_counting_res(res_queue))

    await asyncio.gather(producer_task, consumer_task, sum_task)

    res_counter = sum_task.result()
    return res_counter


def proc_input(input_fasta_file: str, out_dir=".", buffer_size=2 ** 26):
    """
    process input fasta file, convert fasta file to chunks
    Args:
        input_fasta_file: input fasta file
        out_dir: output directory
        buffer_size: size of the chunk, int, 2**26 bytes are 64MB, 2**30 bytes are 1GB, should be less equal than 2**31

    Returns: write chunks as chunk_#.pkl pickle file under "chunks" directory
    """
    assert os.path.exists(out_dir)
    assert 0 < buffer_size <= 2 ** 31
    chunk_dir = os.path.join(out_dir, "chunks")

    if os.path.exists(chunk_dir):
        rm_files(chunk_dir)
    else:
        mk_dir(chunk_dir)

    buffer = Buffer(buffer_size)
    # convert input fasta file into pickle chunks
    convert_input_chunks(input_fasta_file, buffer, out_dir=chunk_dir)


def count_chunk_kmers(kmer_len, out_dir=".", q_size=20):
    """
    Count kmers in the chunks under "chunks" directory.
    The output is a Counter object (dictionary) and saved as .pkl file under the "kmer_counts"
    folder.
    Args:
        kmer_len: kmer len, int, should be 3-31
        q_size: queue size for concurrent processing of chunks, int, maximum number of chunks loaded into memory
        out_dir: output directory
        buffer_size: size of the chunk, int, 2**26 bytes are 64MB, 2**30 bytes are 1GB, should be less equal than 2**31

    Returns: a Counter object (dictionary), key is hash, value is count
    """
    assert os.path.exists(out_dir)
    chunk_dir = os.path.join(out_dir, "chunks")
    counts_dir = os.path.join(out_dir, "kmer_counts")

    assert os.path.exists(chunk_dir)
    mk_dir(counts_dir)

    # count kmers
    res_counter = asyncio.run(count_kmer_producer_consumer_chunk(chunk_dir, kmer_len, q_size))
    res_counter_file = f"k_{kmer_len}.pkl"
    with open(os.path.join(counts_dir, res_counter_file), "wb") as fh:
        pickle.dump(res_counter, fh)

    return res_counter


def count_kmer(input_fasta_file: str, kmer_len, q_size=20, out_dir=".", buffer_size=2 ** 26, rm_chunks_flag=True):
    """
    Count kmers of the input fasta file, which will be converted to chunks saved as .pkl file in folder "chunks" under
    the output directory. The output is a Counter object (dictionary) and saved as .pkl file under the "kmer_counts"
    folder.
    Args:
        input_fasta_file: path to input fasta file, str
        kmer_len: kmer len, int, should be 3-31
        q_size: queue size for concurrent processing of chunks, int, maximum number of chunks loaded into memory
        out_dir: output directory
        buffer_size: size of the chunk, int, 2**26 bytes are 64MB, 2**30 bytes are 1GB, should be less equal than 2**31
        rm_chunks_flag: remove all files under "chunks" folder, bool, True or False

    Returns: a Counter object (dictionary), key is hash, value is count
    """
    assert os.path.exists(out_dir)
    assert 0 < buffer_size <= 2 ** 31
    chunk_dir = os.path.join(out_dir, "chunks")
    counts_dir = os.path.join(out_dir, "kmer_counts")
    mk_dir(chunk_dir)
    mk_dir(counts_dir)

    if rm_chunks_flag:
        rm_files(chunk_dir)
        buffer = Buffer(buffer_size)
        # convert input fasta file into pickle chunks
        convert_input_chunks(input_fasta_file, buffer, out_dir=chunk_dir)

    # count kmers
    res_counter = asyncio.run(count_kmer_producer_consumer_chunk(chunk_dir, kmer_len, q_size))
    res_counter_file = f"k_{kmer_len}.pkl"
    with open(os.path.join(counts_dir, res_counter_file), "wb") as fh:
        pickle.dump(res_counter, fh)

    return res_counter


def merge_revcom(kmer_hash_counter: Counter, kmer_len: int, keep_lower_hash_flag=True):
    """
    merge reverse complements
    Args:
        kmer_hash_counter: Counter object, dictionary, key is kmer's hash, value is its count
        kmer_len: kmer length
        keep_lower_hash_flag: if keeping the lower hash as the key when merging a pair of reverse complements
    Returns:
        a counter object in which reverse complement counts are merged
    """
    uniq_kmer_hash_arr = np.array(list(kmer_hash_counter.keys()))
    revcom_uniq_kmer_hash_arr = get_revcom_hash_arr(uniq_kmer_hash_arr, kmer_len)

    if keep_lower_hash_flag:
        inds = np.where(uniq_kmer_hash_arr < revcom_uniq_kmer_hash_arr)[0]
    else:
        inds = np.where(uniq_kmer_hash_arr > revcom_uniq_kmer_hash_arr)[0]

    kh_arr = uniq_kmer_hash_arr[inds]
    cnt_arr = np.zeros_like(kh_arr, dtype=np.uint32)

    for i, kh in enumerate(kh_arr):
        orig_ind = inds[i]
        rc_kh = revcom_uniq_kmer_hash_arr[orig_ind]
        cnt_arr[i] = kmer_hash_counter[kh] + kmer_hash_counter.get(rc_kh, 0)

    palindrome_inds = np.where(uniq_kmer_hash_arr == revcom_uniq_kmer_hash_arr)[0]
    palindrome_kh_arr = uniq_kmer_hash_arr[palindrome_inds]
    palindrome_cnt_arr = np.zeros_like(palindrome_kh_arr, dtype=np.uint32)
    for i, kh in enumerate(palindrome_kh_arr):
        palindrome_cnt_arr[i] = kmer_hash_counter[kh]

    res_counter = Counter(dict(zip(chain(kh_arr, palindrome_kh_arr), chain(cnt_arr, palindrome_cnt_arr))))
    return res_counter


def cal_hamming_dist(kh_arr: np.ndarray, consensus_kh: np.uint64, kmer_len: int) -> np.ndarray:
    """
    calculate the Hamming distances between each element in kh_arr and the consensus sequence
    Args:
        kh_arr: kmer hash array
        consensus_kh: kmer hash of the consensus sequence
        kmer_len: kmer length
    Returns:
        Hamming distance array, np.ndarray object
    """
    ham_dist_arr = np.empty_like(kh_arr, dtype=np.uint8)
    hash_dtype = get_hash_dtype(kmer_len)
    consensus_kh_arr = np.array([consensus_kh], dtype=hash_dtype)
    hash_arr_size = len(kh_arr)

    if hash_dtype == np.uint32:
        cal_ham_dist_kernel_uint32(kh_arr, consensus_kh_arr, ham_dist_arr, hash_arr_size, kmer_len)
    elif hash_dtype == np.uint64:
        cal_ham_dist_kernel_uint64(kh_arr, consensus_kh_arr, ham_dist_arr, hash_arr_size, kmer_len)
    else:
        raise Exception(f"Unknown kmer hash type hash_dtype={hash_dtype}")
    return ham_dist_arr


def get_hamming_ball(kh_arr: np.ndarray, consensus_kh: np.uint64, kmer_len: int, max_ham_dist: int) -> np.ndarray:
    dist_arr = cal_hamming_dist(kh_arr, consensus_kh, kmer_len)
    hamming_ball_arr = kh_arr[dist_arr <= max_ham_dist]
    return hamming_ball_arr


def is_motif(kh_arr: np.ndarray, consensus_kh: np.uint64,
             kmer_len: int, max_ham_dist: int, revcom_flag=False) -> np.ndarray:
    """
    check if each kmer is a motif kmer given in kh_arr
    Args:
        kh_arr: input kmer hash array to be checked
        consensus_kh: consensus kmer hash
        kmer_len: kmer length
        max_ham_dist: maximum Hamming distance to the consensus, inclusive
        revcom_flag: if distance to the reverse complement of the consensus should be considered
    Returns:
        a logical np.ndarray
    """
    dist_arr = cal_hamming_dist(kh_arr, consensus_kh, kmer_len)

    if not revcom_flag:
        return dist_arr <= max_ham_dist

    rc_hash = revcom_hash(consensus_kh, kmer_len)
    rc_dist_arr = cal_hamming_dist(kh_arr, rc_hash, kmer_len)

    return np.logical_or(dist_arr <= max_ham_dist, rc_dist_arr <= max_ham_dist)


def contain_motif(kh_arr: np.ndarray, kh_len: int,
                  consensus_kh: np.uint64, consensus_kh_len: int, max_ham_dist: int,
                  revcom_flag=False):
    assert kh_len >= consensus_kh_len

    hash_dtype = get_hash_dtype(consensus_kh_len)
    motif_flag_arr = np.full_like(kh_arr, False, dtype=bool)
    mask = hash_dtype((1 << 2 * consensus_kh_len) - 1)

    for offset in range(kh_len - consensus_kh_len + 1):
        tmp_kh_arr = np.right_shift(kh_arr, 2 * offset)
        tmp_kh_arr = np.bitwise_and(tmp_kh_arr, mask)
        tmp_kh_arr = tmp_kh_arr.astype(hash_dtype)
        motif_flag_arr = np.logical_or(motif_flag_arr, is_motif(tmp_kh_arr, consensus_kh, consensus_kh_len,
                                                                max_ham_dist=max_ham_dist, revcom_flag=revcom_flag))
    return motif_flag_arr


def convert_kh_counter(kmer_hash_counter: Counter, kmer_len: int, target_kmer_len: int) -> Counter:
    assert kmer_len > target_kmer_len
    hash_dtype = get_hash_dtype(kmer_len)
    kh_arr = np.array(list(kmer_hash_counter.keys()), dtype=hash_dtype)
    kh_cnt_arr = np.array(list(kmer_hash_counter.values()), dtype=np.uint32)

    valid_inds = kh_cnt_arr > 0
    kh_arr = kh_arr[valid_inds]
    kh_cnt_arr = kh_cnt_arr[valid_inds]

    target_counter = Counter()
    target_hash_dtype = get_hash_dtype(target_kmer_len)
    mask = target_hash_dtype((1 << 2 * target_kmer_len) - 1)

    for i in range(kmer_len - target_kmer_len + 1):
        target_kh_arr = np.bitwise_and(np.right_shift(kh_arr, 2 * i), mask).astype(target_hash_dtype)
        target_counter += Counter(dict(zip(target_kh_arr, kh_cnt_arr)))

    # each target kmer (not the boundary ones) are counted (kmer_len - target_kmer_len + 1) times
    for k in target_counter:
        target_counter[k] //= (kmer_len - target_kmer_len + 1)

    return target_counter


def test_convert_kh_arr():
    seq = "TTTTCGTCCACGACGCTACCTTAAAGCATCCTTCTATGATACCATAGAAGCAGCTCCTTATCGTTTTAGCTTTCGTATTCGTCTAATCGTCTTTTACTCGACGAAAA"
    kmer_len = 8
    from inimotif import KmerCounter
    kc8 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kmer_dict = kc8.scan_seq(seq)
    c8 = Counter(kmer_dict)

    kmer_len = 11
    kc11 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kmer_dict11 = kc11.scan_seq(seq)
    c11 = Counter(kmer_dict11)

    c8_11 = convert_kh_counter(c11, 11, 8)
    print(c8_11 - c8) # should be empty
    print(c8 - c8_11)  # should give boundary elements
    for kh in (c8 - c8_11):
        print(kc8.hash2kmer(kh))

    print(c8_11 - merge_revcom(c8_11, 8))
    for kh in (c8_11 - merge_revcom(c8_11, 8)):
        print(kc8.hash2kmer(kh))


def test_motif():
    def calculate_hamming_dist(str1, str2):
        if len(str1) != len(str2):
            raise ValueError("Input strings must have equal length")

        hamming_dist = 0
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                hamming_dist += 1

        return hamming_dist

    seq = "TTTTCGTCCACGACGCTACCTTAAAGCATCCTTCTATGATACCATAGAAGCAGCTCCTTATCGTTTTAGCTTTCGTATTCGTCTAATCGTCTTTTACTC"
    kmer_len = 8
    from inimotif import KmerCounter
    kc11 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kmer_dict = kc11.scan_seq(seq)
    conseq = seq[0:kmer_len]
    conseq_kh = kc11.kmer2hash(conseq)
    n_all_kmer = len(kmer_dict)
    hash_dtype = get_hash_dtype(kmer_len)
    kh_arr = np.zeros(n_all_kmer, dtype=hash_dtype)
    for i, kh in enumerate(kmer_dict):
        kh_arr[i] = kh

    ham_ball_kh_arr = get_hamming_ball(kh_arr, conseq_kh, kmer_len, max_ham_dist=2)
    for kh in ham_ball_kh_arr:
        print(kc11.hash2kmer(kh))

    flag_arr = is_motif(ham_ball_kh_arr, conseq_kh,  kmer_len, max_ham_dist=2, revcom_flag=False)
    if np.all(flag_arr):
        print("all hamming ball kmers are motif")
    else:
        raise Exception("not all hamming ball kmers are motif")

    ham_ball_kh_arr[-1] = revcom_hash(ham_ball_kh_arr[-1], kmer_len)
    ham_ball_kh_arr[-2] = revcom_hash(ham_ball_kh_arr[-2], kmer_len)
    flag_arr = is_motif(ham_ball_kh_arr, conseq_kh, kmer_len, max_ham_dist=2, revcom_flag=False)
    assert not np.all(flag_arr)
    print("Last two kmer hash changed into reverse complement, now treated as non-motif.")
    flag_arr = is_motif(ham_ball_kh_arr, conseq_kh, kmer_len, max_ham_dist=2, revcom_flag=True)
    assert np.all(flag_arr)
    print("Last two kmer hash changed into reverse complement, treated as motif if revcom_flag=True.")

    # add one some non motif seq
    seqs = ["CATCCTTC", "GCAGCTCC"]
    nm_seqs_kh_arr = np.array([kc11.kmer2hash(seq) for seq in seqs], dtype=hash_dtype)
    flag_arr = is_motif(nm_seqs_kh_arr, conseq_kh, kmer_len, max_ham_dist=2, revcom_flag=False)
    if not np.any(flag_arr):
        print("all non_motif seqs are not motif seqs")
    else:
        raise Exception("some non motif sequences are treated as motif")

    # merge revcom with orig kmer hash
    merge_kmer_dict = merge_revcom(kmer_dict, kmer_len)
    for kh in merge_kmer_dict:
        assert kh <= kc11.revcom_hash(kh)


def test_contain_motif():
    seq = "TTTTCGTCCACGACGCTACCTTAAAGCATCCTTCTATGATACCATAGAAGCAGCTCCTTATCGTTTTAGCTTTCGTATTCGTCTAATCGTCTTTTACTC"
    kmer_len = 8
    from inimotif import KmerCounter
    kc8 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kmer_dict = kc8.scan_seq(seq)
    conseq = seq[0:kmer_len]
    conseq_kh = kc8.kmer2hash(conseq)
    n_all_kmer = len(kmer_dict)
    hash_dtype = get_hash_dtype(kmer_len)
    kh_arr = np.zeros(n_all_kmer, dtype=hash_dtype)
    for i, kh in enumerate(kmer_dict):
        kh_arr[i] = kh

    kmer_len = 10
    kc10 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kmer_dict10 = kc10.scan_seq(seq)
    kh_arr10 = np.zeros(n_all_kmer, dtype=hash_dtype)
    for i, kh in enumerate(kmer_dict10):
        kh_arr10[i] = kh

    flag_arr = contain_motif(kh_arr10, 10, conseq_kh, 8, 2, revcom_flag=False)
    for kh8, kh10 in zip(kh_arr[flag_arr], kh_arr10[flag_arr]):
        print(f"{kc8.hash2kmer(kh8)} {kc10.hash2kmer(kh10)} ")
    print()

    kh_arr8 = kh_arr[flag_arr]
    kh_arr10 = kh_arr10[flag_arr]
    kh_arr10[-2] = revcom_hash(kh_arr10[-2], kmer_len)
    kh_arr10[-1] = revcom_hash(kh_arr10[-1], kmer_len)
    flag_arr = contain_motif(kh_arr10, 10, conseq_kh, 8, 2, revcom_flag=True)
    assert np.all(flag_arr)
    for kh8, kh10 in zip(kh_arr8, kh_arr10):
        print(f"{kc8.hash2kmer(kh8)} {kc10.hash2kmer(kh10)} ")
    print()


def test_cal_hamming_dist():
    def calculate_hamming_dist(str1, str2):
        if len(str1) != len(str2):
            raise ValueError("Input strings must have equal length")

        hamming_dist = 0
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                hamming_dist += 1

        return hamming_dist

    seq = "TTTTCGTCCACGACGCTACCTTAAAGCATCCTTCTATGATACCATAGAAGCAGCTCCTTATCGTTTTAGCTTTCGTATTCGTCTAATCGTCTTTTACTC"
    kmer_len = 8
    from inimotif import KmerCounter
    kc11 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kmer_dict = kc11.scan_seq(seq)
    conseq = seq[0:kmer_len]
    conseq_kh = kc11.kmer2hash(conseq)
    n_all_kmer = len(kmer_dict)
    hash_dtype = get_hash_dtype(kmer_len)
    kh_arr = np.zeros(n_all_kmer, dtype=hash_dtype)
    for i, kh in enumerate(kmer_dict):
        kh_arr[i] = kh

    ham_dist_arr = cal_hamming_dist(kh_arr, conseq_kh, kmer_len)
    print(ham_dist_arr)

    for i, kh in enumerate(kmer_dict):
        assert ham_dist_arr[i] == calculate_hamming_dist(kc11.hash2kmer(kh), conseq)
    else:
        print(f"all calculated Hamming distances are correct.")

    hamming_ball_kh_arr = get_hamming_ball(kh_arr, conseq_kh, kmer_len, max_ham_dist=2)
    print(f"conseq={conseq}, hamming ball hashs are:")
    for kh in hamming_ball_kh_arr:
        print(f"{kc11.hash2kmer(kh)} d={calculate_hamming_dist(kc11.hash2kmer(kh), conseq)}")


def test_buffer():
    # Example usage
    buffer_size = 10
    buffer = Buffer(buffer_size)

    data1 = np.array([1, 2, 3], dtype=np.uint8)
    res1 = buffer.append(data1)

    data2 = np.array([4, 5, 6, 7, 8], dtype=np.uint8)
    res2 = buffer.append(data2)

    data3 = np.array([4, 5, 6, 7, 8], dtype=np.uint8)
    res3 = buffer.append(data3)

    print("Buffer contents:")
    print(buffer.buffer[:buffer.pointer])  # Print only the valid portion of the buffer
    print("Pointer position:", buffer.pointer)
    print(f"{res1= }")
    print(f"{res2= }")
    print(f"{res3= }")

    # this will raise an error
    # buffer.flush()

    # this is fine, need to handle data_to_write first
    buffer.data_to_write = None
    buffer.flush()


def test_gen_rand_fa_file(n_seq=100, min_len=30, max_len=60):
    import random

    def generate_random_sequence(length):
        bases = ['A', 'C', 'G', 'T']
        return ''.join(random.choice(bases) for _ in range(length))

    def generate_random_fasta_file(file_path, num_sequences):
        with open(file_path, 'w') as fasta_file:
            for i in range(num_sequences):
                sequence_length = random.randint(min_len, max_len)
                sequence = generate_random_sequence(sequence_length)
                fasta_file.write(f'>seq{i}\n')
                fasta_file.write(sequence + '\n')

    generate_random_fasta_file('random_sequences.fasta', n_seq)


def test_producer_consumer_chunk():
    input_fasta_file = "random_sequences.fasta"
    kmer_len = 8
    q_size = 20
    buffer_size = 100

    res_counter = count_kmer(input_fasta_file, kmer_len,
                             q_size=q_size, out_dir=".", buffer_size=buffer_size, rm_chunks_flag=False)
    print(res_counter)

    from inimotif import KmerCounter
    kc11 = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    kc11.scan_file(input_fasta_file)
    res1 = Counter(kc11.kmer_dict)

    if res_counter == res1:
        print("result are the same")
    else:
        print("result are not the same")
        print(res_counter - res1)
        invalid_hash = get_invalid_hash(get_hash_dtype(kmer_len))
        print(f"all kmers number in async: {sum(res_counter.values()) - res_counter[invalid_hash]}")
        print(f"all kmers number in ground truth: {sum(res1.values())}")

    uniq_kmer_hash_arr = np.array(list(res_counter.keys()))
    revcom_uniq_kmer_hash_arr = get_revcom_hash_arr(uniq_kmer_hash_arr, kmer_len)
    for kh, rc_kh in zip(uniq_kmer_hash_arr, revcom_uniq_kmer_hash_arr):
        assert kc11.revcom_hash(kh) == rc_kh
    else:
        print("all revcom hash are the same.")

    merged_counter = merge_revcom(res_counter, kmer_len)
    print(merged_counter-res_counter)


def preprocess(input_fasta_file: str, min_kmer_len, max_kmer_len, out_dir=".",
               q_size=20, buffer_size=2**26):
    assert min_kmer_len > 1
    assert max_kmer_len < 32 - 5
    # convert input fasta file into chunks
    proc_input(input_fasta_file, out_dir, buffer_size=buffer_size)

    # count kmers
    for kmer_len in range(min_kmer_len, max_kmer_len+5):
        count_chunk_kmers(kmer_len, out_dir, q_size=q_size)


def find_motif(out_dir: str, kmer_len: int, top_k: int, high_kmer_len: int, max_ham_dist: int, revcom_mode=True) -> List:
    # find consensus
    kmer_counts_file = os.path.join(out_dir, "kmer_counts", f"k_{kmer_len}.pkl")
    with open(kmer_counts_file, "rb") as fh:
        kh_counter = pickle.load(fh)

    high_kmer_counts_file = os.path.join(out_dir, "kmer_counts", f"k_{high_kmer_len}.pkl")
    with open(high_kmer_counts_file, "rb") as fh:
        high_kh_counter = pickle.load(fh)

    if revcom_mode:
        kh_counter = merge_revcom(kh_counter, kmer_len)  # merge reverse complements
        high_kh_counter = merge_revcom(high_kh_counter, high_kmer_len)

    res_conseq_kh_list = []
    high_kh_arr = np.array(list(high_kh_counter.keys()))

    # take out top_k consensuses
    for i in range(top_k):
        consensus_kh = kh_counter.most_common(1)[0][0]  # get consensus sequence

        # do something here, check if filtering criterion is met
        res_conseq_kh_list.append(consensus_kh)

        # remove this consenus motif from high kmer counter
        high_motif_flag_arr = contain_motif(high_kh_arr, high_kmer_len, consensus_kh, kmer_len,
                                       max_ham_dist=max_ham_dist, revcom_flag=revcom_mode)  # substring of high kmer might be revcom
        for kh in high_kh_arr[high_motif_flag_arr]:
            high_kh_counter[kh] = 0
        high_kh_arr = high_kh_arr[np.logical_not(high_motif_flag_arr)]

        # convert high_kh_counter to kmer counter, substring of high kmer might be revcom, need to merge
        kh_counter = convert_kh_counter(high_kh_counter, high_kmer_len, kmer_len)
        if revcom_mode:
            kh_counter = merge_revcom(kh_counter, high_kmer_len)

    return res_conseq_kh_list


def sample_code():
    # set parameters
    # test_gen_rand_fa_file(n_seq=100, min_len=30, max_len=60)
    input_fasta_file = "random_sequences.fasta"

    min_kmer_len = 7
    max_kmer_len = 21
    buffer_size = 100  # need to set to 2**26 in real case
    out_dir = "."

    # preprocess: convert input fasta file, count kmers
    preproc_flag = True
    if preproc_flag:
        preprocess(input_fasta_file, min_kmer_len, max_kmer_len, out_dir=out_dir, buffer_size=buffer_size)

    # find motifs
    kmer_len = 8
    max_ham_dist = 2
    top_k = 3
    revcom_mode = True
    high_kmer_len = 11
    assert min_kmer_len <= kmer_len <= max_kmer_len
    consensus_kh_list = find_motif(out_dir, kmer_len, top_k, high_kmer_len, max_ham_dist, revcom_mode)

    from inimotif import KmerCounter
    print(f"Top {top_k} consensus kmers")
    kc = KmerCounter(kmer_len, revcom_flag=False, unique_kmer_in_seq_mode=False)
    for kh in consensus_kh_list:
        print(kc.hash2kmer(kh))


if __name__ == "__main__":
    # test_gen_rand_fa_file()
    # test_buffer()
    # test_producer_consumer_chunk()
    # test_cal_hamming_dist()
    # test_motif()
    # test_contain_motif()
    # test_convert_kh_arr()
    sample_code()

import re

from inimotif_core import KmerCounter

# a pattern for matching DNA repeats
class RepeatPattern:
    def __init__(self, seq, n_min_rep, revcom_flag=True, n_max_rep=None):
        self.revcom_flag = revcom_flag
        
        if n_max_rep is None:
            n_max_rep = ''
        self.forward_pattern = re.compile(f'({seq}){{{n_min_rep},{n_max_rep}}}')
        
        if self.revcom_flag:
            kmer_counter = KmerCounter(len(seq))
            revcom_seq = kmer_counter.revcom(seq)
            self.revcom_pattern = re.compile(f'({revcom_seq}){{{n_min_rep},{n_max_rep}}}')
    
    def mask(self, in_str):
        def myrepl(m):
            return ('N'*len(m.group(0)))
        in_str = re.sub(self.forward_pattern, myrepl, in_str)
        if self.revcom_flag:
            in_str = re.sub(self.revcom_pattern, myrepl, in_str)
        return in_str

class Motif:
    def __init__(self, seq, n_max_mutation=0, revcom_flag=True):
        assert n_max_mutation<len(seq), f'n_max_mutation={n_max_mutation} is smaller than seq length {len(seq)}!'
        kc = KmerCounter(len(seq))
        self.kc = kc
        self.revcom_flag = revcom_flag
        
        forward_seq_hash = kc.kmer2hash(seq)
        self.forward_hamball = kc.get_hamming_ball(forward_seq_hash,n_max_mutation)
        
        if revcom_flag:
            revcom_seq_hash = kc.revcom_hash(forward_seq_hash)
            self.revcom_hamball = kc.get_hamming_ball(revcom_seq_hash,n_max_mutation)
    
    def mask(self, in_str):
        in_str_list = [c for c in in_str]
        len_str = len(in_str)
        k = self.kc.k

        i=0
        prev_hash = self.kc.dtype(-1)
        while(i<=len_str-k):
            tmpstr = in_str[i:i+k]
            # omit a kmer if it contains "N"
            if "N" in tmpstr:
                i += 1 + max([pos for pos,char in enumerate(tmpstr) if char == 'N'])
                prev_hash = self.kc.dtype(-1)
                continue
            if prev_hash==self.kc.dtype(-1):
                tmphash = self.kc.kmer2hash(tmpstr)
            # reuse hash in previous position
            else:
                tmphash = (prev_hash<< self.kc.dtype(2) ) & self.kc.mask
                tmphash += self.kc.base[ in_str[i+k-1] ]
            prev_hash = tmphash
            forward_flag = tmphash in self.forward_hamball
            revcom_flag = self.revcom_flag and (tmphash in self.revcom_hamball)
            if forward_flag or revcom_flag:
                for ii in range(i,i+k):
                    in_str_list[ii] = 'N' 
                i += k
                continue
            i += 1
        return "".join(in_str_list)

class Masker:
    def __init__(self):
        self.pattern_list = []
    
    def clear(self):
        self.pattern_list = []
        
    def add_rep_pat(self, seq, n_min_rep, revcom_flag):
        self.pattern_list.append( RepeatPattern(seq, n_min_rep, revcom_flag) )
    
    def add_motif_pat(self, seq, n_max_mutation, revcom_flag):
        self.pattern_list.append( Motif(seq, n_max_mutation, revcom_flag) )
        
    def mask(self,in_str):
        in_str = in_str.upper()
        for pat in self.pattern_list:
            in_str = pat.mask(in_str)
        return in_str
    
if __name__=="__main__":
    in_str = 'AAAAAAAAAAACGTGCCCCCGTGGGGGCGTGAAAACACGCCCCCACGTTTTCACGTTTTCACG'
    mask = Masker()
    mask.add_rep_pat('A',4,True)
    mask.add_rep_pat('C',4,True)
    mask.add_motif_pat('CGTG',1,False)
    # mask.add_motif_pat('CGTG',1,True)
    
    out_str = mask.mask(in_str)
    print(in_str)
    print(out_str)
    
    
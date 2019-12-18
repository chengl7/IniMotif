#!/usr/bin/env python3

"""
Author: Anna Price
Email: PriceA35@cardiff.ac.uk
"""

from Bio import SeqIO
import argparse
import csv
from itertools import groupby
import os
import glob

def genomeParser(chro, inputREF):
    # save the reference fasta to biopython Seq object
    
    # build path to the reference fasta
    fastaName = chro + ".fa*"
    findFasta = os.path.join(inputREF, fastaName)
    fastaPath = glob.glob(findFasta)
    
    # use biopython to parse the reference fasta file
    if fastaPath:
        for seq_rec in SeqIO.parse(str(fastaPath[0]), 'fasta'):
            genome = seq_rec.seq
    else:
        print ("Error: reference fasta not found")
        exit()

    return genome

def buildFASTA(chro, genome, seqName, startCoord, endCoord, outputFASTA):
    # build FASTA file with sequence information in the headers
    
    # write fasta file
    outfastaName = chro + "_windows" + ".fasta"
    with open(os.path.join(outputFASTA, outfastaName), 'w') as write_file:
        for name, coord1, coord2 in zip(seqName, startCoord, endCoord):
            FASTAhead = ">{0}|{1}|{2}".format(name, int(coord1), int(coord2))
            write_file.write(FASTAhead + '\n')
            seqWin = genome[int(coord1):int(coord2)]
            seqWin = str(seqWin)
            # convert sequence to uppercase
            seqWin = seqWin.upper()
            # chunk sequence into lengths of 50
            seq50 = [seqWin[i:i+50] for i in range(0, len(seqWin), 50)]
            for elem in seq50:
                write_file.write(elem + '\n')

def concatFASTA(outputFASTA):
    # concatenate the output fasta files into one file

    findFasta = os.path.join(outputFASTA, "*.fasta")
    collectFasta = glob.glob(findFasta)
    with open(os.path.join(outputFASTA, "concat.fasta"), 'w') as outfile:
        for fasta in collectFasta:
            with open(fasta) as infile:
                for line in infile:
                    outfile.write(line)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--input-tsv", dest="inputTSV", required=True, \
                        help="Path to the input chip-seq tsv")
    parser.add_argument("-r", "--input-ref", dest="inputREF", required=True, \
                        help="Path to the directory where the reference chromosomes are")
    parser.add_argument("-o", "--output-dir", dest="outputFASTA", required=True, \
                        help="Output directory path for the generated fasta files")
    parser.add_argument("-i", "--col-name", dest="colName", required=True, \
                        help="Column number for the chromosome name")
    parser.add_argument("-s", "--col-start", dest="colStart", required=True, \
                        help="Column number for the start coordinate of the window")
    parser.add_argument("-e", "--col-end", dest="colEnd", required=True, \
                        help="Column number for the end coordinate of the window")
    parser.add_argument("-c", "--concat-fasta", dest="concat", action="store_true", \
                        help="Concatenate the final fasta files into one file")
                        
    args = parser.parse_args()
    inputTSV = args.inputTSV
    inputREF = args.inputREF
    outputFASTA = args.outputFASTA
    colName = args.colName
    colStart = args.colStart
    colEnd = args.colEnd
    concat = args.concat
    
    # reset column numbers (counts from 0)
    colName = int(colName) - 1
    colStart = int(colStart) - 1
    colEnd = int(colEnd) - 1
    
    # initialise arrays for chromosome names and coordinates
    chrName = []
    startCoord = []
    endCoord = []
    
    # open tsv
    csv_reader = csv.reader(open(inputTSV), delimiter="\t")
    # skip header
    next(csv_reader)
    # split tsv by chromosome
    for chro, rows in groupby(csv_reader, lambda row: row[colName]):
        for row in rows:
            chrName.append(row[colName])
            startCoord.append(row[colStart])
            endCoord.append(row[colEnd])
        genome = genomeParser(chro, inputREF)
        buildFASTA(chro, genome, chrName, startCoord, endCoord, outputFASTA)
        # reset arrays for next chromosome
        chrName = []
        startCoord = []
        endCoord = []

    # if --concat-fasta flag used, then concatenate the output fasta files
    if concat:
        concatFASTA(outputFASTA)

if __name__ == "__main__":
    main()

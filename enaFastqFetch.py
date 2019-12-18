#!/usr/bin/env python3

"""
Author: Anna Price
Email: PriceA35@cardiff.ac.uk
"""

import argparse
from xmlparser import getXML, parseXMLgetFTP, parseFTPgetFASTQ
from reportwriter import writeReport

def main():
    parser = argparse.ArgumentParser()
    parser.set_defaults(method=getXML)
    parser.add_argument("-s", "--search-term", dest="search", required=True, \
                        help = "term you wish to search for, e.g. Mycobacterium, 1763, SRR5188398, SRX2504319, PRJNA360902, SELEX")
    parser.add_argument("-d", "--data-type", dest="dataType", required=True, \
                        help = "datatype you wish to search for, e.g. run, study, experiment, sample")
    parser.add_argument("-n", "--num-downloads", dest="numDown", \
                        help = "number of runs/studies/experiments/samples you wish to download")
    parser.add_argument("-r", "--report-file", dest="reportWrite", action="store_true", \
                        help = "generate a report file")
    args = parser.parse_args() 
    args.method(**vars(args))

    dataType = args.dataType
    numDown = args.numDown
    reportWrite = args.reportWrite
    
    parseXMLgetFTP('ena.xml', dataType, numDown)
    seqType = parseFTPgetFASTQ('fastq.txt')

    if reportWrite:
        writeReport('ena.xml', dataType, numDown, seqType)

if __name__ == "__main__":
    main()

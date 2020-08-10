#!/usr/bin/python
#-*- coding:utf-8 -*-
from __future__ import division
from tinyfasta import FastaParser
import argparse
import collections
import sys
import os
import math

#usage: python ~/github/NanoHub/src/NanoHub/SplitFasta.py --fasta M671-2999 --out temp.fasta --window 10000

def split_fasta(fa_file, out_file, window):
    window = int(window)

    if window == 0:
        print("Please make sure that window value is langer than 0.")
        sys.exit(1)

    out_h = open(out_file, "w")
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.split()[0].lstrip(">")
        seq = str(record.sequence)
        seqLen = len(seq)


        ### get the step
        step = math.ceil(seqLen / window)

        if step <= 2:
            newSeq = seq
            newDesc = "%s-%s" % (desc, "0")
            out_h.write(">%s\n%s\n" % (newDesc, newSeq))
        else:
            for i in range(step-2):
                newSeq = seq[i*window : (i+1)*window]
                newDesc = "%s-%s" % (desc, str(i))
                out_h.write(">%s\n%s\n" % (newDesc, newSeq))

            ### the last two records
            newSeq = seq[(i+1)*window :]
            newDesc = "%s-%s" % (desc, str(i+1))
            out_h.write(">%s\n%s\n" % (newDesc, newSeq))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Split the fasta file to small ones based on window length.")
    parser.add_argument("-f", "--fasta", help="The input fasta  file.")
    parser.add_argument("-w", "--window", default=10000, help="The window length.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    split_fasta(args.fasta, args.out, args.window)




if __name__ == "__main__":
    main()


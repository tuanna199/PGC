#!/usr/bin/python
from tinyfasta import FastaParser
import collections
import sys
import os
import argparse


#usage: python ~/github/NanoHub/src/NanoHub/FastaDiscardRetain.py --fasta ../Samples/11.contigs_cut.fa --name 11.contigs_blast_refGenome_contig_name.txt --out temp.fasta --method retain

def fasta_record(fa_file):
    FaSeq = {}
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.split()[0].lstrip(">")
        seq = str(record.sequence)
        FaSeq[desc] = seq
    return FaSeq


def query_seq_name(in_file):
    in_h = open(in_file, "r")
    Names = set()
    for line in in_h:
        lines = line.strip().split("\t")
        name = lines[0]
        Names.add(name)
    in_h.close()
    return Names


def discard_or_retain_record(fa_file, name_file, out_file, method):
    FaSeq = fasta_record(fa_file)
    SeqIDs = set(list(FaSeq.keys()))

    Names = query_seq_name(name_file)

    method = method.lower()
    if method == "discard":
        target = SeqIDs - Names
    elif method == "retain":
        target = Names
    else:
        print("Please make sure that the method is 'discard' or 'retain'.")
        sys.exit(1)

    nonTarget = Names - SeqIDs
    if len(nonTarget) != 0:
        print("Please check whether that extracted seq ids %s is in file %s." % (nonTarget, fa_file))
        sys.exit(1)

    out_h = open(out_file, "w")
    sortedNames = sorted(list(target))
    for t in target:
        if t in FaSeq:
            seq = FaSeq[t]
            out_h.write(">%s\n%s\n" % (t, seq))
        else:
            print("Please check whether that extracted seq id %s is in file %s." % (t, fa_file))
            sys.exit(1)
    out_h.close()




def main():
    parser = argparse.ArgumentParser(description="Retain or discard the seq based on seq id.")
    parser.add_argument("-f", "--fasta", help="The input fasta file.")
    parser.add_argument("-n", "--name", help="The seq id names.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-m", "--method", help="The method for seqlected seq id, 'retain' or 'discard'.")
    args = parser.parse_args()
    discard_or_retain_record(args.fasta, args.name, args.out, args.method)




if __name__ == "__main__":
    main()


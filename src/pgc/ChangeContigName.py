#!/usr/bin/python
from tinyfasta import FastaParser
import collections
import sys
import os
import argparse

#usage: python ~/github/NanoHub/src/NanoHub/ChangeContigName.py --fasta  /home/wuzhikun/PublicData/common_bean/Contigs/Samples/115.contigs_cut_univec_ref.fa  --out temp.fa  --tag tag_name --method sample


def change_contig_number(fa_file, out_file, tag_file):
    fileName = fa_file.split("/")[-1].split(".")[0]

    out_h = open(out_file, "w")
    tag_h = open(tag_file, "w")
    count = 1
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.split()[0].lstrip(">")
        seq = str(record.sequence)
        newDesc = "%s:contig%d" % (fileName, count)
        out_h.write(">%s\n%s\n" % (newDesc, seq))
        tag_h.write("%s\t%s\n" % (desc, newDesc))
        count += 1
    out_h.close()



def add_contig_name(fa_file, out_file, tag_file):
    """
    tag_file:
    ctg645-35_ctg645-38     M426-1:ctg645-35_ctg645-38
    ctg1119-1_ctg1119-3     M426-1:ctg1119-1_ctg1119-3
    ctg742-35_ctg742-37     M426-1:ctg742-35_ctg742-37
    """
    fileName = fa_file.split("/")[-1].split(".")[0].split("_")[0]

    out_h = open(out_file, "w")
    tag_h = open(tag_file, "w")
    count = 1
    for record in FastaParser(fa_file):
        desc = str(record.description)
        desc = desc.split()[0].lstrip(">")
        seq = str(record.sequence)
        newDesc = "%s:%s" % (fileName, desc)
        out_h.write(">%s\n%s\n" % (newDesc, seq))
        tag_h.write("%s\t%s\n" % (desc, newDesc))
        count += 1
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Change the contig name based on sample information.")
    parser.add_argument("-f", "--fasta", help="The input fasta file.")
    parser.add_argument("-t", "--tag", help="The out file for original and changed tag.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-m", "--method", help="'sample' means only add sample name to original tag, 'number' means add both sample name and contig number.")
    args = parser.parse_args()

    method = args.method.lower()
    if method == "number":
        change_contig_number(args.fasta, args.out, args.tag)
    elif method == "sample":
        add_contig_name(args.fasta, args.out, args.tag)
    else:
        print("Please make sure that the method should be 'number' or 'sample'.")
        sys.exit(1)


if __name__ == "__main__":
    main()







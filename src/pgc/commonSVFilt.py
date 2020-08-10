#!/usr/bin/python
from BaseFunc import parse_genotype_format
import argparse
import collections
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/commonSVFilt.py --vcf Sample_SV_common_CN329.vcf --out temp.vcf > temp

def common_SV_filt(SV_file, out_file):
    unvalueSet = {"./.", "0/0"}

    sv_h = open(SV_file, "r")
    out_h = open(out_file, "w")
    for line in sv_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Format = lines[8]
            Genos = lines[9:]
            genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genos]
            ResidueSet = set(genos) - unvalueSet
            if len(ResidueSet) == 0:
                print("The record is not valid because all genotypes are same to reference:\n%s\n" % line)
            else:
                out_h.write("%s\n" % line)




def main():
    parser = argparse.ArgumentParser(description="Filt the record if the all genotyepes are same to reference.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    common_SV_filt(args.vcf, args.out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
from BaseFunc import Infor_target_values
import argparse
import collections
import bisect

#usage: python ~/github/NanoHub/src/NanoHub/SVFiltTypeChr.py --vcf CN122.vcf --out temp.vcf --types DEL,INS,DUP,INV --chrs 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X > temp

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2020.02.19"


def vcf_filt_SV_type(vcf_file, out_file, types, chrs):
    TYPES = types.split(",")
    TYPES = [t.strip() for t in TYPES]

    CHRS = chrs.split(",")
    CHRS = [c.strip() for c in CHRS]

    out_h = open(out_file, "w")
    in_h = open(vcf_file, "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Chr = lines[0]
            start = lines[1]
            Alt = lines[4]
            Infor = lines[7]
            Type = Infor_target_values(Infor, "SVTYPE")

            if Type in TYPES and Chr in CHRS:
                out_h.write("%s\n" % line)
            else:
                print(line)
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Filt out the SV that existed in the given regions.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--types", default="DEL,INS,DUP,INV,TRA", help="The types of structural variations.")
    parser.add_argument("-c", "--chrs", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X", help="The target chromosomes.")
    args = parser.parse_args()
    vcf_filt_SV_type(args.vcf, args.out, args.types, args.chrs)

if __name__ == "__main__":
    main()


#!/usr/bin/python
from __future__ import division
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
import argparse
import sys
import math

### usage: python ~/github/NanoHub/src/NanoHub/VCFSVSeq.py --vcf M416-0.vcf --out M416-0.fasta

def SV_tag(lines):
    ID = lines[0]
    Chr, Start = lines[:2]
    Infor = lines[7]
    ### change the ID based on the Infor
    Chr2 = Infor_target_values(Infor, "CHR2")
    End = Infor_target_values(Infor, "END")
    SVType = Infor_target_values(Infor, "SVTYPE")
    if "SVLEN" in Infor:
        SVLength = Infor_target_values(Infor, "SVLEN")
    elif "AVGLEN" in Infor:
        SVLength = Infor_target_values(Infor, "AVGLEN")
        
    if SVLength.startswith("-"):
        SVLength = SVLength.lstrip("-")


    if SVType == "BND":
        SVType = "TRA"

    if SVType == "TRA":
        SVLength = 0

    Tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr2, End, SVLength, SVType) ### just one chromosomes
    
    return Tag


def vcf_SV_seq(vcf_file, out_file):
    vcf_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        else:
            lines = line.split("\t")
            Infor = lines[7]
            SVType = Infor_target_values(Infor, "SVTYPE")

            Tag = SV_tag(lines)
            
            if SVType == "INS" or SVType == "DEL":
                if "SEQ" in Infor:
                    Seq = Infor_target_values(Infor, "SEQ")
                    out_h.write(">%s\n%s\n" % (Tag, Seq))
                else:
                    print(Tag)
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the sequence based on vcf file from Sniffles.")
    parser.add_argument("-v", "--vcf", help="The input sv file with vcf format.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    vcf_SV_seq(args.vcf, args.out)




if __name__ == "__main__":
    main()
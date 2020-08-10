#!/usr/bin/python 
from BaseFunc import Infor_target_values
import argparse
import re

#usage: python ~/github/NanoHub/src/NanoHub/NanosvType.py --vcf M625-0.vcf --out M625_type.vcf

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.05.27"

def assign_SV_type(vcf_file, out_file):
    """
    inversion:
    N]1:23456789]  or [1:23456789[N

    deletion:
    N[1:23456789[

    duplications:
    ]1:23456789]N
    """
    vcf_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Chr1 = lines[0]
            Ref = lines[3]
            Alt = lines[4]
            Infor = lines[7]

            ### add END and BND for insertion
            if Alt.startswith("<") and Alt.endswith(">"):
                End = Infor_target_values(Infor, "END")

                BND = Ref + "[" + End + "["
                Infor = Infor + ";CHR2=" + Chr1 +";BND=" + BND
                out_h.write("%s\t%s\t%s\n" % ("\t".join(lines[:7]), Infor, "\t".join(lines[8:])))
            else:
            ### match the chromosome and position of BND
                string = "[\[\]](.+):(\d+)[\[\]]"
                ChrMatch = re.findall(string, Alt)
                if ChrMatch:
                    Chr2, End = ChrMatch[0]
                else:
                    print("Please check and make sure that there are chromosome and end position in the Alt in the record %s." % line)
                    
                if Chr1 == Chr2:
                    if Alt.startswith("[") or Alt.endswith("]"):
                        assignType = "<INV>"
                    elif Alt.endswith("["):
                        assignType = "<DEL>"
                    elif Alt.startswith("]"):
                        assignType = "<DUP>"
                else:
                    assignType = "<TRA>"

                Infor = Infor.replace("SVTYPE=BND", "SVTYPE=" + assignType.lstrip("<").rstrip(">"))
                # out_h.write("%s\t%s\t%s\t%s;CHR2=%s;END=%s;BND=%s\t%s\n" % ("\t".join(lines[:4]), assignType, "\t".join(lines[5:7]), Infor, Chr2, End, Alt, "\t".join(lines[8:])))
                out_h.write("%s\t%s\t%s\t%s;CHR2=%s;BND=%s\t%s\n" % ("\t".join(lines[:4]), assignType, "\t".join(lines[5:7]), Infor, Chr2, Alt, "\t".join(lines[8:])))
    vcf_h.close()
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Convert the BND to SV types for NanoSV results.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    assign_SV_type(args.vcf, args.out)

if __name__ == "__main__":
    main()


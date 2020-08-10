#!/usr/bin/python
from BaseFunc import Infor_target_values
import collections
import argparse
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/SVVCFExclude.py --bed /home/wuzhikun/Project/Population/mapping/minimap2/bed/CN008_SV_overlap.bed --vcf /home/wuzhikun/Project/Population/SVCall/FinalFilt/CN008_filt_centromere.vcf --out temp.vcf --method exclude

def bed_tags(bed_file):
    Tags = []
    bed_h = open(bed_file, "r")
    for line in bed_h:
        line = line.strip()
        if line != "":
            lines = line.split("\t")
            tag = lines[3]
            Tags.append(tag)
    bed_h.close()
    return Tags



def exclude_SV_record(bed_file, vcf_file, out_file, method):
    method = method.lower()

    Tags = bed_tags(bed_file)
    print(Tags)

    sv_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    if method == "exclude":
        for line in sv_h:
            line = line.strip()
            lines = line.split("\t")
            if line.startswith("#"):
                out_h.write("%s\n" % line)
            else:
                Chr, Start = lines[:2]
                Infor = lines[7]
                Chr2, End, SVlength, SVtype = Infor_target_values(Infor, "CHR2,END,SVLEN,SVTYPE")
                SVlength = SVlength.lstrip("-")
                tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr2, End, SVlength, SVtype)
                if tag not in Tags:
                    out_h.write("%s\n" % line)
    elif method == "include":
        for line in sv_h:
            line = line.strip()
            lines = line.split("\t")
            if line.startswith("#"):
                out_h.write("%s\n" % line)
            else:
                Chr, Start = lines[:2]
                Infor = lines[7]
                Chr2, End, SVlength, SVtype = Infor_target_values(Infor, "CHR2,END,SVLEN,SVTYPE")
                SVlength = SVlength.lstrip("-")
                tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr2, End, SVlength, SVtype)
                if tag in Tags:
                    out_h.write("%s\n" % line)
    else:
        print("Please make sure the parameter 'method' is 'inlcude' or 'exclude'.")
        sys.exit(1)
    sv_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Exclude the Tags with sample from the bed file.")
    parser.add_argument("-b", "--bed", help="The input bed file.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-e", "--method", help="The method, such as 'include' or 'exclude'.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    exclude_SV_record(args.bed, args.vcf, args.out, args.method)

if __name__ == "__main__":
    main()

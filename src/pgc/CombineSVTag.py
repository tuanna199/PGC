#!/usr/bin/python
from __future__ import division
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
from BaseFunc import SV_tag
import argparse
import sys
import math
import collections
import operator
from tinyfasta import FastaParser

#usage: python ~/github/NanoHub/src/NanoHub/CombineSVTag.py --vcf Sample_common_SV_convert.vcf --out temp.xls

# def Format_tag(Format, geno):
#     """
#     Format: GT:PSV:LN:DR:ST:TY:CO
#     geno: 1/1:NA:137:1,36:+-:INS,INS:1_136673-1_136891,1_136999-1_137209


#     1/1:NA:90:0,17:+-:INS:1_137053-1_137183;0/1:NA:188:29,19:+-:INS:1_136930-1_137112

#     ['NaN-0-NaN', '1_10169-X_449438-0-TRA']
#     """
#     LN = parse_genotype_format(Format, geno, "LN")
#     TY = parse_genotype_format(Format, geno, "TY")
#     CO = parse_genotype_format(Format, geno, "CO")
#     COs = CO.split(",")
#     TYs = TY.split(",")

#     Tags = []
#     COLen = len(CO)

#     ### set length of translocation as 0
#     if TY == "TRA":
#         LN = "0"

#     if COLen == 1:
#         tag = "%s-%s-%s" % (CO, LN, TY)
#         Tags.append(tag)
#     elif COLen > 1:
#         for c, t in zip(COs, TYs):
#             tag = "%s-%s-%s" % (c, LN, t)
#             Tags.append(tag)

#     TagLen = len(Tags)
#     if TagLen == 1:
#         return Tags[0]
#     else:
#         return ";".join(Tags)


def Format_tag(Format, Geno):
    """
    Format: GT:PSV:LN:DR:ST:TY:CO
    geno: 1/1:NA:137:1,36:+-:INS,INS:1_136673-1_136891,1_136999-1_137209

    geno:
    1/1:NA:90:0,17:+-:INS:1_137053-1_137183;0/1:NA:188:29,19:+-:INS:1_136930-1_137112

    ['NaN-0-NaN', '1_10169-X_449438-0-TRA']
    """
    genos = Geno.split(";")

    Tags = []
    for geno in genos:
        LN = parse_genotype_format(Format, geno, "LN")
        TY = parse_genotype_format(Format, geno, "TY")
        CO = parse_genotype_format(Format, geno, "CO")

        ### set length of translocation as 0
        if TY == "TRA":
            LN = "0"

        tag = "%s-%s-%s" % (CO, LN, TY)
        Tags.append(tag)

    TagLen = len(Tags)
    if TagLen == 1:
        return Tags[0]
    else:
        return ";".join(Tags)





def vcf_position_sample(vcf_file, out_file):
    TagSampleTag = collections.defaultdict(dict)
    SampleTags = collections.defaultdict(dict)

    out_h = open(out_file, "w")
    vcf_h = open(vcf_file, "r")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            lines = line.split("\t")
            samples = lines[9:]
            sampleNames = [s.split("/")[-1].split(".")[0].split("_")[0] for s in samples]
            out_h.write("Tag\t%s\n" % "\t".join(sampleNames))
        else:
            lines = line.split("\t")
            Infor = lines[7]
            Format = lines[8]
            Genos = lines[9:]
            SVType = Infor_target_values(Infor, "SVTYPE")

            # GenoTags = [Format_tag(Format, geno) for geno in Genos]
            GenoTags = [Format_tag(Format, geno) for geno in Genos]
            ### change "NaN-0-NaN" to "-"
            GenoTags = ["-" if g == "NaN-0-NaN" else g for g in GenoTags]

            Tag = SV_tag(lines)

            out_h.write("%s\t%s\n" % (Tag, "\t".join(GenoTags)))


            ### get denovo SV tag record


    out_h.close()




def main():
    parser = argparse.ArgumentParser(description="Filt SV file based on the number of support records.")
    parser.add_argument("-v", "--vcf", help="The input sv file with vcf format.")
    parser.add_argument("-o", "--out", help="The output file.")
    # parser.add_argument("-m", "--method", default="multiple", help="The file containing 'nultiple' of 'single' sample.")
    args = parser.parse_args()
    vcf_position_sample(args.vcf, args.out)





if __name__ == "__main__":
    main()

#!/usr/bin/python
import argparse
import collections
import sys
import os
import numpy

__author__ = "Zhikun Wu"
__date__ = "2019.11.05"
__email__ = "598466208@qq.com"

#usage: python ~/github/NanoHub/src/NanoHub/OverlapSVTypeLength.py --bed /home/wuzhikun/Project/NanoTrio/population/bed/overlap/Sample_DEL_overlap_filt.bed --SVLength /home/wuzhikun/Project/NanoTrio/population/bed/Sample_common_SV_length.xls --out temp.xls

def get_SV_record(bed_file, TypeChrRecord): #, length_file
    bed_h = open(bed_file, "r")
    for line in bed_h:
        lines = line.strip().split("\t")
        Chr, Start, End, Name = lines[:4]
        SVType = Name.split("-")[-1]
        Length = int(float(Name.split("-")[-2]))

        TypeChrRecord[SVType][Chr].append(Length)


        # ### get the length for SV tag
        # try:
        #     Length = SVNameLength[Name]
        #     ##################################################
        #     ##################################################
        #     ### INS should be less than 10 Mb
        #     if SVType == "INS":
        #         if Length < 10000000:
        #             TypeChrRecord[SVType][Chr].append(Length)
        #     else:
        #         TypeChrRecord[SVType][Chr].append(Length)
        # except KeyError:
        #     print("Please check whether the tag name %s is in tag record file %s." % (Name, length_file))
        #     sys.exit(1)
    bed_h.close()



# def SV_length(length_file):
#     """
#     length_file:
#     1_10468_3_198174528_TRA 100000.0
#     1_66288_66527_DEL       284.5
#     1_83940_84004_DEL       74.0
#     """
#     SVNameLength = {}

#     length_h = open(length_file, "r")
#     for line in length_h:
#         lines = line.strip().split("\t")
#         Name, Length = lines[:2]
#         Length = float(Length)
#         SVNameLength[Name] = Length
#     length_h.close()

#     return SVNameLength


def multiple_bed_record_length(bed_files, out_file):
    """
    bed_file:
    1   100922996   100923315   1_100922996_100923315_DEL   1   100922995   100923315   CHM1_chr1-100922996-DEL-320 DEL   320
    1   101253202   101253276   1_101253202_101253276_DEL   1   101253197   101253252   CHM1_chr1-101253198-DEL-55  DEL   55


    out_file:
    Chr     SVType  Number  TotalLength     AveLength
    1       DEL     1120    747112.4        667.1
    10      DEL     883     458618.7        519.4
    11      DEL     643     407550.9        633.8
    12      DEL     846     524027.2        619.4
    """
    # SVNameLength = SV_length(length_file)


    files = bed_files.split(",")
    files = [f.strip() for f in files]

    TypeChrRecord = collections.defaultdict(lambda: collections.defaultdict(list))
    for f in files:
        get_SV_record(f, TypeChrRecord) #SVNameLength, length_file


    ### Chr type record
    ChrTypeRecord = collections.defaultdict(dict)

    ### get all types and chromosomes
    Types = sorted(list(TypeChrRecord.keys()))
    Chrs = set()
    for t in Types:
        chrs = set(list(TypeChrRecord[t].keys()))
        for c in chrs:
            Chrs.add(c)

    ChrList = sorted(list(Chrs))

    for t in Types:
        for c in ChrList:
            if c in TypeChrRecord[t]:
                LenList = TypeChrRecord[t][c]
                LenNum = len(LenList)
                TotalLen = sum(LenList)
                AveLen = numpy.average(LenList)
            else:
                LenNum = 0
                TotalLen = 0
                AveLen = 0

            LenNum = str(LenNum)
            TotalLen = "%.1f" % TotalLen
            AveLen = "%.1f" % AveLen
            ChrTypeRecord[c][t] = [LenNum, TotalLen, AveLen]

    ### output the record
    out_h = open(out_file, "w")
    out_h.write("Chr\tSVType\tNumber\tTotalLength\tAveLength\n")
    for c in ChrList:
        for t in Types:
            record = ChrTypeRecord[c][t]
            out_h.write("%s\t%s\t%s\n" % (c, t, "\t".join(record)))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the SV type and number for each chromosome based on overlap bed file.")
    parser.add_argument("-b", "--bed", help="The input bed files which are separated with ','.")
    parser.add_argument("-o", "--out", help="The output file.")
    # parser.add_argument("-l", "--length", help="The input file contain SV name and length.")
    args = parser.parse_args()
    multiple_bed_record_length(args.bed, args.out)


if __name__ == "__main__":
    main()








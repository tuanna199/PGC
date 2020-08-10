#!/usr/bin/python
from BaseFunc import column_index
import collections
import argparse
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/SVListExclude.py --list Cell2019_Table_S1.txt --bed Cell2019_SV_dechr_DEL.bed --out Cell2019_SV_dechr_DEL_exclude_AKHX.bed --exclude "AK1,HX1"


def exclude_SV_record(SV_list, bed_file, out_file, exclude):
    """
    SV_list:

    #CHROM  POS     END     ID      SVTYPE  SVLEN   CONTIG_SUPPORT  CONTIG_DEPTH    CONTIG  CONTIG_NCBI     CONTIG_START    CONTIG_END      REPEAT_TYPE
    chr1    59598   59599   NA19434_chr1-59599-INS-308      INS     308     3       7       NA19434_chr1-20000-80000-ctg7180000000004       1-20000:0
    chr1    90068   90069   NA19240_chr1-90069-INS-59       INS     59      2       4       NA19240_chr1-80050-100552|ctg7180000000001|quiver/0_45755
    chr1    90178   90179   NA19434_chr1-90179-INS-58       INS     58      2       4       NA19434_chr1-60000-120000-ctg7180000000001      1-60000:0

    bed_file:
    1       136211  136308  NA19240_chr1-136212-DEL-97      DEL     97
    1       180096  180152  HG04217_chr1-180097-DEL-56      DEL     56
    1       180208  180608  HG00514_chr1-180209-DEL-400     DEL     400
    1       180235  180367  NA19434_chr1-180236-DEL-132     DEL     132
    """
    ### get the excluded IDS
    ExcludeIDS = []

    excludes = exclude.split(",")
    excludes = [e.strip() for e in excludes]

    sv_h = open(SV_list, "r")
    headers = sv_h.readline().strip().split("\t")
    IDIndex = column_index(headers, "ID")
    SampleIndex = column_index(headers, "MERGE_SAMPLES")

    for line in sv_h:
        lines = line.strip().split("\t")
        ID = lines[IDIndex]
        Sample = lines[SampleIndex]
        Samples = Sample.split(",")
        Samples = [s.strip() for s in Samples]

        sampleLen =len(Samples)
        if sampleLen == 1:
            if Samples[0] in excludes:
                ExcludeIDS.append(ID)
        elif sampleLen == 2:
            if set(Samples) == set(excludes):
                ExcludeIDS.append(ID)
    sv_h.close()

    out_h = open(out_file, "w")
    bed_h = open(bed_file, "r")
    for line in bed_h:
        line = line.strip()
        lines = line.split("\t")
        Tag = lines[3]
        if Tag in ExcludeIDS:
            print(Tag)
            continue
        else:
            out_h.write("%s\n" % line)
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Exclude the Tags with sample from the bed file.")
    parser.add_argument("-b", "--bed", help="The input bed file.")
    parser.add_argument("-s", "--list", help="The SV list file.")
    parser.add_argument("-e", "--exclude", default="AK1,HX1", help="The exclude sample name.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    exclude_SV_record(args.list, args.bed, args.out, args.exclude)

if __name__ == "__main__":
    main()

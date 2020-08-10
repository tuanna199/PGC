#!/usr/bin/python
import argparse
import sys
import os
import collections


#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/OverlapSVSummary.py --file /home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_DEL_overlap_filt.bed,/home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_INS_overlap_filt.bed,/home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_INV_overlap_filt.bed --out /home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_overlap_summary.txt


def SV_tag(SV_bed):
    TagSet = set()
    bed_h = open(SV_bed, "r")
    for line in bed_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        TagSet.add(tag)
    bed_h.close()
    return TagSet

def multiple_file_SV(fileStr, out_file):
    TypeList = ["DEL", "INS", "DUP", "INV"]
    FileTypeSV = collections.defaultdict(dict)
    files = fileStr.split(",")
    files = [f.strip() for f in files]
    for f in files:
        ff = f.split("/")
        fdir = ff[-2]
        fbase = ff[-1]
        for i in TypeList:
            if i in fbase:
                Type = i
                break

        TagSet = SV_tag(f)
        FileTypeSV[fdir][Type] = TagSet

    out_h = open(out_file, "w")
    out_h.write("Database\tType\tNumber\n")
    sortDirs = sorted(list(FileTypeSV.keys()))
    for d in sortDirs:
        for t in TypeList:
            if t in FileTypeSV[d]:
                tags = FileTypeSV[d][t]
                tagLen = len(tags)
                out_h.write("%s\t%s\t%d\n" % (d, t, tagLen))
            else:
                print("Please check whether there are SV type %s in database %s." % (t, d))

    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Statistics for SV type and number for multiple overlap bed files.")
    parser.add_argument("-f", "--file", help="The input bed files which contained tags.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    multiple_file_SV(args.file, args.out)

if __name__ == "__main__":
    main()



#!/usr/bin/env python
from __future__ import division
import collections
import argparse
import sys
import bisect

#usage: python ~/github/NanoHub/src/NanoHub/bedOverlapRecord.py --bed Sample_common_SV_DEL.bed --ratioThreshold 0.8 --out temp.txt

def region_overlap_length(Start1, End1, Start2, End2):
    region = [Start1, End1]
    startIndex = bisect.bisect(region, Start2)
    endIndex = bisect.bisect(region, End2)
    if startIndex == 1 and endIndex == 2:
        overlap = (Start2, End1)
        overlapLen = End1 - Start2
    elif startIndex == 0 and endIndex == 1:
        overlap = (Start1, End2)
        overlapLen = End2 - Start1
    elif startIndex == 0 and endIndex == 2:
        overlap = (Start1, End1)
        overlapLen = End1 - Start1
    elif startIndex == 1 and endIndex == 1:
        overlap = (Start2, End2)
        overlapLen = End2 - Start2
    else:
        print("Please check whether the region %s is overlap with %s." % ((Start1, End1), (Start2, End2)))
        sys.exit(1)

    return overlapLen


def bed_overlap_record(bed_file, out_file, ratioThreshold, method):
    """
    bed_file:
    1       100180721       100180847       1_100180721_100180847_DEL       1       100180635       100180944       AluYa5
    1       100181666       100181766       1_100181666_100181766_DEL       1       100181508       100181798       AluYd8
    1       100181666       100181766       1_100181666_100181766_DEL       1       100181763       100181858       SVA_E
    1       100922996       100923315       1_100922996_100923315_DEL       1       100923007       100923315       AluYe5
    1       101114416       101114581       1_101114416_101114581_DEL       1       101113749       101114646       L1PA6_3end
    """
    ratioThreshold = float(ratioThreshold)

    bed_h = open(bed_file, "r")
    out_h = open(out_file, "w")
    for line in bed_h:
        line = line.strip()
        lines = line.split("\t")
        Chr1, Start1, End1, Tag1 = lines[:4]
        Start1 = int(Start1)
        End1 = int(End1)
        Chr2, Start2, End2 = lines[4:7]
        Start2 = int(Start2)
        End2 = int(End2)

        Length1 = End1 - Start1
        Length2 = End2 - Start2

        overlapLen = region_overlap_length(Start1, End1, Start2, End2)

        lenRatio1 = overlapLen / Length1
        lenRatio2 = overlapLen / Length2

        # print(Tag1)
        # print(lenRatio1, lenRatio2)

        ### overlap strategy
        method = method.lower()

        if method == "first":
            if lenRatio1 >= ratioThreshold:
                out_h.write("%s\n" % line)
        elif method == "second":
            if lenRatio2 >= ratioThreshold:
                out_h.write("%s\n" % line)
        elif method == "both":
            if lenRatio1 >= ratioThreshold and lenRatio2 >= ratioThreshold:
                out_h.write("%s\n" % line)
        else:
            print("Please make sure that the 'method' should be 'first', 'second' or 'both'.")
            sys.exit(1)
    out_h.close()




def main():
    parser = argparse.ArgumentParser(description="Get the overlap region")
    parser.add_argument("-b", "--bed", help="The input intersect bed file derived from bedtools.")
    parser.add_argument("-r", "--ratioThreshold", default=0.5, help="The threshould of ratio of overlap.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-m", "--method", help="The selected one with length ratio compared to ratioThreshold, such as 'first', 'second', 'both'.")
    args = parser.parse_args()
    bed_overlap_record(args.bed, args.out, args.ratioThreshold, args.method)


if __name__ == "__main__":
    main()


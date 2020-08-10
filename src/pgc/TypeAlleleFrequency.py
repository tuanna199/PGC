#!/usr/bin/python
from __future__ import division
from BaseFunc import column_index
import argparse
import sys
import os
import bisect
import collections

#usage:  python ~/github/NanoHub/src/NanoHub/TypeAlleleFrequency.py --input Sample_SV_genotype_frequency.txt --out temp

def which_region(Regions, AltFreq):
    RegionIndex = []
    for region in Regions:
        AltIndex = bisect.bisect(region, AltFreq)
        RegionIndex.append(AltIndex)
    targetIndex = RegionIndex.index(1)
    targetRegion = Regions[targetIndex]
    return targetRegion


def allele_frequency_Stats(freq_file, out_file):
    """
    freq_file:
    Tag     Ref_freq        Alt_freq        MAF AltCount
    1_10059-1_10060-3242-INS        0.9975  0.0025  0.0025  1
    1_10380-1_10381-321-INS 0.9988  0.0012  0.0012  1
    1_10657-1_10760-78-INS  0.9988  0.0012  0.0012  1



    out_file:
    Type    0-0.05  0.05-0.1        0.1-0.2 0.2-0.3 0.3-0.4 0.4-0.5 0.5-0.6 0.6-0.7 0.7-0.8 0.8-0.9 0.9-1.0
    DEL     45658   2871    3432    2352    1943    1905    1062    285     121     68      38
    DUP     3142    179     144     70      23      24      8       1       0       1       0
    INS     36561   3149    3655    2456    2124    2317    1532    357     141     52      20
    INV     2080    53      40      28      13      17      8       6       4       4       6
    """
    Regions = [[0, 0.01],[0.01, 0.05], [0.05, 0.1], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.6], [0.6, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0]]

    TypeRegionCount = collections.defaultdict(lambda: collections.Counter())
    Types = set()

    in_h = open(freq_file, "r")
    headers = in_h.readline().strip().split("\t")
    RefIndex = column_index(headers, "Ref_freq")
    AltIndex = column_index(headers, "Alt_freq")
    TagIndex = column_index(headers, "Tag")
    # TypeIndex = column_index(headers, "Type")

    for line in in_h:
        lines = line.strip().split("\t")
        # Type = lines[TypeIndex]
        Tag = lines[TagIndex]
        Type = Tag.split("-")[-1]
        Types.add(Type)

        Alt = float(lines[AltIndex])


        if Alt == 1.0:
            targetRegion = [1.0]
        else:
            targetRegion = which_region(Regions, Alt)
        targetRegion = tuple(targetRegion)
        TypeRegionCount[Type][targetRegion] += 1
    in_h.close()

    out_h = open(out_file, "w")
    OutRegions = ["%s-%s" % (str(r[0]), str(r[1])) for r in Regions]
    out_h.write("Type\t%s\t1.0\n" % "\t".join(OutRegions))


    AllTypes = sorted(list(Types))
    for t in AllTypes:
        RegionCount = TypeRegionCount[t]

        Counts = []
        for r in Regions:
            r  = tuple(r)
            if r in RegionCount:
                count = RegionCount[r]
            else:
                count = 0
            Counts.append(count)


        Counts.append(RegionCount[(1.0,)])

        Counts = [str(c) for c in Counts]
        out_h.write("%s\t%s\n" % (t, "\t".join(Counts)))
    out_h.close()




def main():
    parser = argparse.ArgumentParser(description="Frequency for SV allele.")
    parser.add_argument("-i", "--input", help="The input file containing SV type and allele frequency.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    allele_frequency_Stats(args.input, args.out)



if __name__ == "__main__":
    main()

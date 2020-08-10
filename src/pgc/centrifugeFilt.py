#!/usr/bin/python
# from tinyfasta import FastaParser
from __future__ import division
from BaseFunc import column_index
import collections
import sys
import os
import argparse


#usage: python ~/github/NanoHub/src/NanoHub/centrifugeFilt.py --input NY8.contigs_cut_centrifuge.txt --out temp_tag --ratioThreshold 0.9

def filt_centrifuge_record(centrifuge_file, out_file, ratioThreshold):
    """
    centrifuge_file:
    $ head NY8.contigs_cut_centrifuge.txt
    readID  seqID   taxID   score   2ndBestScore    hitLength   queryLength numMatches
    k141_0  unclassified    0   0   0   0   371 1
    k141_2836   unclassified    0   0   0   0   314 1
    k141_2837   unclassified    0   0   0   0   436 1
    k141_1419   unclassified    0   0   0   0   332 1
    k141_2  unclassified    0   0   0   0   337 1
    k141_1420   unclassified    0   0   0   0   489 1
    k141_1  unclassified    0   0   0   0   540 1
    k141_1422   unclassified    0   0   0   0   309 1
    k141_1421   NC_000018.10    9606    64  0   23  373 1
    """
    ratioThreshold = float(ratioThreshold)

    cen_h = open(centrifuge_file, "r")
    headers = cen_h.readline().strip().split("\t")
    headers = [h.lower() for h in headers]
    tagIndex = headers.index("readid")
    hitLengthIndex = headers.index("hitlength")
    queryLengthIndex = headers.index("querylength")
    
    TagRatios = collections.defaultdict(list)


    for line in cen_h:
        line = line.strip()
        lines = line.split("\t")
        tag = lines[tagIndex]
        hitLength = lines[hitLengthIndex]
        queryLength = lines[queryLengthIndex]
        ratio = int(hitLength) / int(queryLength)
        TagRatios[tag].append(ratio)
    cen_h.close()

    ### output the record with ratio larger than ratio threshold
    out_h = open(out_file, "w")
    for tag in TagRatios:
        ratios = TagRatios[tag]
        maxRatio = max(ratios)
        if maxRatio > ratioThreshold:
            out_h.write("%s\n" % tag)
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Filt the record derived from centrifuge.")
    parser.add_argument("-i", "--input", help="The input file of centrifuge.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-r", "--ratioThreshold", default=0.9, help="The ratio threshold.")
    args = parser.parse_args()
    filt_centrifuge_record(args.input, args.out, args.ratioThreshold)



if __name__ == "__main__":
    main()




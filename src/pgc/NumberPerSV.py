#!/usr/bin/python
import collections
import argparse
import os 
import sys
import numpy

def number_per_SV(sv_tag, out_file, list_file):
    """
    sv_tag:
    Tag     CN001   CN002   CN003   CN004   CN005   CN006   CN007   CN008   CN009   CN010   CN011   CN012   CN013   CN014   CN015   CN016   CN017   CN018   CN01
    1_90312-1_90388-74.3-INS        -       -       1_90296-1_90378-58-INS  -       -       1_90378-1_90469-54-INS  1_90282-1_90367-109-INS 1_90371-1_90442-114-
    """
    SVTypeNumber = collections.defaultdict(list)

    sv_h = open(sv_tag, "r")
    headers = sv_h.readline().strip().split("\t")

    for line in sv_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        genos = lines[1:]
        SVType = tag.split("-")[-1]

        ### number per SV
        number = 0
        for g in genos:
            if g != "-":
                number += 1
        SVTypeNumber[SVType].append(number)
    sv_h.close()

    ### output the statistics
    out_h = open(out_file, "w")
    out_h.write("SVType\tSVNumber\tIndividuals\tStdev\n")

    list_h = open(list_file, "w")
    list_h.write("SVType\tNumber\n")

    types = sorted(list(SVTypeNumber.keys()))
    for t in types:
        numList = SVTypeNumber[t]
        listLen = len(numList)
        meanNum = numpy.mean(numList)
        std = numpy.std(numList)
        meanNum = "%.3f" % meanNum
        std = "%.3f" % std
        out_h.write("%s\t%d\t%s\t%s\n" % (t, listLen, meanNum, std))

        for i in numList:
            list_h.write("%s\t%d\n" % (t, i))
    out_h.close()
    list_h.close()



def main():
    parser = argparse.ArgumentParser(description="Statistics of individuals per SV for each type.")
    parser.add_argument("-t", "--tag", help="The input file with tag of population.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-l", "--list", help="The output file containing list of SV and number.")
    args = parser.parse_args()
    number_per_SV(args.tag, args.out, args.list)

if __name__ == "__main__":
    main()

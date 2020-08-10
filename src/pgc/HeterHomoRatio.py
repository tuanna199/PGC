#!/usr/bin/python
from __future__ import division
from BaseFunc import column_index
import collections
import argparse
import sys
import os
import numpy as np

#usage: python ~/github/NanoHub/src/NanoHub/HeterHomoRatio.py --genotype Sample_SV_genotype.txt --out temp

def SV_heter_homo_ratio(geno_file, out_file):
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    samples = headers[6:]

    chrs = set()

    ChrSampleHeter = collections.defaultdict(lambda: collections.Counter())
    ChrSampleHomo = collections.defaultdict(lambda: collections.Counter())
    for line in geno_h:
        lines = line.strip().split("\t")
        Chr = lines[0]
        chrs.add(Chr)

        genotypes = lines[6:]
        for i, g in enumerate(genotypes):
            if g == "0/1":
                ChrSampleHeter[Chr][i] += 1
            elif g == "1/1":
                ChrSampleHomo[Chr][i] += 1
    geno_h.close()

    # print(ChrSampleHeter)
    # print(ChrSampleHomo)


    out_h = open(out_file, "w")
    out_h.write("Chr\t%s\n" % "\t".join(samples))

    TotalValues = []
    Chrs = sorted(chrs)
    for C in CHRS:
        if C in Chrs:
            SampleRatios = []
            SampleValues = []
            for i, s in enumerate(samples):
                if i in ChrSampleHomo[C]:
                    SampleHomos = ChrSampleHomo[C][i]
                else:
                    SampleHomos = 0

                if i in ChrSampleHeter[C]:
                    SampleHeter = ChrSampleHeter[C][i]
                else:
                    SampleHeter = 0

                ### estimate ratio
                if SampleHomos == 0:
                    Ratio = "Inf"
                else:
                    Ratio = SampleHeter / SampleHomos
                    Ratio = "%.3f" % Ratio

                # print(SampleHeter, SampleHomos, Ratio)

                Value = "%d/%d" % (SampleHeter, SampleHomos)
                SampleRatios.append(Ratio)
                SampleValues.append(Value)
            out_h.write("%s\t%s\n" % (C, "\t".join(SampleRatios)))
            print("%s\t%s" % (C, "\t".join(SampleValues)))
            TotalValues.append(SampleValues)

    ### total ratio
    ValueTranspose = total_sample_values(TotalValues)
    SampleTotalRatio = []
    for sample in ValueTranspose:
        TotalHe = 0
        TotalHo = 0
        for v in sample:
            he, ho = v.split("/")
            he = int(he)
            ho = int(ho)
            TotalHe += he
            TotalHo += ho

        if TotalHo != 0:
            TotalRatio = TotalHe / TotalHo
            TotalRatio = "%.3f" % TotalRatio
        else:
            TotalRatio = "Inf"
        SampleTotalRatio.append(TotalRatio)
    out_h.write("Total\t%s\n" % "\t".join(SampleTotalRatio))
    out_h.close()


def total_sample_values(TotalValues):
    """
    TotalValues = [[1,2,3], [4,5,6], [7,8,9], [10,11,12]]

    ValueTranspose
    [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]
    """
    ValueTranspose = [[TotalValues[j][i] for j in range(len(TotalValues))] for i in range(len(TotalValues[0]))]
    return ValueTranspose

def main():
    parser = argparse.ArgumentParser(description="Estimate the heter to homo ratio for SV genotypes.")
    parser.add_argument("-g", "--genotype", help="The input genotype file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    SV_heter_homo_ratio(args.genotype, args.out)




if __name__ == "__main__":
    main()

#!/usr/bin/python
from __future__ import division
from BaseFunc import column_index
import collections
import argparse
import sys
import os
import numpy as np
import bisect


#usage: python ~/github/NanoHub/src/NanoHub/CategoryHeterHomo.py --genotype population/genotype/Sample_SV_genotype.txt --category type --out temp


def SV_heter_homo_ratio(geno_file, out_file, category):
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

    CategorySampleHeter = collections.defaultdict(lambda: collections.Counter())
    CategorySampleHomo = collections.defaultdict(lambda: collections.Counter())

    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    columns = headers[:6]
    columns = [c.lower() for c in columns]
    samples = headers[6:]

    ### target category
    targetCat = ""
    for c in columns:
        if category.lower() in c:
            targetCat = c
            CatIndex = columns.index(targetCat)
            break
    if targetCat == "":
        print("Please check whether that the target category %s is in file %s headers." % (category, geno_file))
        sys.exit(1)



    AllCats = set()
    for line in geno_h:
        lines = line.strip().split("\t")
        Chr = lines[0]
        # SVLength, SVType = lines[4:6]
        ### get target category
        Cat = lines[CatIndex]
        AllCats.add(Cat)
        ### select target Chromosomes
        if Chr in CHRS:
            genotypes = lines[6:]
            for i, g in enumerate(genotypes):
                if g == "0/1":
                    CategorySampleHeter[Cat][i] += 1
                elif g == "1/1":
                    CategorySampleHomo[Cat][i] += 1
    geno_h.close()


    out_h = open(out_file, "w")
    out_h.write("Category\t%s\n" % "\t".join(samples))


    ### length category
    LengthCatSum = collections.defaultdict(list)
    LengList = [100, 200, 500, 1000, 2000, 5000]
    Regions = ["50-100", "100-200", "200-500", "500-1000", "1000-2000", "2000-5000", ">5000"]
    LengListLen = len(LengList)

    TotalValues = []
    SortCats = sorted(AllCats)
    SortCatsLen = len(SortCats)
    for C in SortCats:
        SampleRatios = []
        SampleValues = []
        for i, s in enumerate(samples):
            if i in CategorySampleHomo[C]:
                SampleHomos = CategorySampleHomo[C][i]
            else:
                SampleHomos = 0

            if i in CategorySampleHeter[C]:
                SampleHeter = CategorySampleHeter[C][i]
            else:
                SampleHeter = 0

            ### estimate ratio
            if SampleHomos == 0:
                Ratio = "Inf"
            else:
                Ratio = SampleHeter / SampleHomos
                Ratio = "%.3f" % Ratio

            Value = "%d/%d" % (SampleHeter, SampleHomos)
            SampleRatios.append(Ratio)
            SampleValues.append(Value)
        ### total values
        TotalValues.append(SampleValues)


        if SortCatsLen < 10:
            out_h.write("%s\t%s\n" % (C, "\t".join(SampleRatios)))
            LengthCatSum = ""
        else:
            C = float(C)
            numIndex = bisect.bisect(LengList, C)
            region = Regions[numIndex]
            LengthCatSum[region].append(SampleValues)

    if LengthCatSum != "":
        for r in Regions:
            regionValues = LengthCatSum[r]
            regionTatio = accumulate_values(regionValues)
            # accuNum = len(regionValues)
            # regionTatio = [str(float(r) / accuNum) for r in regionTatio]
            out_h.write("%s\t%s\n" % (r, "\t".join(regionTatio)))


    ### total ratio
    SampleTotalRatio = accumulate_values(TotalValues)
    out_h.write("Total\t%s\n" % "\t".join(SampleTotalRatio))
    out_h.close()


def accumulate_values(TotalValues):
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
    return SampleTotalRatio 

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
    parser.add_argument("-c", "--category", help="The category which is in header of genotype file, such as chr, svlength, type.")
    args = parser.parse_args()
    SV_heter_homo_ratio(args.genotype, args.out, args.category)




if __name__ == "__main__":
    main()

#!/usr/bin/python
from __future__ import division
from BaseFunc import column_index
import collections
import argparse
import sys
import os
import numpy as np
import bisect


#usage: python ~/github/NanoHub/src/NanoHub/SVBurden.py --genotype  /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype_DUP.txt --out temp  --metafile /home/wuzhikun/Project/Population/meta_information.txt --category sex

#or:

#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/SVBurden.py --genotype /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype.txt --out /home/wuzhikun/Project/Population/population/genotype/Sample_SV_burden_age.txt --metafile /home/wuzhikun/Project/Population/meta_information.txt --category age

def category_samples(meta_file):
    """
    meta_file:
    Sample  Group   Type    Sex     Age
    M416-0  M416    Proband Unknown 0
    M416-1  M416    Parent  Male    0
    M425-0  M425    Proband Male    26
    """
    # SexSamples = collections.defaultdict(list)
    # AgeSamples = collections.defaultdict(list)
    SampleSex = {}
    SampleAge = {}
    meta_h = open(meta_file, "r")
    headers = meta_h.readline().strip().split("\t")
    sexIndex = column_index(headers, "sex")
    ageIndex = column_index(headers, "age")

    for line in meta_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        age = lines[ageIndex]
        sex = lines[sexIndex]

        sex = sex.lower()
        if sex == "male" or sex == "female":
            # SexSamples[sex].append(sample)
            SampleSex[sample] = sex

        age = int(float(age))
            # AgeSamples[age].append(sample)
        SampleAge[sample] = age
    meta_h.close()
    return SampleAge, SampleSex



def SV_heter_homo_ratio(geno_file, out_file, meta_file, category):
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

    SampleAge, SampleSex = category_samples(meta_file)



    SampleHeter = collections.Counter()
    SampleHomo = collections.Counter()

    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    columns = headers[:6]
    columns = [c.lower() for c in columns]
    samples = headers[6:]



    for line in geno_h:
        lines = line.strip().split("\t")
        Chr = lines[0]

        ### select target Chromosomes
        if Chr in CHRS:
            genotypes = lines[6:]
            for i, g in enumerate(genotypes):
                if g == "0/1":
                    SampleHeter[i] += 1
                elif g == "1/1":
                    SampleHomo[i] += 1
    geno_h.close()


    out_h = open(out_file, "w")
    out_h.write("Category\tSample\tHeter\tHomo\tTotal\n")


    # GroupSamples = collections.defaultdict(list)
    ### samples of category
    category = category.lower()
    if category == "sex":
        SampleGroup = SampleSex
    elif category == "age":
        SampleGroup = {}
        for s in SampleAge:
            a = SampleAge[s]
            if a == 0:
                continue
            elif a <= 30:
                SampleGroup[s] = "<30"
            elif a >= 40:
                SampleGroup[s] = ">40"
    else:
        print("Please select the category, such as sex of age")
        sys.exit(1)


    for i, s in enumerate(samples):
        ### get group
        if s in SampleGroup:
            group = SampleGroup[s]

            if i in SampleHeter:
                heterCount = SampleHeter[i]
            else:
                heterCount = 0

            if i in SampleHomo:
                homoCount = SampleHomo[i]
            else:
                homoCount = 0

            Total = heterCount + homoCount

            out_h.write("%s\t%s\t%d\t%d\t%d\n" % (group, s, heterCount, homoCount, Total))
    out_h.close()





def main():
    parser = argparse.ArgumentParser(description="Estimate the heter to homo ratio for SV genotypes.")
    parser.add_argument("-g", "--genotype", help="The input genotype file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-c", "--category", help="The category, such as sex, age.")
    parser.add_argument("-m", "--metafile", help="The meta information file.")
    args = parser.parse_args()
    SV_heter_homo_ratio(args.genotype, args.out, args.metafile, args.category)




if __name__ == "__main__":
    main()

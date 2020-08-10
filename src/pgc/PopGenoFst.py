#!/usr/bin/python
from __future__ import division
from BaseFunc import column_index
import collections
import argparse
import sys
import os
import numpy as np

#usage: python ~/github/NanoHub/src/NanoHub/PopGenoFst.py --structure /home/wuzhikun/Project/Population/population/admixture/Sample_assign_population.xls --genotype Sample_SV_genotype.test --out temp --target P1,P2

def population_group(pop_file):
    """
    geno    P1  P2  group
    CN001   0.852993    0.147007    P1
    CN002   0.921790    0.078210    P1
    CN003   0.913900    0.086100    P1
    """
    PopSamples = collections.defaultdict(list)
    SamplePop = {}

    pop_h = open(pop_file, "r")
    headers = pop_h.readline().strip().split("\t")
    try:
        genoIndex = column_index(headers, "geno")
        groupIndex = column_index(headers, "group")
    except ValueError:
        print("Please check whether the target columns exists in file %s." % pop_file)
        sys.exit(1)

    for line in pop_h:
        lines = line.strip().split("\t")
        if lines != []:
            geno = lines[genoIndex]
            group = lines[groupIndex]
            PopSamples[group].append(geno)
            SamplePop[geno] = group
    pop_h.close()
    return PopSamples, SamplePop


def geno_frequency(genos):
    genoLen = len(genos)
    pCount = 0
    qCount = 0
    valueLen = 0
    for g in genos:
        if g == "0/0":
            pCount += 2
            valueLen += 1
        elif g == "0/1":
            pCount += 1
            qCount += 1
            valueLen += 1
        elif g == "1/1":
            qCount += 2
            valueLen += 1
        else:
            continue

    if valueLen == 0:
        print("Please check whether all genotypes are missing.")
        return []
    else:
        pFreq = pCount / (valueLen * 2)
        qFreq = qCount / (valueLen * 2)
        Heter = pFreq * qFreq * 2
        return pFreq, qFreq, Heter




def population_geno_Fst(pop_file, geno_file, out_file, targetPop):
    PopSamples, SamplePop = population_group(pop_file)
    AllPop = sorted(list(PopSamples.keys()))



    ### parse the index of samples
    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    Samples = headers[6:]

    ### get the genotype index for each group 
    GroupIndex = collections.defaultdict(list)

    for i, s in enumerate(Samples):
        # try:
        #     pop = SamplePop[s]
        #     GroupIndex[pop].append(i)
        # except KeyError:
        #     print("Please check whether the sample %s is in population structure file %s." % (s, pop_file))
        #     sys.exit(1)

        if s in SamplePop:
            pop = SamplePop[s]
            GroupIndex[pop].append(i)

    GenoGroupList = sorted(list(GroupIndex.keys()))


    ### check whetehr the target population is population structure file.
    targetPops = targetPop.split(",")
    targetPops = [t.strip() for t in targetPops]
    for t in targetPops:
        if t not in AllPop:
            print("Please check whether the target population %s is in population structure file %s." % (targetPops, pop_file))
            sys.exit(1)

        if t not in GenoGroupList:
            print("Please check whether the target population %s is in genotype file %s." % (targetPops, geno_file))
            sys.exit(1)

    out_h = open(out_file, "w")
    for line in geno_h:
        line = line.strip()
        lines = line.split("\t")
        genotypes = lines[6:]

        GroupAlleleFreqs = []
        for t in targetPops:
            targetIndex = GroupIndex[t]

            targetGeno = [genotypes[i] for i in targetIndex]
            AlleleFreqs = geno_frequency(targetGeno)
            GroupAlleleFreqs.append(AlleleFreqs)

        ### estimate Fst
        if [] not in GroupAlleleFreqs:
            allHeter = [s[2] for s in GroupAlleleFreqs]
            allP = [s[0] for s in GroupAlleleFreqs]
            allQ = [s[1] for s in GroupAlleleFreqs]
            aveHeter = np.mean(allHeter)
            aveP = np.mean(allP)
            aveQ = np.mean(allQ)
            TotalHeter = aveP * aveQ * 2
            if TotalHeter != 0:
                Fst = (TotalHeter - aveHeter) / TotalHeter
                Fst = "%.3f" % Fst
                out_h.write("%s\t%s\n" % ("\t".join(lines[:6]), Fst))
        else:
            continue
    geno_h.close()
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Estimate Fst for sub-population.")
    parser.add_argument("-s", "--structure", help="The input sample structure file.")
    parser.add_argument("-g", "--genotype", help="The input genotype file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--target", help="The target sub-populations which are separated with ','.")
    args = parser.parse_args()
    population_geno_Fst(args.structure, args.genotype, args.out, args.target)

if __name__ == "__main__":
    main()






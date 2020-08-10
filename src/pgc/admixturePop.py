#!/usr/bin/python
import re
import argparse
import collections
import operator
import sys

#usage: python ~/github/NanoHub/src/NanoHub/admixturePop.py --cv /home/wuzhikun/Project/NanoTrio/population/admixture/admixture_cv.txt --ped /home/wuzhikun/Project/NanoTrio/population/plink/Sample_SV_geno.ped --population /home/wuzhikun/Project/NanoTrio/population/admixture/Sample_SV_geno.2.Q  --out temp_sample_pop.xls

def parse_admixture_CV(cv_file):
    """
    CV error (K=10): 0.22266
    CV error (K=1): 0.17613
    CV error (K=2): 0.18591
    CV error (K=3): 0.18467
    CV error (K=4): 0.19199
    CV error (K=5): 0.19697
    CV error (K=6): 0.20239
    CV error (K=7): 0.20794
    CV error (K=8): 0.21182
    CV error (K=9): 0.21926
    """
    PopCV = {}
    cv_h = open(cv_file, "r")
    for line in cv_h:
        line = line.strip()
        matchStr = re.findall("\(K=(\d+)\)", line)
        if matchStr:
            popNum = matchStr[0]
            ### population starts from at least two
            if popNum != "1":
                lines = line.split()
                cv = lines[-1]
                PopCV[popNum] = cv
    cv_h.close()

    ### select the population based on min cv value
    sortedCV = sorted(PopCV.items(), key=operator.itemgetter(1))
    Pop = sortedCV[0][0]
    print("The population should be devided to %s." % Pop)
    return Pop


def sample_assign_population(pop_file):
    """
    0.082050 0.917940 0.000010
    0.000010 0.999980 0.000010
    0.000010 0.999980 0.000010

    SamplePop = {1:1, 2:1}
    """
    IndexPop = {}
    IndexRecord = {}

    pop_h = open(pop_file, "r")
    sampleID = 1
    for line in pop_h:
        line = line.strip()
        values = line.split()
        values = [float(v) for v in values]
        ### index of max value
        maxValue = max(values)
        maxIndex = values.index(maxValue)
        IndexPop[sampleID] = "P" + str(maxIndex + 1)
        IndexRecord[sampleID] = line
        sampleID += 1
    pop_h.close()

    return IndexPop, IndexRecord

def ped_ID_sample(ped_file):
    IndexSample = {}
    ped_h = open(ped_file, "r")
    sampleInd = 1
    for line in ped_h:
        lines = line.strip().split("\t")
        ### get ind ID and name
        ind = lines[1]
        IndexSample[sampleInd] = ind
        sampleInd += 1
    ped_h.close()

    return IndexSample

def ped_sample_population(cv_file, pop_file, ped_file, out_file, number):
    Pop = parse_admixture_CV(cv_file)

    number = int(number)
    if number != 0:
        Pop = number
        Pop = str(Pop)
    print(Pop)

    ### index and sample
    IndexSample = ped_ID_sample(ped_file)
    ### get true name based on Pop
    ### Sample_SV_geno.3.P
    strings = pop_file.split(".")
    truePop = ".".join(strings[:-2]) + "." + Pop + "." + strings[-1]
    print(truePop)
    IndexPop, IndexRecord = sample_assign_population(truePop)

    out_h = open(out_file, "w")
    Pop = int(Pop)
    populations = ["P" + str(p+1) for p in range(Pop)]
    out_h.write("geno\t%s\tgroup\n" % "\t".join(populations))
    sortedIndex = sorted(list(IndexPop.keys()))
    for i in sortedIndex:
        pop = IndexPop[i]
        record = IndexRecord[i]
        try:
            sample = IndexSample[i]
        except KeyError:
            print("Please check whether the sample in line %d is in ped file %s." % (i, ped_file))
            sys.exit(1)
        new_record = "\t".join(record.split())
        out_h.write("%s\t%s\t%s\n" % (sample, new_record, pop))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Assign population for sample based on admixture.")
    parser.add_argument("-c", "--cv", help="The cv result from admixture.")
    parser.add_argument("-p", "--ped", help="The genotype file with ped format which containing sample name.")
    parser.add_argument("-q", "--population", help="The population file derived from admixture.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-n", "--number", default=0, help="The number of population, if not given, it detect by itself.")
    args = parser.parse_args()
    ped_sample_population(args.cv, args.population, args.ped, args.out, args.number)

if __name__ == "__main__":
    main()


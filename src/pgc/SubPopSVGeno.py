#!/usr/bin/python
from BaseFunc import Infor_target_values
from BaseFunc import Infor_substution_value
import collections
import os
import sys
import argparse


#usage: python ~/github/NanoHub/src/NanoHub/SubPopSVGeno.py --vcf population/Merge/Sample_SV_common.vcf --out temp.vcf --sample CN

def sub_population_SV_genotype(pop_geno, out_file, sample):
    samples = sample.split(",")

    pop_h = open(pop_geno, "r")
    out_h = open(out_file, "w")
    for line in pop_h:
        line = line.strip()
        if line.startswith("##"):
            out_h.write("%s\n" % line)
        elif line.startswith("#"):
            lines = line.split("\t")
            header = "\t".join(lines[:9])
            SAMPLES = lines[9:]
            ### get sample index
            targetSamples = []
            targetIndex = []
            for j in range(len(SAMPLES)):
                s = SAMPLES[j]
                for i in samples:
                    if i in s:
                        targetIndex.append(j)
                        targetSamples.append(s)
            out_h.write("%s\t%s\n" % (header, "\t".join(targetSamples)))
        elif line != "":
            lines = line.split("\t")
            value = "\t".join(lines[:9])
            GENOS = lines[9:]
            genotypes = [GENOS[i] for i in targetIndex]
            ### replace information
            Infor = lines[7]
            Supp = Infor_target_values(Infor, "SUPP_VEC")
            NewSupp = [Supp[i] for i in targetIndex]
            NewSuppRecord = "".join(NewSupp)
            NewInfor = Infor_substution_value(Infor, "SUPP_VEC", NewSuppRecord)

            genoCount = 0
            suppCount = 0
            for g in genotypes:
                # if g != "./.:NaN:0:0,0:--:NaN:NaN":
                gt = g.split(":")[0]
                if gt != "./.":
                    suppCount += 1
                    if gt != "0/0":
                        genoCount += 1
            NewInfor2 = Infor_substution_value(NewInfor, "SUPP", str(suppCount))
            if genoCount != 0:
                out_h.write("%s\t%s\t%s\t%s\n" % ("\t".join(lines[:7]), NewInfor2, lines[8], "\t".join(genotypes)))
    pop_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="The the count number of sliding window across the genome.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-s", "--sample", help="The sample name or partial name.")
    args = parser.parse_args()
    sub_population_SV_genotype(args.vcf, args.out, args.sample)

if __name__ == "__main__":
    main()

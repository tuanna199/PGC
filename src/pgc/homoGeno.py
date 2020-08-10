#!/usr/bin/python
import collections
import argparse
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/homoGeno.py --genotype Sample_SV_genotype.txt --tag temp_tag.xls --filt temp_file.txt

def extract_homo_geno(genos):
    homos = []
    genoLen = len(genos)
    for i in range(genoLen):
        geno = genos[i]
        if geno == "1/1":
            homos.append(i)
    return homos


def target_sample(samples, homoIndex):
    targetSample = []
    for i in homoIndex:
        s = samples[i]
        targetSample.append(s)
    return targetSample


def homo_genotype(geno_file, tag_file, filt_file):
    geno_h = open(geno_file, "r")
    header = geno_h.readline().strip()
    headers = header.split("\t")
    samples = headers[6:]

    tag_h = open(tag_file, "w")
    tag_h.write("Tag\tSamples\n")
    filt_h = open(filt_file, "w")
    filt_h.write("%s\n" % header)
    for line in geno_h:
        line = line.strip()
        lines = line.split("\t")
        tagInfor = lines[:6]
        genos = lines[6:]
        homoIndex = extract_homo_geno(genos)

        tag = "%s_%s-%s_%s-%s-%s" % tuple(tagInfor)
        if homoIndex != []:
            targetSample = target_sample(samples, homoIndex)
            filt_h.write("%s\n" % line)
            tag_h.write("%s\t%s\n" % (tag, ",".join(targetSample)))
    geno_h.close()
    filt_h.close()
    tag_h.close()


def main():
    parser = argparse.ArgumentParser(description='Fit the homo genotype SV record.')
    parser.add_argument('-g', '--genotype', help='The input genotype file.')
    parser.add_argument('-t', '--tag', help='The output tag and sample file.')
    parser.add_argument('-f', '--filt', help='The output filted file with homo genotypes.')
    args = parser.parse_args()
    homo_genotype(args.genotype, args.tag, args.filt)

    
if __name__ == '__main__':
    main()

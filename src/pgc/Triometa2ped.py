#!/usr/bin/python
#-*- coding-utf-8 -*-
from parseMetaInfor import MetaInfor
import argparse
import collections
import sys

#usage: python ~/github/NanoHub/src/NanoHub/meta2ped.py --meta /home/wuzhikun/Project/NanoTrio/meta_information.txt --genotype temp.txt --map test.map --ped test.ped

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.08.02"

def get_target_index(headers, target):
    """
    headers = ["Sample", "Group", "Type", "Sex", "Age"]
    """
    try:
        targetIndex = headers.index(target)
        return targetIndex
    except ValueError:
        print("please check whether the target %s is in the headers." % target)
        sys.exit(1)

def decode_sex(sex):
    if sex == "male":
        sexcode = "1"
    elif sex == "female":
        sexcode = "2"
    else:
        sexcode = "0"
    return sexcode

def meta_to_ped(meta_file):
    """
    meta_file:
    Sample  Group   Type    Sex Age
    M425-0  M425    Proband Male    26
    M425-1  M425    Parent  Male    60
    M425-2  M425    Parent  Female  51

    out_file:
    M425    M425-0  M425-1  M425-2  1   2
    M425    M425-1  0   0   1   1
    M425    M425-2  0   0   2   1

    SampleInfor (dict): Sample -> outfile information
    """
    SampleInfor = {}

    groupParents = collections.defaultdict(dict)
    sampleInfor = {}
    groupProband = {}

    in_h = open(meta_file, "r")
    headers = in_h.readline().strip().split("\t")
    headers = [h.lower() for h in headers]
    sampleIndex = get_target_index(headers, "sample")
    groupIndex = get_target_index(headers, "group")
    typeIndex = get_target_index(headers, "type")
    sexIndex = get_target_index(headers, "sex")

    for line in in_h:
        lines = line.strip().split("\t")
        Sample = lines[sampleIndex]
        Group = lines[groupIndex]
        Type = lines[typeIndex].lower()
        Sex = lines[sexIndex].lower()

        sampleInfor[Sample] = [Group, Type, Sex]


        if Type == "proband":
            groupProband[Group] = Sample
        elif Type == "parent":
            if Sex == "male":
                groupParents[Group]["father"] = Sample
            elif Sex == "female":
                groupParents[Group]["mother"] = Sample
            else:
                print("Please check the sex of sample in the record: %s" % line)
        else:
            print("Please check the affected status in the record: %s" % line)
            sys.exit(1)

    ### out to ped format
    # out_h = open(out_file, "w")
    for group in groupParents:
        try:
            proband = groupProband[group]
        except KeyError:
            print("Please check whether the group %s have the proband record in file %s." % (group, meta_file))
            sys.exit(1)

        sexProband = sampleInfor[proband][-1]
        sexProbandCode = decode_sex(sexProband)

        ### proband
        father = groupParents[group]["father"]
        mother = groupParents[group]["mother"]

        probandInfor = "%s\t%s\t%s\t%s\t%s\t2" % (group, proband, father, mother, sexProbandCode)
        SampleInfor[proband] = probandInfor

        ### parents
        sexFather = sampleInfor[father][-1]
        sexFatherCode = decode_sex(sexFather)
        sexMother = sampleInfor[mother][-1]
        sexMotherCode = decode_sex(sexMother)

        fatherInfor = "%s\t%s\t0\t0\t%s\t1" % (group, father, sexFatherCode)
        motherInfor = "%s\t%s\t0\t0\t%s\t1" % (group, mother, sexMotherCode)
        SampleInfor[father] = fatherInfor
        SampleInfor[mother] = motherInfor

    return SampleInfor, groupParents, groupProband


def decode_geno(geno):
    """
    geno: '0/0', '0/1', '1/1'
    genocode: "1 1", "1 2", "2 2", using number not 0, 0 means missing
    """
    if geno == "0/0":
        genocode = "1 1"
    elif geno == "0/1":
        genocode = "1 2"
    elif geno == "1/1":
        genocode = "2 2"
    else:
        print("Normal genotype is '0/0', '0/1', '1/1', please check the genotype %s." % geno)
        sys.exit(1)
    return genocode


def add_sample_geno(sampleNames, genotypes, SampleGenos, line):
    """
    sampleNames: ["M628-0",  "M628-1",  "M628-2"]
    genotypes: ["0/0", "1/1", "0/0"]
    sampleGenos (dict):  sample -> [geno1, geno2]
    """
    if len(sampleNames) == len(genotypes):
        for sample, geno in zip(sampleNames, genotypes):
            genocode = decode_geno(geno)
            SampleGenos[sample].append(genocode)
    else:
        print("Please check whether the numbers of samples and genotypes are identical for line %s." % line)
        sys.exit(1)


def genotype_to_ped(meta_file, geno_file, map_out, ped_out):
    """
    Chr1    Pos1    Chr2    Pos2    Type    M628-0  M628-1  M628-2
    1   10331   X   449317  BND 0/0 1/1 0/0
    1   10469   11  175280  BND 1/1 1/1 0/0
    1   10469   3   198174528   BND 1/1 0/0 0/0
    """
    SampleGenos = collections.defaultdict(list)

    map_h = open(map_out, "w")

    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    sampleNames = headers[5:]
    for line in geno_h:
        lines = line.strip().split("\t")
        
        ### deal with marker information
        Chr1, Pos1 = lines[:2]
        marker = "_".join(lines[:5])
        ### chromosome X codes as 23, Y codes as 24 
        if Chr1.upper() == 'X':
            Chr1 = "23"
        elif Chr1.upper() == "Y":
            Chr1 = "24"
        map_h.write("%s\t%s\t0\t%s\n" % (Chr1, marker, Pos1))

        ### deal with genotype information
        genotypes = lines[5:]
        add_sample_geno(sampleNames, genotypes, SampleGenos, line)
    geno_h.close()


    print(list(SampleGenos.keys()))
    ### deal with sample information and genotypes
    SampleInfor, groupParents, groupProband = meta_to_ped(meta_file)
    print(SampleInfor)
    print(groupParents)
    print(groupProband)

    ped_h = open(ped_out, "w")
    for group in groupProband:
        proband = groupProband[group]
        father = groupParents[group]["father"]
        mother = groupParents[group]["mother"]
        if proband in SampleGenos and father in SampleGenos  and mother in SampleGenos:
            proInfor = SampleInfor[proband]
            faInfor = SampleInfor[father]
            moInfor = SampleInfor[mother]
            proGeno = SampleGenos[proband]
            faGeno = SampleGenos[father]
            moGeno = SampleGenos[mother]
            ped_h.write("%s\t%s\n" % (proInfor, "\t".join(proGeno)))
            ped_h.write("%s\t%s\n" % (faInfor, "\t".join(faGeno)))
            ped_h.write("%s\t%s\n" % (moInfor, "\t".join(moGeno)))
    ped_h.close()


def main():
    parser = argparse.ArgumentParser(description="Convert the meta information file to ped format.")
    parser.add_argument("-i", "--meta", help="The input meta information file.")
    parser.add_argument("-g", "--genotype", help="The input genotype file.")
    parser.add_argument("-m","--map", help="The output map file.")
    parser.add_argument("-p", "--ped", help="The output ped file.")
    args = parser.parse_args()
    genotype_to_ped(args.meta, args.genotype, args.map, args.ped)

if __name__ == "__main__":
    main()

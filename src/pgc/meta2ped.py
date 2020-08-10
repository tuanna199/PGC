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
    """
    1=male; 2=female; other=unknown
    """
    sex = sex.lower()
    if sex == "male":
        sexcode = "1"
    elif sex == "female":
        sexcode = "2"
    elif sex == "unknown":
        sexcode = "0"
    else:
        sexcode = "0"
        print("Please check the sex status: %s, it should be 'male', 'female' or 'unknown'." % sex)
    return sexcode

def decode_affect(affect):
    """
    -9 missing 
     0 missing
     1 unaffected
     2 affected
    """
    affect = affect.lower()
    if affect in ["proband", "case", "unnormal"]:
        affectCode = "2"
    elif affect in ["parent", "control", "normal"]:
        affectCode = "1"
    elif affect == "unknown":
        affectCode = "0"
    else:
        affectCode = "0"
        print("The affect status name %s is not in the given list, please chech!" % affect)
    return affectCode




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
        Affect = lines[typeIndex]
        Sex = lines[sexIndex]

        ### get the sex and affect code
        SexCode = decode_sex(Sex)
        AffectCode = decode_affect(Affect)

        ### get the sample information
        sampleInfor[Sample] = "\t".join([Sample, Sample, "0", "0", SexCode, AffectCode])
    in_h.close()

    return sampleInfor




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

    map_out:
    1   1_10002_3_10918_TRA 0   10002
    1   1_10031_3_10854_TRA 0   10031
    1   1_10022_11_176424_TRA   0   10022
    1   1_10002_11_176543_TRA   0   10002
    1   1_10098_X_449431_TRA    0   10098

    ped_out:
    M425-1  M425-1  0       0       1       1       1 1     1 1     2 2     1 1     1 1     1 1     1 1     1 1     1 1     1 1
    M426-1  M426-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 2     1 1     1 1     1 1     1 1
    M446-1  M446-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1
    M452-1  M452-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     2 2     2 2
    M462-1  M462-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1
    M470-1  M470-1  0       0       1       1       1 1     1 1     1 1     1 1     2 2     1 1     1 1     1 1     1 1     1 1
    M473-1  M473-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 1     1 2     1 1     1 1     1 1
    M489-1  M489-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1     1 1
    M502-1  M502-1  0       0       1       1       1 1     1 1     1 1     1 1     1 1     1 1     1 1     2 2     1 1     1 1
    """
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"] #, "X", "Y"

    SampleGenos = collections.defaultdict(list)

    map_h = open(map_out, "w")

    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    sampleNames = headers[6:]
    for line in geno_h:
        lines = line.strip().split("\t")
        
        ### deal with marker information
        Chr1, Pos1 = lines[:2]
        marker = "%s_%s-%s_%s-%s-%s" % tuple(lines[:6])
        
        ### select the chromosome 
        if Chr1 in CHRS:
            # ### chromosome X codes as 23, Y codes as 24 
            # if Chr1.upper() == 'X':
            #     Chr1 = "23"
            # elif Chr1.upper() == "Y":
            #     Chr1 = "24"
                
            map_h.write("%s\t%s\t0\t%s\n" % (Chr1, marker, Pos1))

            ### deal with genotype information
            genotypes = lines[6:]
            add_sample_geno(sampleNames, genotypes, SampleGenos, line)
    geno_h.close()


    ### deal with sample information and genotypes
    SampleInfor = meta_to_ped(meta_file)


    ### out to ped file
    ped_h = open(ped_out, "w")

    for sample in SampleGenos:
        if sample in SampleInfor:
            Infor = SampleInfor[sample]
            genos = SampleGenos[sample]
            ped_h.write("%s\t%s\n" % (Infor, "\t".join(genos)))
        else:
            print("Please check wether the sample information of %s is in the mata file %s." % (sample, meta_file))
            sys.exit(1)

    map_h.close()
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

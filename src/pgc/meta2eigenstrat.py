#!/usr/bin/python
from __future__ import division
import argparse
import numpy
import sys

#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/meta2eigenstrat.py --meta /home/wuzhikun/Project/NanoTrio/meta_information.txt --genotype /home/wuzhikun/Project/NanoTrio/population/Sample_SV_genotype.txt --ind /home/wuzhikun/Project/NanoTrio/population/eigenstrat/Sample_SV.ind --geno /home/wuzhikun/Project/NanoTrio/population/eigenstrat/Sample_SV.geno --snp /home/wuzhikun/Project/NanoTrio/population/eigenstrat/Sample_SV.snp

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.08.05"

def decode_geno(geno):
    """
    geno: '0/0', '0/1', '1/1'
    genocode: "1 1", "1 2", "2 2", using number not 0, 0 means missing
    """
    if geno == "0/0":
        genocode = "0"
    elif geno == "0/1":
        genocode = "1"
    elif geno == "1/1":
        genocode = "2"
    else:
        print("Normal genotype is '0/0', '0/1', '1/1', please check the genotype %s." % geno)
        sys.exit(1)
    return genocode

def random_sex(sexList):
    """
    sexList = ["F", "M"]
    """
    sexIndex = numpy.random.randint(2, size=1)
    sexcode = sexList[sexIndex[0]]
    return sexcode

def decode_sex(sex):
    """
    1=male; 2=female; other=unknown
    """
    sexList = ["F", "M"]

    sex = sex.lower()
    if sex == "male":
        sexcode = "M"
    elif sex == "female":
        sexcode = "F"
    ### random sex if tag as unknown
    elif sex == "unknown":
        sexcode = random_sex(sexList)
    else:
        sexcode = random_sex(sexList)
        print("Please check the sex status: %s, it should be 'male', 'female' or 'unknown'." % sex)
    return sexcode

def decode_affect(affect):
    """
    -9 missing 
     0 missing
     1 unaffected
     2 affected
    """
    affectList = ["Case", "Control"]

    affect = affect.lower()
    if affect in ["proband", "case", "unnormal"]:
        affectCode = "Case"
    elif affect in ["parent", "control", "normal"]:
        affectCode = "Control"
    ### random affect if tag as unknown
    elif affect == "unknown":
        affectCode = random_sex(affectList)
    else:
        affectCode = random_sex(affectList)
        print("The affect status name %s is not in the given list, please chech!" % affect)
    return affectCode


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

def ramdom_two_diff_base():
    baseList = ["A", "T", "C", "G"]
    baseIndex = numpy.random.randint(4, size=2)
    bases = [baseList[b] for b in baseIndex]
    return bases

def meta_geno_to_eigenstrat(meta_file, geno_file, snp_out, geno_out, ind_out):
    """
    meta_file:
    Sample  Group   Type    Sex Age
    M416-0  M416    Proband Unknown 0
    M416-1  M416    Parent  Male    0
    M416-2  M416    Parent  Female  0

    geno_file:
    Chr1    Pos1    Chr2    Pos2    Type    M425-1  M426-1  M446-1  M452-1  M462-1  M470-1  M473-1  M489-1  M502-1  M509-1  M530-
    1       10002   3       10918   TRA     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     1/1
    1       10031   3       10854   TRA     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       10022   11      176424  TRA     1/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0


    geno_out: (one line per sample)
    11100
    01212
    21101

    ind_out:
    SAMPLE0 F       Case
    SAMPLE1 M       Case
    SAMPLE2 F    Control

    snp_out:
    rs0000  11        0.000000               0 A C
    rs1111  11        0.001000          100000 A G
    rs2222  11        0.002000          200000 A T
    """
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"] #, "X", "Y"

    geno_h = open(geno_file, "r")
    snp_h = open(snp_out, "w")
    out_h = open(geno_out, "w")

    headers = geno_h.readline().strip().split("\t")
    SAMPLES = headers[6:]
    for line in geno_h:
        lines = line.strip().split("\t")
        
        ### deal with marker information
        Chr1, Pos1 = lines[:2]
        marker = "_".join(lines[:6])

        ### select the chromosome 
        if Chr1 in CHRS:
            genePos = int(Pos1) / 10000000
            genePos = "%.6f" % genePos

            randombase = ramdom_two_diff_base()
            snp_h.write("%s\t%s\t%s\t%s\t%s\n" % (marker, Chr1, genePos, Pos1, "\t".join(randombase)))

            ### output genotypes
            genos = lines[6:]
            genotypes = [decode_geno(g) for g in genos]
            out_h.write("%s\n" % "".join(genotypes))
    snp_h.close()
    geno_h.close()
    out_h.close()


    ### deal with meta information
    SampelInfors = {}

    meta_h = open(meta_file, "r")
    headers = meta_h.readline().strip().split("\t")
    headers = [h.lower() for h in headers]
    sampleIndex = get_target_index(headers, "sample")
    groupIndex = get_target_index(headers, "group")
    typeIndex = get_target_index(headers, "type")
    sexIndex = get_target_index(headers, "sex")

    for line in meta_h:
        lines = line.strip().split("\t")
        Sample = lines[sampleIndex]
        Group = lines[groupIndex]
        Affect = lines[typeIndex]
        Sex = lines[sexIndex]

        ### get the sex and affect code
        SexCode = decode_sex(Sex)
        AffectCode = decode_affect(Affect)
        infor = "%s\t%s\t%s" % (Sample, SexCode, AffectCode)
        SampelInfors[Sample] = infor
    meta_h.close()

    ### output meta information
    ind_h = open(ind_out, "w")
    for s in SAMPLES:
        if s in SampelInfors:
            infor = SampelInfors[s]
        else:
            print("Please check whether the sample name %s in genotype file %s is in the meta file %s." % (s, geno_file, meta_file))
        ind_h.write("%s\n" % infor)
    ind_h.close()

def main():
    parser = argparse.ArgumentParser(description="Convert the meta information file to ped format.")
    parser.add_argument("-i", "--meta", help="The input meta information file.")
    parser.add_argument("-g", "--genotype", help="The input genotype file.")
    parser.add_argument("-d","--ind", help="The output ind file.")
    parser.add_argument("-s", "--snp", help="The output snp file.")
    parser.add_argument("-o", "--geno", help="The output geno file.")
    args = parser.parse_args()
    meta_geno_to_eigenstrat(args.meta, args.genotype, args.snp, args.geno, args.ind)


if __name__ == "__main__":
    main()








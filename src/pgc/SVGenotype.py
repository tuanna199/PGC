#!/usr/bin/python
#-*- coding:utf-8 -*-
from __future__ import division
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
from BaseFunc import parse_genotype_format_mul
import argparse
import collections
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/SVGenotype.py --vcf M628_common.vcf --out temp.txt --frequency frequency.txt

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.08.02"

def parse_sv_genotype(Format, GenoInfor):
    """
    Format:
    GT:PSV:LN:DR:ST:TY:CO

    GenoInfor:
    ["./.:NaN:0:0,0:--:NaN:NaN", "1/1:NA:438986:0,2:+-:TRA:1_10331-X_449317", "./.:NaN:0:0,0:--:NaN:NaN"]
    """
    ############################## just get first genotype ##########################
    ### genotype
    genotypes = [parse_genotype_format_mul(Format, genoinfor, "GT") for genoinfor in GenoInfor]
    genotypes = [g if g != "./." else "0/0" for g in genotypes]

    ### reference and alternative count
    AltCounts = [parse_genotype_format_mul(Format, genoinfor, "DR") for genoinfor in GenoInfor]
    Counts = [[int(cc) for cc in c.split(",")] for c in AltCounts]
    CountsSum = [sum(a) for a in Counts]
    ### alterbative frequency
    AltFreq = []
    CountFrequency = []
    for csum, count in zip(CountsSum, Counts):
        if csum == 0:
            freq = 0
        else:
            freq = count[1] / csum
        AltFreq.append(freq)

        freq = "%.3f" % freq

        if freq == "0.000":
            freq = "0"
        if freq == "1.000":
            freq = "1"
        
        contfreq = ",".join(map(str, count)) + ":" + freq
        CountFrequency.append(contfreq)
    return genotypes, CountFrequency




def merge_vcf_to_genotype(vcf_file, out_file, freq_file):
    """
    vcf_file:
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M628-0.bam     /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M628-1.bam     /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M628-2.bam
    1       10331   TRA000SUR       N       N[X:449317[     .       PASS    SUPP=1;SUPP_VEC=010;AVGLEN=100000;SVTYPE=BND;SVMETHOD=SURVIVORv2;CHR2=X;END=449317;CIPOS=0,0;CIEND=0,0;STRANDS=+-       GT:PSV:LN:DR:ST:TY:CO   ./.:NaN:0:0,0:--:NaN:NaN        1/1:NA:438986:0,2:+-:TRA:1_10331-X_449317       ./.:NaN:0:0,0:--:NaN:NaN
    1       10469   TRA001SUR       N       N[11:175280[    .       PASS    SUPP=2;SUPP_VEC=110;AVGLEN=100000;SVTYPE=BND;SVMETHOD=SURVIVORv2;CHR2=11;END=175280;CIPOS=0,0;CIEND=0,0;STRANDS=+-      GT:PSV:LN:DR:ST:TY:CO   1/1:NA:164811:0,8:+-:TRA:1_10469-11_175280      1/1:NA:164811:0,4:+-:TRA:1_10469-11_175280      ./.:NaN:0:0,0:--:NaN:NaN
    1       10469   TRA002SUR       N       N[3:198174528[  .       PASS    SUPP=1;SUPP_VEC=100;AVGLEN=100000;SVTYPE=BND;SVMETHOD=SURVIVORv2;CHR2=3;END=198174528;CIPOS=0,0;CIEND=0,0;STRANDS=++    GT:PSV:LN:DR:ST:TY:CO   1/1:NA:198164059:0,5:++:TRA:1_10469-3_198174528 ./.:NaN:0:0,0:--:NaN:NaN        ./.:NaN:0:0,0:--:NaN:NaN


    out_file:
    Chr1    Pos1    Chr2    Pos2    Type    M628-0  M628-1  M628-2
    1   10331   X   449317  BND 0/0 1/1 0/0
    1   10469   11  175280  BND 1/1 1/1 0/0
    1   10469   3   198174528   BND 1/1 0/0 0/0

    freq_file:
    Chr1    Pos1    Chr2    Pos2    Type    M628-0  M628-1  M628-2
    1   10331   X   449317  BND 0,0,0.000   0,2,1.000   0,0,0.000
    1   10469   11  175280  BND 0,8,1.000   0,4,1.000   0,0,0.000
    1   10469   3   198174528   BND 0,5,1.000   0,0,0.000   0,0,0.000    
    """
    vcf_h = open(vcf_file, "r")
    out_h = open(out_file, "w")
    freq_h = open(freq_file, "w")

    for line in vcf_h:
        line = line.strip()
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            headers = line.split("\t")
            ### sampels = ["/home/wuzhikun/Project/NanoTrio/mapping/minimap2/M628-0.bam", ]
            samples = headers[9:]
            sampleNames = [s.split("/")[-1].rstrip(".bam").split("_")[0].split(".")[0] for s in samples]

            ### output headers
            out_h.write("Chr1\tPos1\tChr2\tPos2\tSVlength\tType\t%s\n" % "\t".join(sampleNames))
            freq_h.write("Chr1\tPos1\tChr2\tPos2\tSVlength\tType\t%s\n" % "\t".join(sampleNames))
        else:
            lines = line.split("\t")
            Chr1, Start = lines[:2]
            AltType = lines[4]
            ### get information of another breakpoint
            Infor = lines[7]
            SVType, Chr2, End= Infor_target_values(Infor, "SVTYPE,CHR2,END")
            SVlength = Infor_target_values(Infor, "AVGLEN")
            # SVlength = str(abs(float(SVlength)))
            SVlength = str(abs(int(SVlength)))

            
            ### get the genotype information for samples
            Format = lines[8]
            GenoInfor = lines[9:]

            genotypes, CountFrequency = parse_sv_genotype(Format, GenoInfor)

            ### output the genotypes
            out_h.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Chr1, Start, Chr2, End, SVlength, SVType, "\t".join(genotypes)))

            ### output the count and frequency to another file
            freq_h.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Chr1, Start, Chr2, End, SVlength, SVType, "\t".join(CountFrequency)))
    vcf_h.close()
    out_h.close()
    freq_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the genotype based on the merged SV derived from SURVIVOR.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-f", "--frequency", help="The output count and frequency file.")
    args = parser.parse_args()
    merge_vcf_to_genotype(args.vcf, args.out, args.frequency)

if __name__ == "__main__":
    main()


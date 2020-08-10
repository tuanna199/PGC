#!/usr/bin/python
from __future__ import division
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
import argparse
import sys
import math

#usage: python ~/github/NanoHub/src/NanoHub/SVFiltReads.py --vcf M671-0.vcf --quality  ~/Project/NanoTrio/QualityControl/Samples_quality_summary.xls --ratioThreshold 0.2  --column Clean_total_base --out temp.vcf

def sample_sequence_depth(quality_file, column):
    SampleReads = {}
    quality_h = open(quality_file, "r")
    headers = quality_h.readline().strip().split("\t")
    try:
        readsIndex = headers.index(column)
    except KeyError:
        print("Please check whether the target column name %s exists in header of file %s." % (column, quality_file))
        sys.exit(1)
    for line in quality_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        reads = lines[readsIndex]
        SampleReads[sample] = reads
    quality_h.close()
    return SampleReads




def filt_SV_base_reads(quality_file, sv_file, out_file, ratioThreshold, column, support):
    support = int(support)

    ratioThreshold = float(ratioThreshold)

    # SampleReads = sample_sequence_depth(quality_file, column)

    # ### get the sample base
    # sampleName = sv_file.split("/")[-1].split(".")[0].split("_")[0]
    # try:
    #     readBase = int(SampleReads[sampleName])
    # except KeyError:
    #     print("Please check whther the sample name %s exists in the quality file %s." % (sampleName, quality_file))
    #     sys.exit(1)

    # depth = readBase / 3000000000 


    TypeList = ["INS", "DEL", "INV", "DUP"]
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

    in_h = open(sv_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Chr = lines[0]

            Infor = lines[7]
            # ### filt based on support read number
            # RE = Infor_target_values(Infor, "RE")
            # RE = int(RE)
            # AF = Infor_target_values(Infor, "AF")
            # AF = float(AF)
            SVType = Infor_target_values(Infor, "SVTYPE")
            if Chr in CHRS and SVType in TypeList:
                Format = lines[8]
                sampleValue = lines[9]
                if "AD" in Format:
                    AD = parse_genotype_format(Format, sampleValue, "AD")
                    DR, DV = AD.split(",")
                elif "DR" in Format:
                    DR, DV = parse_genotype_format(Format, sampleValue, "DR,DV")
                DR = sum([int(i) for i in DR.split(",")])
                DV = sum([int(i) for i in DV.split(",")])
                totalReads = DR + DV

                RE = DV
                AF = DV / totalReads
                # if RE >= ratioThreshold * totalReads * depth / 10 and totalReads >= depth * 0.4 and AF >= 0.2:
                if RE >= support:  #totalReads >= depth * 0.4 and 
                    out_h.write("%s\n" % line)
                else:
                    print(line)

    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Filt SV file based on the number of support records.")
    parser.add_argument("-v", "--vcf", help="The input sv file with vcf format.")
    parser.add_argument("-q", "--quality", help="The input sample read quality file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-r", "--ratioThreshold", default=0.2, help="The threshold for number of SV support record.")
    parser.add_argument("-c", "--column", default="Clean_total_base", help="The target column of quality file.")
    parser.add_argument("-s", "--support", default=2,  help="The read support number of SV.")
    args = parser.parse_args()
    filt_SV_base_reads(args.quality, args.vcf, args.out, args.ratioThreshold, args.column, args.support)



if __name__ == "__main__":
    main()


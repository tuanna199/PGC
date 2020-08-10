#!/usr/bin/env python
from BaseFunc import Infor_target_values
import argparse
import sys

#usage: python ~/github/TrioWGS/src/TrioWGS/SVvcf2bed.py --vcf M625_denovo_SV_proband.vcf --out M625_denovo_SV_proband.bed --distance 0


__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.04.14"

def SV_vcf_to_bed(sv_vcf, out_file, distance, variantType, Select, method):
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

    method = method.lower()

    distance = int(distance)
    in_h = open(sv_vcf, "r")
    out_h = open(out_file, "w")

    for line in in_h:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        else:
            lines = line.split("\t")
            Chr, Start, SVid = lines[:3]
            Start = int(Start)
            Ref = lines[3]
            ### Type is alternative allele in SNV and InDel calling
            Type = lines[4]
            Type = Type.lstrip("<").strip(">")
            Infor = lines[7]
            if variantType.lower() == "sv":
                SVType, Chr2, End= Infor_target_values(Infor, "SVTYPE,CHR2,END")
                End = int(End)

                if "SVLEN" in Infor:
                    SVLen = Infor_target_values(Infor, "SVLEN")
                elif "AVGLEN" in Infor:
                    SVLen = Infor_target_values(Infor, "AVGLEN")
                SVLen = SVLen.lstrip("-")

                ### it will go wrong then conduct IGV when "/" is in the name
                if "/" in SVType:
                    SVType = SVType.replace("/", "_")
                ### skip if the structural variantion is in different chromosomes

                ### just select chromosome
                if Select.lower() == "chr" or  Select.lower() == "chrs":
                    if Chr in CHRS:
                        if Chr == Chr2:
                            newStart = Start - distance
                            if newStart < 0:
                                newStart = 0

                        


                            if method == "region":
                                newEnd = End + distance
                            elif method == "point":
                                ### insertion is a point
                                if "INS" in SVType:
                                    newEnd = Start + distance
                                else:
                                    newEnd = End + distance
                            else:
                                print("Please make sure method is 'region' or 'point'.")
                                sys.exit(1)

                            # name = "_".join([Chr, str(newStart), str(newEnd), SVType])
                            name = "%s_%s-%s_%s-%s-%s" % (Chr, str(Start), Chr2, str(End), SVLen, SVType)

                            out_h.write("%s\t%d\t%d\t%s\n" % (Chr, newStart,  newEnd, name))
                    else:
                        continue
                else:
                    if Chr == Chr2:
                        newStart = Start - distance
                        if newStart < 0:
                            newStart = 0
                        newEnd = End + distance

                        # name = "_".join([Chr, str(newStart), str(newEnd), SVType])
                        # name = "_".join([Chr, str(Start), str(End), SVType])
                        name = "%s_%s-%s_%s-%s-%s" % (Chr, str(Start), Chr2, str(End), SVLen, SVType)
                        out_h.write("%s\t%d\t%d\t%s\n" % (Chr, newStart,  newEnd, name))

            elif variantType.lower() == "cnv":
                SVType, End= Infor_target_values(Infor, "SVTYPE,END")
                End = int(End)

                ### skip if the structural variantion is in different chromosomes
                newStart = Start - distance
                if newStart < 0:
                    newStart = 0
                newEnd = End + distance

                # name = "_".join([Chr, str(newStart), str(newEnd), SVType])
                name = "_".join([Chr, str(Start), str(End), SVType])

                out_h.write("%s\t%d\t%d\t%s\n" % (Chr, newStart,  newEnd, name))
            elif variantType.lower() == 'snv' or variantType.lower() == 'indel' or variantType.lower() == 'snp':
                if distance == 0:
                    distance = 100

                newStart = Start - distance
                if newStart < 0:
                    newStart = 0
                newEnd = Start + distance
                # ### name with extend distance
                # name = "_".join([Chr, str(newStart), str(newEnd), variantType])

                altLen = max([len(Ref), len(Type)])
                name = "_".join([Chr, str(Start), str(Start + altLen - 1), variantType])

                out_h.write("%s\t%d\t%d\t%s\n" % (Chr, newStart,  newEnd, name))
            else:
                sys.exit("Please make sure that the variant type is snp, snv, indel or sv.")

    out_h.close()




def main():
    parser = argparse.ArgumentParser(description="Get the main information based on the vcf file derived from Sniffles.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-d", "--distance", default=0, help="The distance to extend for the region.")
    parser.add_argument("-t", "--variant", default="SV", help="The variant type, such as SV, SNV, SNP, INDEL.")
    parser.add_argument("-s", "--select", default="None", help="Whether select chromosomes or alternative contigs.")
    parser.add_argument("-m", "--method", default="point", help="The method, such as 'point' or 'region' for insertion.")
    args = parser.parse_args()
    SV_vcf_to_bed(args.vcf, args.out, args.distance, args.variant, args.select, args.method)



if __name__ == "__main__":
    main()





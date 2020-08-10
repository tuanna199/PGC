#!/usr/bin/python
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
import collections
import argparse
import sys
import os
import operator

#usage: python ~/github/NanoHub/src/NanoHub/MergeSVOriginal.py --merge /home/wuzhikun/Project/NanoTrio/Multiple/M416-1/M416-1_common_SV.vcf --sniffle /home/wuzhikun/Project/NanoTrio/Multiple/SVFilt/M416-1_sniffles.vcf --nanovar /home/wuzhikun/Project/NanoTrio/Multiple/SVFilt/M416-1_nanovar.vcf --out temp.vcf

def nearst_region(Format, Start, End, geno):

    GenoDistance = {}
    genos = geno.split(";")
    genoLen = len(genos)
    if genoLen == 1:
        Tag1 = parse_genotype_format(Format, genos[0], "CO")
        Length = parse_genotype_format(Format, genos[0], "LN")
        Length = Length.lstrip("-")
        Type = parse_genotype_format(Format, genos[0], "TY")
    else:
        for g in genos:
            CO = parse_genotype_format(Format, g, "CO")
            regions = CO.split("-")
            Start2 =  int(regions[0].split("_")[-1])
            End2 =  int(regions[1].split("_")[-1])
            distance = abs(Start - Start2) + abs(End - End2)
            GenoDistance[g] = distance
        sortedDistance = sorted(GenoDistance.items(), key=operator.itemgetter(1))
        targetGeno = sortedDistance[0][0]
        Tag1 = parse_genotype_format(Format, targetGeno, "CO")
        Length = parse_genotype_format(Format, targetGeno, "LN")
        Length = Length.lstrip("-")
        Type = parse_genotype_format(Format, targetGeno, "TY")

    Tag = "%s-%s-%s" % (Tag1, Length, Type)
    return Tag


def tag_record(sv_file):
    TagRecord = {}

    sv_h = open(sv_file, "r")
    for line in sv_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start = lines[:2]
            Infor = lines[7]
            Chr2, End, Length, Type = Infor_target_values(Infor, "CHR2,END,SVLEN,SVTYPE")
            Length = Length.lstrip("-")
            Tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr2, End, Length, Type)
            TagRecord[Tag] = line
    sv_h.close()
    return TagRecord



def merge_SV_original(merged_SV, sniffle_file, nanovar_file, out_file, method):
    SniffleRecord = tag_record(sniffle_file)


    out_h = open(out_file, "w")
    merge_h = open(merged_SV, "r")

    method = method.lower()
    if method == "common":
        
        for line in merge_h:
            line = line.strip()
            if line.startswith("##"):
                out_h.write("%s\n" % line)
            elif line.startswith("#"):
                out_h.write("%s\n" % "\t".join(line.split("\t")[:-2]))
            else:
                lines = line.split("\t")
                Chr = lines[0]
                Start = int(lines[1])
                Infor = lines[7]
                End = int(Infor_target_values(Infor, "END"))
                Format = lines[8]
                Genotypes = lines[9:]
                genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genotypes]
                zeroNum = genos.count("./.")
                nonzeroNum = len(genos) - zeroNum

                if nonzeroNum == len(genos):
                    Tag = nearst_region(Format, Start, End, Genotypes[0])
                    try:
                        record = SniffleRecord[Tag]
                    except KeyError:
                        print("Please check whether the tag %s is in file %s." % (Tag, sniffle_file))
                        sys.exit(1)

                    out_h.write("%s\n" % record)

    elif method == "dominant":

        ### nanovar record
        NanovarRecord = tag_record(nanovar_file)
        
        for line in merge_h:
            line = line.strip()
            if line.startswith("##"):
                out_h.write("%s\n" % line)
            elif line.startswith("#"):
                out_h.write("%s\n" % "\t".join(line.split("\t")[:-2]))
            else:
                lines = line.split("\t")
                Chr = lines[0]
                Start = int(lines[1])
                Infor = lines[7]
                End = int(Infor_target_values(Infor, "END"))
                Format = lines[8]
                Genotypes = lines[9:]
                genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genotypes]
                zeroNum = genos.count("./.")
                nonzeroNum = len(genos) - zeroNum

                if nonzeroNum >= 2:
                    if genos[0] != "./.":
                        Tag = nearst_region(Format, Start, End, Genotypes[0])
                        try:
                            record = SniffleRecord[Tag]
                        except KeyError:
                            print("Please check whether the tag %s is in file %s." % (Tag, sniffle_file))
                            sys.exit(1)
                    elif genos[1] != "./.":
                        Tag = nearst_region(Format, Start, End, Genotypes[1])

                        try:
                            record = NanovarRecord[Tag]
                        except KeyError:
                            print("Please check whether the tag %s is in file %s." % (Tag, nanovar_file))
                            sys.exit(1)
                    else:
                        print("Please check whether the genotypes of first two of merged SV are valid")
                        sys.exit(1)

                    out_h.write("%s\n" % record)

    elif method == "denovo":
        for line in merge_h:
            line = line.strip()
            if line.startswith("##"):
                out_h.write("%s\n" % line)
            elif line.startswith("#"):
                out_h.write("%s\n" % "\t".join(line.split("\t")[:-2]))
            else:
                lines = line.split("\t")
                Chr = lines[0]
                Start = int(lines[1])
                Infor = lines[7]
                End = int(Infor_target_values(Infor, "END"))
                Format = lines[8]
                Genotypes = lines[9:]
                genos = [parse_genotype_format(Format, g.split(";")[0], "GT") for g in Genotypes]
                zeroNum = genos.count("./.")
                nonzeroNum = len(genos) - zeroNum

                if nonzeroNum == 1 and genos[0] != "./.":
                    Tag = nearst_region(Format, Start, End, Genotypes[0])
                    try:
                        record = SniffleRecord[Tag]
                    except KeyError:
                        print("Please check whether the tag %s is in file %s." % (Tag, sniffle_file))
                        sys.exit(1)

                    out_h.write("%s\n" % record)
    else:
        print("Please make sure the method is 'common', 'dominant' or 'denovo'.")
        sys.exit(1)

    merge_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Convert the BND to SV types for NanoSV results.")
    parser.add_argument("-v", "--merge", help="The input merged vcf file.")
    parser.add_argument("-s", "--sniffle", help="The input sniffle SV file.")
    parser.add_argument("-n", "--nanovar", help="The input nanovar SV file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-m", "--method", help="The method, such as 'common', 'dominant' or 'denovo'.")
    args = parser.parse_args()
    merge_SV_original(args.merge, args.sniffle, args.nanovar, args.out, args.method)

if __name__ == "__main__":
    main()








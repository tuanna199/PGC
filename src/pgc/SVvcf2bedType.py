#!/usr/bin/env python
from BaseFunc import Infor_target_values
import argparse
import collections
import sys

#usage: python ~/github/NanoHub/src/NanoHub/SVvcf2bedType.py --vcf M671-2.vcf --out temp --distance 0  --insDistance 50 --select chr


__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.10.14"


def output_SV_type_region(TypeRecords, SVType, SVLength, Chr, Start, Chr2, End, distance, insDistance, del_h, ins_h, dup_h, inv_h, tra_h, other_h):

    newStart = Start - distance
    if newStart < 0:
        newStart = 0
    newEnd = End + distance

    # name = "_".join([Chr, str(Start), str(End), SVType])
    name = "%s_%s-%s_%s-%s-%s" % (Chr, str(Start), Chr2, str(End), str(SVLength), SVType)

    if SVType == "TRA":
        # name = "_".join([Chr, str(Start), str(Chr2), str(End), SVType])
        record = "%s\t%d\t%s\t%d\t%s\n" % (Chr, newStart, Chr2, newEnd, name)
    elif SVType == "INS":
        ### extand the distance for insertion
        ### end may be smaller than start, so select start as end for ins
        # record = "%s\t%d\t%d\t%s\n" % (Chr, newStart - insDistance, newStart + insDistance, name)
        newEnd = newStart + int(float(SVLength))
        record = "%s\t%d\t%d\t%s\n" % (Chr, newStart, newEnd, name)
    else:
        record = "%s\t%d\t%d\t%s\n" % (Chr, newStart, newEnd, name)
    TypeRecords[SVType].append(record)

    return name





def out_record(sortedRecords, out_h):
    for record in sortedRecords:
        out_h.write(record)



def SV_vcf_to_bed(sv_vcf, out_pre, distance, insDistance, Select, Length_out):
    insDistance = int(insDistance)
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"] #, "Y"

    TypeRecords = collections.defaultdict(list)

    distance = int(distance)
    in_h = open(sv_vcf, "r")
    inv_h = open(out_pre + "_INV.bed", "w")
    dup_h = open(out_pre + "_DUP.bed", "w")
    del_h = open(out_pre + "_DEL.bed", "w")
    ins_h = open(out_pre + "_INS.bed", "w")
    tra_h = open(out_pre + "_TRA.bed", "w")
    other_h = open(out_pre + "_OTHER.bed", "w")

    length_h = open(Length_out, "w")

    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start, SVid = lines[:3]
            Start = int(Start)
            Ref = lines[3]
            ### Type is alternative allele in SNV and InDel calling
            Type = lines[4]
            # Type = Type.lstrip("<").strip(">")
            Infor = lines[7]
            Type = Infor_target_values(Infor, "SVType")

            if "SVLEN" in Infor:
                SVLength = Infor_target_values(Infor, "SVLEN")
            elif "AVGLEN" in Infor:
                SVLength = Infor_target_values(Infor, "AVGLEN")
            else:
                print("Please check whether the tag name 'SVLEN' or 'AVGLEN' in the record %s." % line)
                sys.exit(1)
            SVLength = SVLength.lstrip("-")
            # SVLength = int(float(SVLength))
            # SVLength = "%.0f" % SVLength


            SVType, Chr2, End= Infor_target_values(Infor, "SVTYPE,CHR2,END")
            End = int(End)

            # ### it will go wrong then conduct IGV when "/" is in the name
            # if "/" in SVType:
            #     SVType = SVType.replace("/", "_")
            ### skip if the structural variantion is in different chromosomes

            ### just select chromosome
            if Select.lower() == "chr" or  Select.lower() == "chrs":
                if Chr in CHRS and Chr2 in CHRS:
                    SVName = output_SV_type_region(TypeRecords, SVType, SVLength, Chr, Start, Chr2, End, distance, insDistance, del_h, ins_h, dup_h, inv_h, tra_h, other_h)
            else:
                SVName = output_SV_type_region(TypeRecords, SVType, SVLength, Chr, Start, Chr2, End, distance, insDistance, del_h, ins_h, dup_h, inv_h, tra_h, other_h)


            ### output the SV name and length
            length_h.write("%s\t%s\n" % (SVName, SVLength))
    length_h.close()

    ### output
    for SVType in TypeRecords:
        sortedRecords = sorted(TypeRecords[SVType])

        if SVType == "DEL":
            out_record(sortedRecords, del_h)
        elif SVType == "INS":
            out_record(sortedRecords, ins_h)
        elif SVType == "DUP":
            out_record(sortedRecords, dup_h)
        elif SVType == "INV":
            out_record(sortedRecords, inv_h)
        elif SVType == "TRA":
            out_record(sortedRecords, tra_h)
        else:
            out_record(sortedRecords, other_h)
    del_h.close()
    ins_h.close()
    dup_h.close()
    inv_h.close()
    tra_h.close()
    other_h.close()




def main():
    parser = argparse.ArgumentParser(description="Get the main information based on the vcf file derived from Sniffles.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-d", "--distance", default=0, help="The distance to extend for the region.")
    parser.add_argument("-i", "--insDistance", default=50, help="The extended distance for insertion.")
    parser.add_argument("-s", "--select", default="None", help="Whether select chromosomes or alternative contigs.")
    parser.add_argument("-l", "--lengthOut", help="The output file contain SV name and length.")
    args = parser.parse_args()
    SV_vcf_to_bed(args.vcf, args.out, args.distance, args.insDistance, args.select, args.lengthOut)



if __name__ == "__main__":
    main()





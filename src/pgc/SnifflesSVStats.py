#!/usr/bin/env python
import sys
import argparse
import collections

#usage: python ~/github/TrioWGS/src/TrioWGS/SnifflesSVStats.py --vcf M625-0_test,vcf --summary M625-0_summary.xls --record M625-0_record.xls --sample M625-0 --position SV_position.xls

__author__ = "Zhikun Wu"
__date__ = "2019.04.12"
__email__ = "598466208@qq.com"

def Infor_target_values(Infor, target):
    """
    Infor = "IMPRECISE;SVMETHOD=Snifflesv1.0.10;CHR2=1;END=181364;STD_quant_start=26.000000;STD_quant_stop=79.501572;Kurtosis_quant_start=-2.000000;Kurtosis_quant_stop=-1.999842;SVTYPE=DEL;RNAMES=a5c2a7ee-ce33-4dd3-8fac-3a18286ce747,b0fdea87-ced4-44a6-b6dc-b9071015fac0;SUPTYPE=AL;SVLEN=-182;STRANDS=+-;RE=2"

    target = "END"
    return:
    '181364'

    target = "END,SVTYPE"
    return:
    ["181364", "DEL"]
    """
    TagValues = {}
    Infors = Infor.split(";")
    for f in Infors:
        fs = f.split("=") 
        if len(fs) == 1:
            if fs[0] == "PRECISE":
                TagValues["PRECISE"] = True
            else:
                TagValues["PRECISE"] = False
        elif len(fs) == 2:
            tag, value = fs
            TagValues[tag] = value
        else:
            print("Please check and make sure the items of tag is no more than two.")
            sys.exit(1)
    ### get the target items
    Items = []
    targets = target.split(",")
    targets = [t.strip() for t in targets]
    for t in targets:
        try:
            value = TagValues[t]
            Items.append(value)
        except KeyError:
            print("Please check and make sure the given tag %s is in the Infor record %s." % (t, Infor))
            Items.append("0")
    if len(Items) == 1:
        return Items[0]
    else:
        return Items


def output_type_number(TypeCount, Type, targets, type_h):
    """
    TypeCount: {"Type": {"SVLEN": ["5", "10"], "RE": [2, 3]}}
    """

    targetLen = len(targets)
    try:
        TargetCount = TypeCount[Type]
        try:
            twoCounts = [TargetCount[t] for t in targets]
            VarCounts = len(twoCounts[0])
            for i in range(VarCounts):
                counts = [twoCounts[j][i] for j in range(targetLen)]
                counts = [str(c) for c in counts]
                type_h.write("%s\t%s\n" % (Type, "\t".join(counts)))
        except KeyError:
            print("Please check whether the targets %s is in the dict %s." % (targets, TypeCount))
            sys.exit(1)
    except KeyError:
        print("Please check the type %s is in the dict TypeCount %s." % (Type, TypeCount))
        sys.exit(1)
    return VarCounts


def Sniffles_SV_Stats(SV_file, out_file, length_file, targetItem, SV_pos, sampleName="Sample", TRA=None):
    """
    We just output record summary of the simple types, such as DEL, DUP, INS, INV, TRA.

    Types:
    DEL
    DEL/INV
    DUP
    DUP/INS
    INS
    INV
    INVDUP
    INV/INVDUP
    TRA
    """
    TypeCount = collections.defaultdict(lambda: collections.defaultdict(list))
    in_h = open(SV_file, "r")

    if TRA:
        tra_h = open(TRA, "w")
    
    pos_h = open(SV_pos, "w")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Pos = lines[:2]
            Type = lines[4]
            Type = Type.lstrip("<").rstrip(">")
            Infor = lines[7]
            ### get the values of tags
            Chr2, End, RE = Infor_target_values(Infor, "CHR2,END,RE")
            if "SVLEN" in Infor:
                SVLen = Infor_target_values(Infor, "SVLEN")
                SVLen = int(SVLen.lstrip("-"))
            elif "AVGLEN" in Infor:
                SVLen = Infor_target_values(Infor, "AVGLEN")
                SVLen = float(SVLen.lstrip("-"))



            if SVLen != "NA":
                TypeCount[Type]["SVLEN"].append(SVLen)
                TypeCount[Type]["RE"].append(int(RE))

                ### output the chr, start and end
                # if Chr == Chr2:
                pos_h.write("%s\t%s\t%s\t%s\n" % (Chr, Pos, End, Type))

            ### output TRA
            if Chr != Chr2:
                if TRA:
                    tra_h.write("%s\t%s\t%s\t%s\n" % (Chr, Pos, Chr2, End))

    pos_h.close()
    if TRA:
        tra_h.close()
    ### output each record for types and tags
    type_h = open(length_file, "w")
    targetItems = targetItem.split(",")
    targetItems = [t.strip() for t in targetItems]
    type_h.write("Type\t%s\n" % "\t".join(targetItems))
    DELNum = output_type_number(TypeCount, "DEL", targetItems, type_h)
    DUPNum = output_type_number(TypeCount, "DUP", targetItems, type_h)
    INSNum = output_type_number(TypeCount, "INS", targetItems, type_h)
    INVNum = output_type_number(TypeCount, "INV", targetItems, type_h)
    TRANum = output_type_number(TypeCount, "TRA", targetItems, type_h)
    INVDUPNum = output_type_number(TypeCount, "INVDUP", targetItems, type_h)
    type_h.close()

    ### output summary file
    out_h = open(out_file, "w")
    out_h.write("SampleID\tDEL\tDUP\tINS\tINV\tTRA\tINVDUP\tTotal\n")
    Total = DELNum + DUPNum + INSNum + INVNum + TRANum + INVDUPNum
    out_h.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (sampleName, DELNum, DUPNum, INSNum, INVNum, TRANum, INVDUPNum, Total))

def main():
    parser = argparse.ArgumentParser(description="Get the summary for structural variants called from Sniffles.")
    parser.add_argument("-v", "--vcf", help="the input vcf file.")
    parser.add_argument("-s", "--summary", help="Get the summary of different type SV.")
    parser.add_argument("-r", "--record", help="Output the record of length and read evidence.")
    parser.add_argument("-t", "--target", default="SVLEN,RE", help="The target tags which are separated with ','.")
    parser.add_argument("-p", "--position", help="The position record for SV.")
    parser.add_argument("-n", "--sample", default="Sample", help="The sample name.")
    parser.add_argument("-a", "--tra", help="The output record of tra.")
    args = parser.parse_args()
    if args.tra:
        Sniffles_SV_Stats(args.vcf, args.summary, args.record, args.target, args.position, sampleName=args.sample, TRA=args.tra)
    else:
        Sniffles_SV_Stats(args.vcf, args.summary, args.record, args.target, args.position, sampleName=args.sample, TRA=None)

if __name__ == "__main__":
    main()

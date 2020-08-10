#!/usr/bin/python
from __future__ import division
from BaseFunc import Infor_target_values
from BaseFunc import parse_genotype_format
import argparse
import sys
import math
import collections

### usage: python ~/github/NanoHub/src/NanoHub/SVFiltLength.py --vcf M416-0.vcf --out temp_out.vcf --length length.txt --number number.txt

def SV_tag(lines):
    ID = lines[0]
    Chr, Start = lines[:2]
    Infor = lines[7]
    ### change the ID based on the Infor
    Chr2 = Infor_target_values(Infor, "CHR2")
    End = Infor_target_values(Infor, "END")
    SVType = Infor_target_values(Infor, "SVTYPE")
    if "SVLEN" in Infor:
        SVLength = Infor_target_values(Infor, "SVLEN")
    elif "AVGLEN" in Infor:
        SVLength = Infor_target_values(Infor, "AVGLEN")
        
    if SVLength.startswith("-"):
        SVLength = SVLength.lstrip("-")


    if SVType == "BND":
        SVType = "TRA"

    if SVType == "TRA":
        SVLength = 0

    Tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr2, End, SVLength, SVType) ### just one chromosomes
    
    return Tag, float(SVLength)


def primer_type(SVTypes):
    """
    1033213 <DEL>
        820 <DEL/INV>
      21978 <DUP>
        445 <DUP/INS>
    1131560 <INS>
      17939 <INV>
      30658 <INVDUP>
         25 <INV/INVDUP>
      56824 <TRA>
    """
    Types = ["DEL", "INS", "INV", "DUP", "INVDUP"]
    type1, type2 = SVTypes
    typeIndex1 = Types.index(type1)
    typeIndex2 = Types.index(type2)
    if typeIndex1 < typeIndex2:
        return type1
    else:
        return type2


def vcf_filt_length(vcf_file, out_file, length_file, number_file, delThreshold, invThreshold):
    """
    >>> theArray = [['a','b','c'],['d','e','f'],['g','h','i']]
    >>> [*zip(*theArray)]
    [('a', 'd', 'g'), ('b', 'e', 'h'), ('c', 'f', 'i')]
    """
    delThreshold = int(delThreshold)
    invThreshold = int(invThreshold)

    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]

    vcf_h = open(vcf_file, "r")
    out_h = open(out_file, "w")

    ChrTypeLens = collections.defaultdict(lambda: collections.defaultdict(list))

    TypeSet = set()

    for line in vcf_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Chr = lines[0]
            if Chr in CHRS:
                Infor = lines[7]
                SVType = Infor_target_values(Infor, "SVTYPE")

                ### multiple types for one record
                ### change to one type
                if "/" in SVType:
                    SVTypes = SVType.split("/")
                    if len(SVTypes) == 2:
                        SVType = primer_type(SVTypes)

                        infors = Infor.split(";")
                        for i, v in enumerate(infors):
                            if "SVTYPE" in v:
                                typeIndex = i
                        newInfor = "%s;SVTYPE=%s;%s" % (";".join(infors[:typeIndex]), SVType, ";".join(infors[typeIndex+1:]))
                        line = "%s\t<%s>\t%s\t%s\t%s" % ("\t".join(lines[:4]), SVType, "\t".join(lines[5:7]), newInfor, "\t".join(lines[8:]))
                    else:
                        print("Please make sure that there are at most two types for record %s." % line)
                        sys.exit(1)



                TypeSet.add(SVType)

                Tag, SVLength = SV_tag(lines)
                
                if SVType == "INS" or SVType == "DEL":
                    if SVLength < delThreshold:
                        out_h.write("%s\n" % line)
                        ChrTypeLens[Chr][SVType].append(SVLength)
                    else:
                        print(Tag)
                # elif SVType == "DUP" or SVType == "INV":
                elif SVType == "INV" or SVType == "DUP" or SVType == "INVDUP":
                    if SVLength < invThreshold:
                        out_h.write("%s\n" % line)
                        ChrTypeLens[Chr][SVType].append(SVLength)
                    else:
                        print(Tag)
                else:
                    out_h.write("%s\n" % line)
                    ChrTypeLens[Chr][SVType].append(SVLength)
    out_h.close()


    ### output length
    length_h = open(length_file, "w")
    number_h = open(number_file, "w")
    sortedChr = sorted(list(ChrTypeLens.keys()))

    TypeList = list(TypeSet)

    length_h.write("Chr\t%s\tTotal\n" % "\t".join(TypeList))
    number_h.write("Chr\t%s\tTotal\n" % "\t".join(TypeList))

    LENGTHTOTAL = []
    NUMTOTAL = []

    for c in sortedChr:
        LENGTH = []
        NUM = []
        for t in TypeList:
            if t in ChrTypeLens[c]:
                length = ChrTypeLens[c][t]
            else:
                length = []

            lengthSum = sum(length)
            LENGTH.append(lengthSum)

            number = len(length)
            NUM.append(number)

        Total = sum(LENGTH)
        LENGTHTOTAL.append(LENGTH)
        LENGTH = [str(int(l)) for l in LENGTH]
        length_h.write("%s\t%s\t%d\n" % (c, "\t".join(LENGTH), Total))

        TotalNum = sum(NUM)
        NUMTOTAL.append(NUM)
        NUM = [str(n) for n in NUM]
        number_h.write("%s\t%s\t%d\n" % (c, "\t".join(NUM), TotalNum))

    ### list
    ### [list(i) for i in zip(*theArray)]

    ### tuple
    TRANSLENGTH = [*zip(*LENGTHTOTAL)]
    TRANSNUM = [*zip(*NUMTOTAL)]
    TRANSLENGTHSUM = [int(sum(list(n))) for n in TRANSLENGTH]
    TRANSLENGTHSTR = [str(n) for n in TRANSLENGTHSUM]
    TRANSNUMSUM = [int(sum(list(n))) for n in TRANSNUM]
    TRANSNUMSTR = [str(n) for n in TRANSNUMSUM]
    print(TRANSNUMSTR)


    length_h.write("Total\t%s\t%d\n" % ("\t".join(TRANSLENGTHSTR), sum(TRANSLENGTHSUM)))
    number_h.write("Total\t%s\t%d\n" % ("\t".join(TRANSNUMSTR), sum(TRANSNUMSUM)))
    length_h.close()
    number_h.close()


def main():
    parser = argparse.ArgumentParser(description="Filt the SV based on the length.")
    parser.add_argument("-v", "--vcf", help="The input sv file with vcf format.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-l", "--length", help="Output SV length summary file.")
    parser.add_argument("-n", "--number", help="Output SV number summary file.")
    parser.add_argument("-d", "--delThreshold", default=2000000, help="The length threshold for deletion and insertion.")
    parser.add_argument("-i", "--invThreshold", default=20000000, help="The length threshold for inversion and duplication.")
    args = parser.parse_args()
    vcf_filt_length(args.vcf, args.out, args.length, args.number, args.delThreshold, args.invThreshold)




if __name__ == "__main__":
    main()
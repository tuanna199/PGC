#!/usr/bin/python 
from __future__ import division
from BaseFunc import Infor_target_values
import argparse
import re
import collections


#usage:  python /home/wuzhikun/github/NanoHub/src/NanoHub/SVAlleleCalculate.py --vcf /home/wuzhikun/Project/Population/population/Merge/Sample_common_SV_20191231.vcf  --type /home/wuzhikun/Project/Population/population/Stats/Sample_type_freq.xls --poly /home/wuzhikun/Project/Population/population/Stats/Sample_SV_poly.xls --sample /home/wuzhikun/Project/Population/population/Stats/Sample_SV_types.xls --outPrefix /home/wuzhikun/Project/Population/population/Stats/Sample_SV_record

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.05.27"

# def supp_freq(SUPP):
#     """
#     100110
#     """
#     SuppList = [int(s) for s in SUPP]
#     SuppSum = sum(SuppList)
#     SuppLength = len(SuppList)
#     SuppFreq = SuppSum / SuppLength
#     return SuppFreq, SuppSum



def sample_add_SV(SUPP, SampleID, ID):
    SuppList = [int(s) for s in SUPP]
    for s in range(len(SuppList)):
        ss = SuppList[s]
        if ss == 1:
            SampleID[s].append(ID)


def Tag_freq_class(tag_file):
    TagPolyClass = {}
    tag_h = open(tag_file, "r")
    header = tag_h.readline()
    for line in tag_h:
        lines = line.strip().split("\t")
        Class = lines[0]
        Tag = lines[1]
        tags = Tag.split(",")
        for t in tags:
            TagPolyClass[t] = Class
    tag_h.close()
    return TagPolyClass



def Type_freq_number(TagPolyClass, TypeID, IDRecord):
    """
    TypePolyNumber:
    {'DEL': {'Single': 1, 'Poly': 1}, 'INS': {'Poly': 4, 'Single': 1}}

    IDPolyType:
    {'DEL000SUR': 'Single', 'DEL001SUR': 'Poly', 'INS002SUR': 'Poly', 'INS003SUR': 'Poly', 'INS004SUR': 'Single', 'INS005SUR': 'Poly', 'INS006SUR': 'Poly'}

    PolyRecords:
    {'Single': ['1\t66288\tDEL000SUR\tN\t<DEL>\t.\tPASS\tSUPP=1',]}
    """
    IDPolyType = {}
    TypePolyNum = collections.defaultdict(lambda: collections.defaultdict(list))
    TypePolyNumber = collections.defaultdict(dict)
    PolyRecords = collections.defaultdict(list)


    for Type in TypeID:
        IDs = TypeID[Type]
        for ID in IDs:
            Poly = TagPolyClass[ID]
            TypePolyNum[Type][Poly].append(ID)
            IDPolyType[ID] = Poly



    for Type in TypePolyNum:
        for Poly in TypePolyNum[Type]:
            numberList = TypePolyNum[Type][Poly]
            number = len(numberList)
            TypePolyNumber[Type][Poly] = number

    ### get the records for different poly type
    for ID in IDPolyType:
        Type = IDPolyType[ID]
        record = IDRecord[ID]
        PolyRecords[Type].append(record)

    return TypePolyNumber, IDPolyType, PolyRecords


def sample_type_number(IDPolyType, IDType, SampleID):
    """
    SampleTypes:
    {110: Counter({'DEL': 1}), 18: Counter({'DEL': 1}), }

    SamplePolys:
    {110: Counter({'Single': 1}), 18: Counter({'Poly': 1}), }
    """
    SampleTypes = collections.defaultdict(lambda: collections.Counter())
    SamplePolys = collections.defaultdict(lambda: collections.Counter())
    for sample in SampleID:
        IDs = SampleID[sample]
        for ID in IDs:
            Poly = IDPolyType[ID]
            Type = IDType[ID]
            SampleTypes[sample][Type] += 1
            SamplePolys[sample][Poly] += 1
    return SampleTypes, SamplePolys











def assign_SV_type(tag_category, vcf_file, type_out, poly_out, sample_out, recordPrefix):
    """
    inversion:
    N]1:23456789]  or [1:23456789[N

    deletion:
    N[1:23456789[

    duplications:
    ]1:23456789]N

    type_out:
    Sample  DEL DUP INS INV TRA Total
    M509-1  17726   396 10615   449 2410    31596
    M509-2  8715    321 9901    511 1724    21172

    poly_out:
    Sample  Single  Poly    Major   Common  Total
    M509-1  17451   3291    9386    1468    31596
    M509-2  7424    2959    9321    1468    21172

    sample_out:
    Type    Single  Poly    Major   Common
    DEL 534535  3634    5918    517
    DUP 10286   138 124 7
    INS 105397  3856    5988    897
    INV 86812   92  139 9
    TRA 89040   601 550 38
    """
    TagPolyClass = Tag_freq_class(tag_category)

    PolyList = ["Singleton", "Rare", "Low", "Common"]
    TypeList = ["DEL", "INS", "DUP", "INV"]

    Types = set()

    SVFrequency = {}
    SVCount = {}
    IDType = {}
    TypeID = collections.defaultdict(list)
    SampleID = collections.defaultdict(list)
    IndexSample = {}

    IDRecord = {}

    HEADERS = []

    vcf_h = open(vcf_file, "r")
    for line in vcf_h:
        line = line.strip()
        if line.startswith("##"):
            HEADERS.append(line)
        elif line.startswith("#"):
            HEADERS.append(line)
            headers = line.strip().split("\t")
            samples = headers[9:]
            samples = [s.split("/")[-1].split(".")[0] for s in samples]
            ### index for sample
            for s in range(len(samples)):
                ss = samples[s]
                IndexSample[s] = ss
        else:
            lines = line.split("\t")
            Chr1 = lines[0]
            Start = lines[1]
            # ID = lines[2]
            Ref = lines[3]
            Alt = lines[4]
            Infor = lines[7]

            # AltType = Alt.lstrip("<").rstrip(">")

            Chr2, SVType, SVLen, End = Infor_target_values(Infor, "CHR2,SVTYPE,AVGLEN,END")
            SVLen = SVLen.lstrip("-")

            Types.add(SVType)

            SUPP = Infor_target_values(Infor, "SUPP_VEC")

            # SuppFreq, SuppSum = supp_freq(SUPP)


            Tag = "%s_%s-%s_%s-%s-%s" % (Chr1, Start, Chr2, End, SVLen, SVType)
            ID = Tag

            IDType[ID] = SVType
            IDRecord[ID] = line
            # SVFrequency[ID] = SuppFreq
            # SVCount[ID] = SuppSum
            ##################### select SV  type #################
            # if SVType in TypeList:
            #     TypeID[SVType].append(ID)

            TypeID[SVType].append(ID)
            ### sample add SV ID
            sample_add_SV(SUPP, SampleID, ID)
    vcf_h.close()

    TypePolyNumber, IDPolyType, PolyRecords = Type_freq_number(TagPolyClass, TypeID, IDRecord)
    # TypePolyNumber, IDPolyType, PolyRecords = Type_freq_number(TypeID, SVFrequency, SVCount, IDRecord)

    SampleTypes, SamplePolys = sample_type_number(IDPolyType, IDType, SampleID)




    ### output file
    type_h = open(type_out, "w")
    
    type_h.write("Type\t%s\n" % "\t".join(PolyList))


    print(sorted(list(TypePolyNumber.keys())))
    ## ['DEL', 'DUP', 'INS', 'INV', 'INVDUP', 'TRA']

    for Type in sorted(list(TypePolyNumber.keys())):

        Numbers = []
        for poly in PolyList:
            if poly in TypePolyNumber[Type]:
                num = TypePolyNumber[Type][poly]
            else:
                num = 0
            Numbers.append(num)

        Numbers = [str(n) for n in Numbers]
        type_h.write("%s\t%s\n" % (Type, "\t".join(Numbers)))

    poly_h = open(poly_out, "w")
    poly_h.write("Sample\t%s\tTotal\n" % "\t".join(PolyList))
    for sampleIndex in sorted(list(SamplePolys.keys())):
        sample = IndexSample[sampleIndex]
        Numbers = [SamplePolys[sampleIndex][poly] for poly in PolyList]
        numberSUM = sum(Numbers)
        Numbers = [str(n) for n in Numbers]
        poly_h.write("%s\t%s\t%s\n" % (sample, "\t".join(Numbers), str(numberSUM)))
    

    sortedTypes = sorted(list(Types))
    sample_h = open(sample_out, "w")
    sample_h.write("Sample\t%s\tTotal\n" % "\t".join(sortedTypes))
    for sampleIndex in sorted(list(SampleTypes.keys())):
        sample = IndexSample[sampleIndex]
        Numbers = [SampleTypes[sampleIndex][t] for t in sortedTypes]
        numberSUM = sum(Numbers)
        Numbers = [str(n) for n in Numbers]
        sample_h.write("%s\t%s\t%s\n" % (sample, "\t".join(Numbers), str(numberSUM)))

    type_h.close()
    poly_h.close()
    sample_h.close()


    ### output records for different poly types
    for t in PolyList:
        r_h = open(recordPrefix + "_" + t + ".vcf", "w")
        records = PolyRecords[t]
        r_h.write("%s\n%s\n" % ("\n".join(HEADERS), "\n".join(records)))
        r_h.close()



def main():
    parser = argparse.ArgumentParser(description="Convert the BND to SV types for NanoSV results.")
    parser.add_argument("-c", "--category", help="The tag and category file.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-t", "--type", help="The output file with type and nmber.")
    parser.add_argument("-p", "--poly", help="The output file with poly and nmber.")
    parser.add_argument("-s", "--sample", help="The output file with sample and type nmber.")
    parser.add_argument("-o", "--outPrefix", help="The output vcf file with this prefix name.")
    args = parser.parse_args()
    assign_SV_type(args.category, args.vcf, args.type, args.poly, args.sample, args.outPrefix)

if __name__ == "__main__":
    main()


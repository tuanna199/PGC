#!/usr/bin/env python
from __future__ import division
import sys
import collections
import bisect
import re

###################### meta target ########################
def meta_target_file(targetList, meta_file, pair_name):
    tumorName = PairInfor(meta_file).pair_tumor(pair_name)
    normalName = PairInfor(meta_file).pair_normal(pair_name)
    for target in targetList:
        if tumorName in target:
            tumorFile = target
        elif normalName in target:
            normalFile = target
    return normalFile, tumorFile



###################### random genotype ###############
def random_genotype(in_file, out_file):
    """
    in_file:
    Chr1    Pos1    Chr2    Pos2    SVlength        Type    CN001   CN002   CN003   CN004   CN005   CN006   CN007   CN008   CN009   CN010   CN011   CN012
    1       66288   1       66527   239.0   DEL     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       67910   1       68341   431.0   DEL     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       83968   1       84057   89.0    INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       88684   1       88831   147.0   DEL     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       88889   1       88956   66.0    INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       90312   1       90388   74.3    INS     0/0     0/0     0/1     0/0     0/0     0/1     0/1     0/1     0/0     0/0     0/0     0/0     0/0
    1       95077   1       95176   101.0   INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       136348  1       137277  928.5   DEL     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       136512  1       136612  95.0    INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    """
    import random
    import collections
    GenoList = ["0/0", "0/1", "1/1"]
    SampleGenoCount = collections.defaultdict(lambda: collections.Counter())
    in_h = open(in_file, "r")
    headers = in_h.readline().strip().split("\t")
    samples = headers[6:]
    sampleLem = len(samples)
    out_h = open(out_file, "w")
    out_h.write("%s\n" % "\t".join(samples))
    for line in in_h:
        lines = line.strip().split("\t")
        t = lines[:6]
        tt = "_".join(t)
        genotypes = []
        for s in range(sampleLem):
            sample = samples[s]
            geno = random.choice(GenoList)
            genotypes.append(geno)
            SampleGenoCount[sample][geno] += 1
        out_h.write("%s\t%s\n" % ("\t".join(t), "\t".join(genotypes)))
    in_h.close()
    out_h.close()
    SampleCount = []
    for s in samples:
        counts = []
        for c in GenoList:
            if c in SampleGenoCount[s]:
                count = SampleGenoCount[s][c]
            else:
                count = 0
            counts.append(count)
        SampleCount.append(counts)
    print(SampleCount)





###################### category tags ##############
def category_tags(category_file):
    """
    category_file:
    Category        Tags
    common  10_101728712-10_101729253-540.3-INS,10_116505797-10_116506120-321.9-DEL,10_124330938-10_124331032-91.7-INS,10_126417903-10_126418217-314.2-
    major   10_100093578-10_100093710-119.6-
    """
    TagCategory = {}
    CategoryCount = collections.Counter()

    in_h = open(category_file, "r")
    header = in_h.readline()
    for line in in_h:
        lines = line.strip().split("\t")
        cate, tag = lines[:2]
        tags = tag.split(",")
        for tag in tags:
            TagCategory[tag] = cate
            CategoryCount[cate] += 1
    in_h.close()
    return TagCategory, CategoryCount




def group_samples(meta_file, target):
    import collections
    import sys
    
    meta_h = open(meta_file, "r")
    headers = meta_h.readline().strip().split("\t")
    headers = [h.lower() for h in headers]
    try:
        groupIndex = headers.index("group")
        sampleIndex = headers.index("sample")
    except KeyError:
        print("Please check whether the header of file %s contain the column 'group' and 'sample'." % meta_file)
        sys.exit(1)

    GroupSamples = collections.defaultdict(list)

    for line in meta_h:
        lines = line.strip().split("\t")
        group = lines[groupIndex]
        sample = lines[sampleIndex]
        GroupSamples[group].append(sample)
    meta_h.close()

    try:
        samples = GroupSamples[target]
    except KeyError:
        print("Please check whether the target group %s is in the file %s." % (target, meta_file))
        sys.exit(1)
    return samples




def match_location(location):
    """
    a = "intron27"
    b = re.findall("(\D*)(\d*)", a)
    b = [('intron', '27'), ('', '')]

    a = "txStart"
    b = re.findall("(\D*)(\d*)", a)
    b = [('txStart', ''), ('', '')]
    """
    locations = location.split("-")
    locLen = len(locations)
    if locLen == 2:
        loc1, loc2 = locations
        locMatch1 = re.findall("(\D*)(\d*)", loc1)
        locMatch2 = re.findall("(\D*)(\d*)", loc2)
        locNum1 = locMatch1[0][1]
        locNum2 = locMatch2[0][1]
        if locNum1 == '' and locNum2 == '':
            newloc1 = loc1.lstrip("tx")
            newloc2 = loc2.lstrip("tx")
        elif locNum1 == '' and locNum2 != '':
            newloc1 = loc1.lstrip("tx")
            newloc2 = locMatch2[0][0]
        elif locNum1 != '' and locNum2 == '':
            newloc1 = locMatch1[0][0]
            newloc2 = loc2.lstrip("tx")
        elif locNum1 != '' and locNum2 != '':
            newloc1 = locMatch1[0][0]
            newloc2 = locMatch2[0][0]
        return newloc1, newloc2
    else:
        print("Please check whether the location of %s has start and end positions." % location)
        sys.exit(1)




def region_overlap_length(Start1, End1, Start2, End2):
    region = [Start1, End1]
    startIndex = bisect.bisect(region, Start2)
    endIndex = bisect.bisect(region, End2)
    if startIndex == 1 and endIndex == 2:
        overlap = (Start2, End1)
        overlapLen = End1 - Start2
    elif startIndex == 0 and endIndex == 1:
        overlap = (Start1, End2)
        overlapLen = End2 - Start1
    elif startIndex == 0 and endIndex == 2:
        overlap = (Start1, End1)
        overlapLen = End1 - Start1
    elif startIndex == 1 and endIndex == 1:
        overlap = (Start2, End2)
        overlapLen = End2 - End1
    else:
        print("Please check whether the region %s is overlap with %s." % ((Start1, End1), (Start2, End2)))
        sys.exit(1)

    return overlapLen




def column_index(records, column):
    column = column.lower()
    records = [r.lower() for r in records]
    try:
        columnIndex = records.index(column)
        return columnIndex
    except ValueError:
        print("Please check whether the column name %s is in record %s." % (column, records))
        sys.exit(1)



def Format_tag(Format, Geno):
    """
    Format: GT:PSV:LN:DR:ST:TY:CO
    geno: 1/1:NA:137:1,36:+-:INS,INS:1_136673-1_136891,1_136999-1_137209

    geno:
    1/1:NA:90:0,17:+-:INS:1_137053-1_137183;0/1:NA:188:29,19:+-:INS:1_136930-1_137112

    ['NaN-0-NaN', '1_10169-X_449438-0-TRA']
    """
    genos = Geno.split(";")

    Tags = []
    for geno in genos:
        LN = parse_genotype_format(Format, geno, "LN")
        TY = parse_genotype_format(Format, geno, "TY")
        CO = parse_genotype_format(Format, geno, "CO")

        ### set length of translocation as 0
        if TY == "TRA":
            LN = "0"

        tag = "%s-%s-%s" % (CO, LN, TY)
        Tags.append(tag)

    TagLen = len(Tags)
    if TagLen == 1:
        return Tags[0]
    else:
        return ";".join(Tags)



############################# SV tag of vcf from Sniffles ########################
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
    
    return Tag

#########################################################################

############################## overlap with chromosome regions ##################
def overlap_region(Start, End, region):
    startIndex = bisect.bisect(region, Start)
    endIndex = bisect.bisect(region, End)
    if startIndex == 1 or endIndex == 1 or (startIndex == 0 and endIndex == 2) or (startIndex == 2 and endIndex == 0):
        is_overlap = True
    else:
        is_overlap = False
    return is_overlap

def overlap_chrom_regions(Chr, Start, End, ChrRegions):
    Start = int(Start)
    End = int(End)
    is_exist = False
    if Chr in ChrRegions:
        regions = ChrRegions[Chr] 
        for region in regions:
            is_overlap = overlap_region(Start, End, region)
            if is_overlap == True:
                is_exist = True
                break
            else:
                continue
    return is_exist
#####################################################################################

def parse_bed_regions(region_file):
    """
    1   121700000   125100000   acen
    10  38000000    41600000    acen
    """
    ChrRegions = collections.defaultdict(list)

    in_h = open(region_file, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        Chr, Start, End = lines[:3]
        Start = int(Start)
        End = int(End)
        ChrRegions[Chr].append([Start, End])
    in_h.close()
    return ChrRegions



def parse_genotype_format(Format, Value, target):
    """
    Format = "GT:AD:DP:GQ:PL"
    Value = "1/1:0,2:2:6:90"
    target = "AD"
    return '0,2'

    target = "AD,PL"
    ['0,2', '90']
    """
    Formats = Format.split(":")
    Values = Value.split(":")
    if len(Formats) != len(Values):
        print("The number of ids and values is not identical for %s and %s." % (Format, Value))
        return None
    else:
        targets = target.split(",")
        targetValues = []
        for t in targets:
            try:
                valueIndex = Formats.index(t)
            except ValueError:
                print("Please check whether the target %s in %s." % (t, Format))
                sys.exit(1)
            v = Values[valueIndex]
            targetValues.append(v)
        ### return the values
        if len(targetValues) == 1:
            return targetValues[0]
        else:
            return targetValues

def parse_genotype_format_mul(Format, Value, target):
    """
    Format = "GT:AD:DP:GQ:PL"
    Value = "1/1:0,2:2:6:90"
    target = "AD"
    return '0,2'

    target = "AD,PL"
    ['0,2', '90']
    """
    if ";" in Value:
        Value = Value.split(";")[0]

    
    Formats = Format.split(":")
    Values = Value.split(":")
    if len(Formats) != len(Values):
        print("The number of ids and values is not identical for %s and %s." % (Format, Value))
        return None
    else:
        targets = target.split(",")
        targetValues = []
        for t in targets:
            try:
                valueIndex = Formats.index(t)
            except ValueError:
                print("Please check whether the target %s in %s." % (t, Format))
                sys.exit(1)
            v = Values[valueIndex]
            targetValues.append(v)
        ### return the values
        if len(targetValues) == 1:
            return targetValues[0]
        else:
            return targetValues



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



def Infor_substution_value(Infor, targetTag, targetValue):
    TagValues = {}
    Infors = Infor.split(";")

    Tags = []
    for i in Infors:
        tag, value = i.split("=")
        Tags.append(tag)
        TagValues[tag] = value

    if targetTag in Tags:
        TagValues[targetTag] = targetValue
    else:
        print("Please check whether the target tag %s is in the information record %s." % (targetTag, Infor))
        sys.exit(1)

    Record = []
    for t in Tags:
        v = TagValues[t]
        record = t + "=" + str(v)
        Record.append(record)
    newInfor = ";".join(Record)
    return newInfor

#!/usr/bin/python
import argparse
import collections
import sys
import os
import numpy

#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/OverlapCategory.py --category /home/wuzhikun/Project/Population/population/Category/Sample_SV_category_tag.xls --bed /home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_DEL_overlap_filt.bed,/home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_INS_overlap_filt.bed,/home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_INV_overlap_filt.bed --out /home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_SV_overlap_category_stats.xls > /home/wuzhikun/Project/Population/log/categoryTag2.log 2>&1

def category_tags(cat_file):
    TagCategory = {}
    cat_h = open(cat_file, "r")
    header = cat_h.readline().strip()

    for line in cat_h:
        lines = line.strip().split("\t")
        cat = lines[0]
        tags = lines[1].split(",")
        for t in tags:
            TagCategory[t] = cat
    cat_h.close()
    return TagCategory




def bed_tag_category(TagCategory, bed_file, TypeCateNumber, TypeCateLength, Type, Categories):

    bed_h = open(bed_file, "r")
    ### get uniq tags
    UniqTags = set()

    for line in bed_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        if tag in UniqTags:
            continue
        else:
            UniqTags.add(tag)
            if tag in TagCategory:
                cat = TagCategory[tag]
                Categories.add(cat)
                ### number
                TypeCateNumber[Type][cat] += 1

                ###length
                length = tag.split("-")[2]
                length = int(float(length))
                TypeCateLength[Type][cat] += length
    bed_h.close()


def target_type(file, TYPES):
    file = file.lower()
    TYPES = [t.lower() for t in TYPES]
    targets = []
    for t in TYPES:
        if t in file:
            targets.append(t)
    if len(targets) == 1:
        target = targets[0]
    else:
        print("Please make sure there is only one item of list %s is in file names %s." % (TYPES, file))

    return target


def multiple_bed_tag_category(cat_file, bedFiles, out_file):
    """
    cat_file:
    Category        Tags
    common  10_101728712-10_101729253-540.3-INS,10_116505797-10_116506120-321.9-DEL,10_124330938-10_124331032-91.7-INS,10_126417903-10_12641
    major   10_100093578-10_100093710-119.6-INS,10_100315059-10_100315195-145.0-INS,10_1010085-10_1010162-110.5-DEL,10_101288328-10_10128903
    poly    10_100277480-10_100277533-56.0-DEL,10_100277902-10_100278030-128.3-DEL,10_100314965-10_100315036-73.2-DEL,10_10033475-10_1003352
    single  10_100093517-10_100093595-78.0-DEL,10_100276307-10_100276359-52.0-INS,10_100276329-10_100276418-89.0-DEL,10_100277884-10_100278

    bed_file:
    1   100922995   100923315   1_100922995-1_100923315-319.7-DEL   1   100922995   100923315   CHM1_chr1-100922996-DEL-320 DEL 320
    1   101253202   101253276   1_101253202-1_101253276-72.8-DEL    1   101253197   101253252   CHM1_chr1-101253198-DEL-55  DEL 55
    1   101258573   101258649   1_101258573-1_101258649-74.8-DEL    1   101258572   101258648   CHM1_chr1-101258573-DEL-76  DEL 76

    out_file:
    Type    Category        Number  Length
    del     common  258     115560
    del     major   4856    2707766
    del     poly    9593    5221666
    del     single  2361    798011
    """

    TYPES = ["DEL", "INS", "DUP", "INV"]
    TYPES = [t.lower() for t in TYPES]

    TagCategory = category_tags(cat_file)


    TypeCateNumber = collections.defaultdict(lambda: collections.Counter())
    TypeCateLength = collections.defaultdict(lambda: collections.Counter())

    Categories = set()

    bedfiles = bedFiles.split(",")
    for f in bedfiles:
        f = f.strip()
        Type = target_type(f, TYPES)
        bed_tag_category(TagCategory, f, TypeCateNumber, TypeCateLength, Type, Categories)

    ### out number
    out_h = open(out_file, "w")
    out_h.write("Type\tCategory\tNumber\tLength\n")

    CatList = sorted(list(Categories))
    print(CatList)
    print(TypeCateNumber)
    for t in TYPES:
        if t in TypeCateNumber:
            Numbers = []
            Lengths = []
            for c in CatList:
                if c in TypeCateNumber[t]:
                    number = TypeCateNumber[t][c]
                    length = TypeCateLength[t][c]
                else:
                    number = 0
                    length = 0
                Numbers.append(str(number))
                Lengths.append(str(length))
            ### output number for type and category
            for cc, nn, ll in zip(CatList, Numbers, Lengths):
                t = t.upper()
                cc = cc.capitalize()
                out_h.write("%s\t%s\t%s\t%s\n" % (t, cc, nn, ll))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Statistics for overlapped SV of category and type .")
    parser.add_argument("-c", "--category", help="The input category file.")
    parser.add_argument("-b", "--bed", help="The input bed file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    multiple_bed_tag_category(args.category, args.bed, args.out)

if __name__ == "__main__":
    main()















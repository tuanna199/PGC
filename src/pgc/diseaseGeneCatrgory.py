#!/usr/bin/python
import argparse
import os
import sys
import collections

#usage: python ~/github/NanoHub/src/NanoHub/diseaseGeneCatrgory.py --annotation OMIM_Disease_gene_tag_genotypes.txt,GWAS_gene_tag_genotypes.txt --name OMIM,GWAS --number 405 --out temp.txt

def gene_category(category_file, gene_anno, total_number):
    """
    gene_anno:
    ADD1    hypertension    4_2337064-4_3678059-1340995-DEL 1       CN279,0/1
    ADORA2A adenocarcinoma,adenoma,adrenoleukodystrophy,aids,albinism,alopecia,alzheimer disease,amelogenesis imperfecta,amyloidosis,aneurysm,anomalies,arr
    AGT     hypertension    1_230703950-1_230706819-2869-DEL        1       CN277,0/1
    """
    TagCategory, CategoryC = category_tags(category_file)

    CategoryCount = collections.Counter()
    gene_h = open(gene_anno, "r")
    for line in gene_h:
        lines = line.strip().split("\t")
        gene = lines[0]
        tag = lines[2]
        # number = lines[3].split()[0]
        # number = int(number)
        # ratio = number / total_number
        # if number == 1:
        #     cat = "Singleton"
        # elif number != 1 and ratio<=0.01:
        #     cat = "Rare"
        # elif ratio > 0.01 and ratio <= 0.05:
        #     cat = "Low"
        # elif ratio > 0.05:
        #     cat = "Common"
        cat = TagCategory[tag]
        CategoryCount[cat] += 1
    gene_h.close()

    return CategoryCount


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
            CategoryCount["All"] += 1
    in_h.close()
    return TagCategory, CategoryCount
    



def multiple_disease_category(anno_file, anno_name, category_file, out_file,  total_number):
    """
    out_file:
    Database        Category        Number
    GWAS    Singleton       637
    GWAS    Rare    170
    GWAS    Low     138
    GWAS    Common  169
    OMIM    Singleton       44
    OMIM    Rare    12
    OMIM    Low     9
    OMIM    Common  18
    """
    Categories = ["Singleton", "Rare", "Low", "Common"]
    total_number = int(total_number)
    if total_number == 0:
        print("Please make sure the total number is not zero.")
        sys.exit(1)

    anno_files = anno_file.split(",")
    anno_files = [f.strip() for f in anno_files]

    anno_names = anno_name.split(",")
    anno_names = [n.strip() for n in anno_names]

    out_h = open(out_file, "w")

    diseasCatCount = {}
    if len(anno_files) == len(anno_names):
        for f, n in zip(anno_files, anno_names):
            CategoryCount = gene_category(category_file, f, total_number)
            diseasCatCount[n] = CategoryCount
    else:
        print("Please check whether the length of annotation files and names is identical.")
        sys.exit(1)

    out_h.write("Database\tCategory\tNumber\n")
    dd = sorted(list(diseasCatCount.keys()))
    for d in dd:
        catCount = diseasCatCount[d]
        for c in Categories:
            if c in catCount:
                n = catCount[c]
            else:
                n = 0
            out_h.write("%s\t%s\t%d\n" % (d, c, n))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Statistics for disease database and category.")
    parser.add_argument("-a", "--annotation", help="The input annotation file")
    parser.add_argument("-n", "--name", help="The names of annotation files.")
    parser.add_argument("-c", "--category", help="The category file.")
    parser.add_argument("-s", "--number", help="The number of total samples.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    multiple_disease_category(args.annotation, args.name, args.category, args.out, args.number)



if __name__ == "__main__":
    main()





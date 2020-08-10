#!/usr/bin/python
from BaseFunc import category_tags
import argparse
import collections
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/annoCategoryRegion.py --category  Category/Sample_SV_category_tag.xls --annotation /home/wuzhikun/Project/Population/population/bed/database/Sample_all_novel_tag_overlap.bed --out temp --targetCategory "major, common"  --targetRegion "CDS"

def select_anno_category_region(category_file, anno_file, out_file, category, region):
    """
    category_file:
    10      100981149       100981693       10_100981149-10_100981693-544.0-DEL     10      100981167       100981229       SEMA4G  ENST00000210633.3       CDS
    10      100981149       100981693       10_100981149-10_100981693-544.0-DEL     10      100981229       100983304       SEMA4G  ENST00000210633.3       Intron
    10      100981149       100981693       10_100981149-10_100981693-544.0-DEL     10      100981229       100984489       SEMA4G  ENST00000517724.5       Intron
    10      100981149       100981693       10_100981149-10_100981693-544.0-DEL     10      100981375       100981462       MRPL43  ENST00000299179.9       UTR5
    10      100981149       100981693       10_100981149-10_100981693-544.0-DEL     10      100981376       100981490       MRPL43  ENST00000342071.5       UTR5
    """
    Category = category.split(",")
    Category = [c.strip().lower() for c in Category]
    print(Category)

    Region = region.split(",")
    Region = [r.strip().lower() for r in Region]
    print(Region)

    ### get the tag and corresponding category
    TagCategory, CategoryCount = category_tags(category_file)
    
    AllCategory = list(CategoryCount.keys())
    for c in Category:
        if c not in AllCategory:
            print("Please check whether the target category %s is in file %s." % (c, category_file))
            sys.exit(1)


    TargetGenes = set()

    anno_h = open(anno_file, "r")
    for line in anno_h:
        line = line.strip()
        lines = line.split("\t")
        tag = lines[3]
        gene = lines[7]
        region = lines[9]
        region = region.lower()

        ### get category for tag
        try:
            Cat = TagCategory[tag]
            Cat = Cat.lower()
        except KeyError:
            print("Please check whether the tag %s is in category file %s." % (tag, category_file))
            sys.exit(1)

        if Cat in Category and region in Region:
            TargetGenes.add(gene)

    ### output the file
    out_h = open(out_file, "w")
    Genes = sorted(TargetGenes)
    for g in Genes:
        out_h.write("%s\n" % g)
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the annotation gene for target category and region of genes.")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    parser.add_argument("-c", "--category", help="The input category file.")
    parser.add_argument("-t", "--targetCategory", default="major,common", help="The target category.")
    parser.add_argument("-r", "--targetRegion", default="cds", help="The target region of genes.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    select_anno_category_region(args.category, args.annotation, args.out, args.targetCategory, args.targetRegion)

if __name__ == "__main__":
    main()




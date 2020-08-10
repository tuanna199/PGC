#!/usr/bin/python
import argparse
import collections
import os


#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/AnnoCategory.py --category /home/wuzhikun/Project/Population/population/Category/Sample_SV_category_tag.xls --annotation /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_modify.tsv --out /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_anno

def annotation_category(cate_file, anno_file, outPrefix):
    """
    cate_file:
    Category        Tags
    common  10_101728712-10_101729253-540.3-INS,10_116505797-10_116506120-321.9-DEL,10_124330938-10_124331032-91.7-INS,10_126417903-10_12641
    major   10_100093578-10_100093710-119.6-INS,10_100315059-10_100315195-145.0-INS,10_1010085-10_1010162-110.5-DEL,10_101288328-10_10128903
    poly    10_100277480-10_100277533-56.0-DEL,10_100277902-10_100278030-128.3-DEL,10_100314965-10_100315036-73.2-DEL,10_10033475-10_1003352
    single  10_100093517-10_100093595-78.0-DEL,10_100276307-10_100276359-52.0-INS,10_100276329-10_100276418-89.0-DEL,10_100277884-10_100278
    """
    ### get dict of catogery
    TagCategory = {}
    Categories = []
    cat_h = open(cate_file, "r")
    header = cat_h.readline().strip()
    for line in cat_h:
        lines = line.strip().split("\t")
        cat = lines[0]
        Categories.append(cat)
        tags = lines[1].split(",")
        for t in tags:
            TagCategory[t] = cat

    ### dict of annotation
    CategoryRecord = collections.defaultdict(list)
    anno_h = open(anno_file, "r")
    annoHeader = anno_h.readline().strip()
    for line in anno_h:
        line = line.strip()
        lines = line.split("\t")
        Tag = lines[0]
        try:
            category = TagCategory[Tag]
        except KeyError:
            print("Please check whether there is tag %s in category file %s." % (Tag, cate_file))
            sys.exit(1)
        CategoryRecord[category].append(line)
    anno_h.close()


    ### output the file
    outdir = "/".join(outPrefix.split("/")[:-1])
    outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for c in Categories:
        ch = open("%s_%s.xls" % (outPrefix, c.capitalize()), "w")
        ch.write("%s\n" % annoHeader)
        Records = CategoryRecord[c]
        for r in Records:
            ch.write("%s\n" % r)
        ch.close()


def main():
    parser = argparse.ArgumentParser(description="Categories for annotation file.")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    parser.add_argument("-c", "--category", help="The input categories.")
    parser.add_argument("-o", "--out", help="The output prefix.")
    args = parser.parse_args()
    annotation_category(args.category, args.annotation, args.out)


if __name__ == "__main__":
    main()

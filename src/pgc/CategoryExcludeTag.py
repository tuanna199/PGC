#!/usr/bin/python
from BaseFunc import category_tags
import collections
import argparse


def exclude_tags(tag_file):
    TAGS = {}
    tag_h = open(tag_file, "r")
    for line in tag_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        TAGS[tag] = 1
    tag_h.close()
    return TAGS


def bed_exclude_tag_record(tag_file, category_file, bed_file, out_file, target):
    
    TAGS = exclude_tags(tag_file)

    TagCategory, CategoryCount = category_tags(category_file)


    ### targets: common,major
    targets = target.split(",")
    targets = [t.strip() for t in targets]

    print(targets)

    out_h = open(out_file, "w")
    bed_h = open(bed_file, "r")
    for line in bed_h:
        line = line.strip()
        lines = line.split("\t")
        Tag = lines[3]
        if Tag in TAGS:
            continue
        else:
            category = TagCategory[Tag]
            print(category)
            if category in targets:
                out_h.write("%s\n" % line)
    bed_h.close()
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the gene annotation for translocation.")
    parser.add_argument("-t", "--tag", help="The input tag file.")
    parser.add_argument("-b", "--bed", help="The input bed file.")
    parser.add_argument("-c", "--category", help="The input category file.")
    parser.add_argument("-g", "--target", default="common,major", help="The target category name.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    bed_exclude_tag_record(args.tag, args.category, args.bed, args.out, args.target)

if __name__ == "__main__":
    main()

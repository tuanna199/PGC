#!/usr/bin/python
import argparse
import sys
import collections


#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/categoryTag.py --files /home/wuzhikun/Project/Population/population/Category/Sample_SV_record_Common_SV_tag.xls,/home/wuzhikun/Project/Population/population/Category/Sample_SV_record_Major_SV_tag.xls,/home/wuzhikun/Project/Population/population/Category/Sample_SV_record_Poly_SV_tag.xls,/home/wuzhikun/Project/Population/population/Category/Sample_SV_record_Single_SV_tag.xls --output /home/wuzhikun/Project/Population/population/Category/Sample_SV_category_tag.xls

Categories = ["common", "major", "poly", "single"]

def tag_category(tag_file, CategoryTag):
    tag_h = open(tag_file, "r")
    header = tag_h.readline().strip()

    tag_file = tag_file.split("/")[-1].lower()

    ### get the category for file
    targets = []
    for c in Categories:
        if c in tag_file:
            targets.append(c)

    if len(targets) == 1:
        categ = targets[0]
    else:
        print("Please check whether the item in list %s is and only one in file %s." % (Categories, tag_file))
        sys.exit(1)


    for line in tag_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        CategoryTag[categ].append(tag)
    tag_h.close()


def multiple_file_category(fileStr, out_file):
    CategoryTag = collections.defaultdict(list)

    out_h = open(out_file, "w")
    out_h.write("Category\tTags\n")

    files = fileStr.split(",")
    for f in files:
        f = f.strip()
        tag_category(f, CategoryTag)

    for c in Categories:
        tags = CategoryTag[c]
        sortTags = ",".join(sorted(tags))
        out_h.write("%s\t%s\n" % (c, sortTags))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the tags for different categories.")
    parser.add_argument("-f", "--files", help="The input files which are separated with ','.")
    parser.add_argument("-o", "--output", help="The output file.")
    args = parser.parse_args()
    multiple_file_category(args.files, args.output)

if __name__ == "__main__":
    main()



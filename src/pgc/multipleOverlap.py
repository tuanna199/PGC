#!/usr/bin/python
import sys
import argparse
import collections
import os

def multiple_file_tags(fileStr, out_file):
    AllTags = set()
    FileTags = collections.defaultdict(list)
    files = fileStr.split(",")
    files = [f.strip() for f in files]
    for f in files:
        in_h = open(f, "r")
        for line in in_h:
            lines = line.strip().split("\t")
            tag = lines[0]
            AllTags.add(tag)
            FileTags[f].append(tag)
        in_h.close()

    out_h = open(out_file, "w")
    sortedTags = sorted(list(AllTags))
    sortedFiles = sorted(list(FileTags.keys()))
    out_h.write("%s\n" % "\t".join(sortedFiles))
    for t in sortedTags:
        values = []
        for f in sortedFiles:
            if t in FileTags[f]:
                v = t
            else:
                v = "NA"
            values.append(v)
        out_h.write("%s\n" % "\t".join(values))


def main():
    parser = argparse.ArgumentParser(description="Tag matrix for multiple files.")
    parser.add_argument("-f", "--files", help="The input files.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    multiple_file_tags(args.files, args.out)



if __name__ == "__main__":
    main()


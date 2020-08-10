#1/usr/bin/python
import collections
import argparse
import sys
import os
import random

#usage: python ~/github/NanoHub/src/NanoHub/shuffleGroup.py --input geno_group.txt --out temp.txt


def shuffle_sample_group(group_file, out_file):
    """
    group_file:
    geno    group
    CN001   South
    CN002   South
    """
    Samples = []
    Groups = []
    group_h = open(group_file, "r")
    header = group_h.readline().strip()
    for line in group_h:
        lines = line.strip().split("\t")
        sample, group = lines
        Samples.append(sample)
        Groups.append(group)
    group_h.close()

    ### shuffle the groups
    random.shuffle(Groups)


    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)
    if len(Samples) == len(Groups):
        for s, g in zip(Samples, Groups):
            out_h.write("%s\t%s\n" % (s, g))
    else:
        print("Please check whether that the length of samples and groups is identical.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Shuffle the group name for sample")
    parser.add_argument("-i", "--input", help="The input file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    shuffle_sample_group(args.input, args.out)



if __name__ == "__main__":
    main()

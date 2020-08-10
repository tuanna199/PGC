#!/usr/bin/python
from BaseFunc import column_index
import collections
import sys
import argparse

#usage: python ~/github/NanoHub/src/NanoHub/PCAaddGroup.py --structure ../admixture/Sample_assign_population.xls --pca Sample_PCA-1.txt --out temp.xls



def sample_group(group_file, column):
    """
    geno    P1      P2      group
    CN001   0.852993        0.147007        P1
    CN002   0.921790        0.078210        P1
    CN003   0.913900        0.086100        P1
    CN004   0.962459        0.037541        P1
    """
    SampleGroup = {}

    in_h = open(group_file, "r")
    headers = in_h.readline().strip().split("\t")
    groupIndex = column_index(headers, column)

    for line in in_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        group = lines[groupIndex]
        SampleGroup[sample] = group
    in_h.close()
    return SampleGroup


def add_group(group_file, pca_file, out_file, column):
    """
    PCA_file:
    FID     IID     PCA1    PCA2
    CN001   CN001   -0.0405 0.0088
    CN002   CN002   -0.0529 0.0058
    CN003   CN003   -0.0481 0.0094

    out_file:
    FID     IID     PCA1    PCA2    Group
    CN001   CN001   -0.0405 0.0088  P1
    CN002   CN002   -0.0529 0.0058  P1
    CN003   CN003   -0.0481 0.0094  P1
    """
    SampleGroup = sample_group(group_file, column)

    in_h = open(pca_file, "r")
    header = in_h.readline().strip()

    out_h = open(out_file, "w")
    out_h.write("%s\tGroup\n" % header)

    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        sample = lines[0]
        if sample in SampleGroup:
            group = SampleGroup[sample]
            out_h.write("%s\t%s\n" % (line, group))
        else:
            print("Please check whether there is group name for sample %s in file %s." % (sample, group_file))
            sys.exit(1)
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="The the count number of sliding window across the genome.")
    parser.add_argument("-s", "--structure", help="The input structure file containing sample and group.")
    parser.add_argument("-c", "--column", help="The target column of group file.")
    parser.add_argument("-p", "--pca", help="The pca file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    add_group(args.structure, args.pca, args.out, args.column)

if __name__ == "__main__":
    main()





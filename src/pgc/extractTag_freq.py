#!/usr/bin/python
import argparse
import collections
import os
import sys

def bed_target_tag(bed_file):
    TargetTags = []

    bed_h = open(bed_file, "r")
    for line in bed_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        TargetTags.append(tag)
    bed_h.close()

    return TargetTags


def extract_tag_SV_frequency(bed_file, tag_freq, out_file):
    """
    bed_file:


    tag_freq:
    Tag Ref_freq    Alt_freq    MAF
    1_66288-1_66527-239.0-DEL   0.9974  0.0026  0.0026
    1_67910-1_68341-431.0-DEL   0.9974  0.0026  0.0026
    1_83968-1_84057-89.0-INS    0.9974  0.0026  0.0026

    # Chr1    Pos1    Chr2    Pos2    SVlength    Type    Ref_freq    Alt_freq
    # 1   66288   1   66527   239.0   DEL 0.997   0.003
    # 1   67910   1   68341   431.0   DEL 0.997   0.003
    # 1   83968   1   84057   89.0    INS 0.997   0.003

    """
    TargetTags = bed_target_tag(bed_file)

    in_h = open(tag_freq, "r")
    header = in_h.readline().strip()

    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        # Chr1, Pos1, Chr2, Pos2, Length, SVType, RefFreq, AltFreq = lines
        # Tag = "%s_%s-%s_%s-%s-%s" % (Chr1, Pos1, Chr2, Pos2, Length, SVType)
        Tag = lines[0]
        if Tag in TargetTags:
            out_h.write("%s\n" % line)
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the target SV tag and frequency.")
    parser.add_argument("-b", "--bed", help="The input bed file.")
    parser.add_argument("-t", "--freq", help="The input tag and frequency file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    extract_tag_SV_frequency(args.bed, args.freq, args.out)

if __name__ == "__main__":
    main()



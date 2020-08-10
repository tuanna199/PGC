#!/usr/bin/python
import collections
import os
import sys
import argparse

def tag_record(bed_file):
    TagRecord = {}
    bed_h = open(bed_file, "r")
    for line in bed_h:
        line = line.strip()
        lines = line.split("\t")
        tag = lines[3]
        TagRecord[tag] = line
    bed_h.close()

    return TagRecord


def sig_loci_tag_annotation(bed_file, annotation_file, sig_loci, bed_out, anno_out, tagColumn):
    """
    bed_file:
    1   10059   13301   1_10059-1_10060-3242-INS
    1   10380   10701   1_10380-1_10381-321-INS
    1   10657   10735   1_10657-1_10760-78-INS

    annotation file:
    1   67910   68341   1_67910-1_68341-433-DEL 1   68090   69090   OR4F5   ENST00000335137.3   UpStream
    1   88684   88831   1_88684-1_88831-147-DEL 1   88294   89294   RP11-34P13.7    ENST00000466430.5   DownStream
    1   88684   88831   1_88684-1_88831-147-DEL 1   88550   89550   RP11-34P13.8    ENST00000495576.1   DownStream
    """
    tagColumn = int(tagColumn)
    tagIndex = tagColumn - 1

    TagRecord = tag_record(bed_file)
    AnnoRecord = tag_record(annotation_file)

    sig_h = open(sig_loci, "r")
    bed_oh = open(bed_out, "w")
    anno_oh = open(anno_out, "w")
    for line in sig_h:
        lines = line.strip().split("\t")
        tag = lines[tagIndex]
        if tag in TagRecord:
            region = TagRecord[tag]
            bed_oh.write("%s\n" % region)
        else:
            print("Please check whetehr the tag %s is in file %s and %s." % (tag, bed_file, annotation_file))
            sys.exit(1)

        if tag in AnnoRecord:
            annotation = AnnoRecord[tag]
            anno_oh.write("%s\t%s\n" % (tag, annotation))
    sig_h.close()
    bed_oh.close()
    anno_oh.close()



def main():
    parser = argparse.ArgumentParser(description='Get the annotation record for significants of GWAS.')
    parser.add_argument('-b', '--bed', help='The input bed file.')
    parser.add_argument('-a', '--annotation', help='The input annotation file.')
    parser.add_argument('-s', '--significant', help='The input significant loci file.')
    parser.add_argument("-bo", "--bedOut", help="The output file for region record.")
    parser.add_argument("-ao", "--annoOut", help="The output file for annotation record.")
    parser.add_argument("-c", "--tagColumn", default = 2, help="The colum number for the tags.")
    args = parser.parse_args()
    sig_loci_tag_annotation(args.bed, args.annotation, args.significant, args.bedOut, args.annoOut, args.tagColumn)

    
if __name__ == '__main__':
    main()




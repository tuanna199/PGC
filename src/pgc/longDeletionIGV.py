#!/usr/bin/python
import collections
import os
import argparse
import sys

#usage: python ~/github/NanoHub/src/NanoHub/longDeletionIGV.py --tag /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion.txt --out /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion_igv.batch --IGVDir /home/wuzhikun/Project/Population/IGVRegion/LongDEL

def long_deletion_IGV(deletion_file, out_file, igv_dir):
    """
    deletion_file:
    1_13249651-1_13395879-146232-DEL        CN012,CN044,CN079,CN099,CN113,CN119,CN148,CN162,CN180,CN214,CN216,CN262,CN288,CN333,CN338,M426-1,M561-1,M579-2
    1_54504649-1_54752462-247813-DEL        M489-1
    1_61814548-1_61980438-165890-DEL        CN139
    """
    del_h = open(deletion_file, "r")
    out_h = open(out_file, "w")
    for line in del_h:
        lines = line.strip().split("\t")
        Tag, Sample = lines[:2]
        samples = Sample.split(",")
        ### just select the sample startswith "CN"
        if len(samples) == 1 and samples[0].startswith("M"):
            continue

        tags = Tag.split("-")
        SVlength = int(tags[-2])
        Chr, Start = tags[0].split("_")
        End = tags[1].split("_")[-1]
        Start = int(Start)
        End = int(End)

        extendDistance = SVlength * 0.5
        newStart = Start - extendDistance
        if newStart < 0:
            newStart = 0
        newEnd = End + extendDistance

        selectSample = samples[0]
        out_h.write("new\n")
        out_h.write("genome /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa\n")
        out_h.write("snapshotDirectory %s\n" % igv_dir)
        out_h.write("load /home/wuzhikun/Project/NewChinese/mapping/minimap2/%s.bam\n" % selectSample)
        out_h.write("maxPanelHeight 700\n")
        # out_h.write("group SUPPLEMENTARY\n")
        out_h.write("goto %s:%d-%d\n" % (Chr, newStart, newEnd))
        out_h.write("snapshot %s\n\n" % Tag)
    out_h.write("exit")
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the igv batch file for long deletion.")
    parser.add_argument("-t", "--tag", help="The input tag file.")
    parser.add_argument("-o", "--out", help="The output batch file.")
    parser.add_argument("-d", "--IGVDir", help="The IGV directory of the picture files..")
    args = parser.parse_args()
    long_deletion_IGV(args.tag, args.out, args.IGVDir)




if __name__ == "__main__":
    main()

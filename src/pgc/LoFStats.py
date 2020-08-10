#!/usr/bin/python
from __future__ import division
import argparse
import collections
import sys
import os
import statistics
import numpy


#usage: python ~/github/NanoHub/src/NanoHub/LoFStats.py --matrix /home/wuzhikun/Project/Population/population/Merge/Sample_common_SV_Tag.xls --tag /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed --out temp.xls


def lof_tag_stats(tag_file, tag_matrix, out_file, name, SV_sample):    
    """
    tag_file:
    1       939439  939444  1_939439-1_939575-122-INS       1       939271  939460  SAMD11  ENST00000622503.4       CDS
    1       939439  939444  1_939439-1_939575-122-INS       1       939274  939460  SAMD11  ENST00000342066.7       CDS


    tag_matrix:
    Tag     CN001   CN002   CN003   CN004   CN005   CN006   CN007   CN008   CN009   CN010   CN011   CN012   CN013   CN014   CN015   CN016   CN017   CN018
    1_10059-1_10060-3242-INS        -       -       -       -       -       -       -       -       -       -       -       -       -       -       -
    """
    tagSets = set()
    TagLength = {}

    GeneTags = collections.defaultdict(set)
    TagGenes = collections.defaultdict(set)

    LongDEL = set()

    tag_h = open(tag_file, "r")
    for line in tag_h:
        lines = line.strip().split("\t")
        if lines[0] != "":
            tag = lines[3]
            gene = lines[7]
            tagSets.add(tag)
            tags = tag.split("-")
            SVlength, SVtype = tags[-2:]
            SVlength = int(SVlength)
            TagLength[tag] = SVlength

            GeneTags[gene].add(tag)
            TagGenes[tag].add(gene)

            ### get the long deletion
            # if SVtype == "DEL" and SVlength > 100000:
            LongDEL.add(tag)
    tag_h.close()

    SVNumber = len(tagSets)
    LengthList = []
    for a in TagLength:
        length = int(TagLength[a])
        LengthList.append(length)
    print(LengthList)
    medianLength = statistics.median(LengthList)

    ### SV per gene
    geneSVNum = []
    for g in GeneTags:
        tag = GeneTags[g]
        tagNum = len(tag)
        geneSVNum.append(tagNum)

    geneAVESV = numpy.mean(geneSVNum)


    ### genes per SV
    SVGeneNum = []
    for t in TagGenes:
        genes = TagGenes[t]
        geneNum = len(genes)
        SVGeneNum.append(geneNum)
    SVGeneNumAVE = numpy.mean(SVGeneNum)

    medianLength = "%.1f" % (medianLength / 1000)
    geneAVESV = "%.1f" % geneAVESV
    SVGeneNumAVE = "%.1f" % SVGeneNumAVE
    #print(SVNumber)
    #print(medianLength)
    #print(geneAVESV)
    #print(SVGeneNumAVE)

    matrix_h = open(tag_matrix, "r")
    sample_h = open(SV_sample, "w")
    headers = matrix_h.readline().strip().split("\t")
    samples = headers[1:]
    allValues = []
    TagValueNum = []
    for line in matrix_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        values = lines[1:]
        if tag in tagSets:
            newValues = [1 if v != "-" else 0 for v in values]
            allValues.append(newValues)

            ### sum for value of each tag
            valueSum = sum(newValues)
            TagValueNum.append(valueSum)

        ### output the sample of long deletion
        if tag in LongDEL:
            targetSample = []
            for i in range(len(values)):
                value = values[i]
                if value != "-":
                    sample = samples[i]
                    targetSample.append(sample)
            sample_h.write("%s\t%s\n" % (tag, ",".join(targetSample)))

    matrix_h.close()
    sample_h.close()



    ### transpose two-dim array
    ### SV number per individual
    transMatrix = [*zip(*allValues)]
    IndividualSVs = [sum(v)  for v in transMatrix]
    individualAverage = numpy.mean(IndividualSVs)

    individualAverage = "%.1f" % individualAverage

    out_h = open(out_file, "w")
    out_h.write("Type\tNumber\tMedianLength\tSVsPerGene\tGenePerSV\tSVsPerIndividual\n")
    out_h.write("%s\t%d\t%s\t%s\t%s\t%s\n" % (name, SVNumber, medianLength, geneAVESV, SVGeneNumAVE, individualAverage))



def main():
    parser = argparse.ArgumentParser(description="Statistics for SVs of individual.")
    parser.add_argument("-t", "--tag", help="The input tag file.")
    parser.add_argument("-m", "--matrix", help="The matrix of tag in individuals")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-n", "--name", default="LoF", help="The name of the database.")
    parser.add_argument("-s", "--sample", help="The output file contain long deletion and the samples.")
    args = parser.parse_args()
    lof_tag_stats(args.tag, args.matrix, args.out, args.name, args.sample)



if __name__ == "__main__":
    main()



#!/usr/bin/python
#-*- coding:utf-8 -*-
import pysam
import re
import os
import argparse
import collections
from tinyfasta import FastaParser
import sys


#usage: python ~/github/NanoHub/src/NanoHub/BAMAlignValue.py --bam M605-2_15_minimap2.bam --out M605-2_15_unmapped_target.fatsa --sequence ../SubReads/M605-2_15/Assembly_split.fasta --query M605-2_15_unmapped_target.query


def cigar_values(CIGAR):
    """
    CIGAR:
    [(0, 9), (2, 2), (0, 10), (2, 1), (0, 7), (2, 5), (0, 12), (1, 1), (0, 16), (2, 1), (0, 7), (2, 3), (0, 1), (2, 3), (0, 31), (1, 1), (0, 19), (2, 1), (0, 1), (2, 1), (0, 10), (2, 4), (0, 46), (2, 2), (0, 31), (2, 1), (0, 5), (1, 1), (0, 13), (2, 3), (0, 13), (2, 1), (0, 14), (1, 1), (0, 5), (1, 3), (0, 10), (1, 2), (0, 29), (2, 1), (0, 4), (2, 1), (0, 8), (2, 1), (0, 25), (2, 1), (0, 2), (2, 1), (0, 15), (2, 1), (0, 1), (2, 2), (0, 24), (2, 1), (0, 11), (2, 1), (0, 2), (1, 1), (0, 47), (2, 1), (0, 8), (2, 1), (0, 34), (2, 1), (0, 3), (2, 2), (0, 13), (2, 1), (0, 8), (2, 1), (0, 1), (2, 1), (0, 16), (2, 1), (0, 64), (2, 1), (0, 5), (2, 1), (0, 43), (2, 1), (0, 8), (2, 1), (0, 46), (2, 1), (0, 25), (2, 1), (0, 7), (2, 1), (0, 15), (2, 1), (0, 7), (2, 1), (0, 7), (2, 1), (0, 4), (2, 1), (0, 1), (2, 1), (0, 11), (2, 1), (0, 25), (2, 1), (0, 32), (2, 1), (0, 2), (2, 1), (0, 15), (2, 1), (0, 8), (2, 1), (0, 16), (1, 2), (0, 18), (2, 1), (0, 6), (2, 2), (0, 1), (2, 1), (0, 14), (2, 1), (0, 15), (2, 1), (0, 1), (2, 1), (0, 9), (2, 1), (0, 28), (4, 44)]

    return:
    Counter({0: 944, 2: 73, 4: 44, 1: 12})
    """
    CigarValus = collections.Counter()
    for i in CIGAR:
        t, v = i
        v = int(v)
        CigarValus[t] += v
    return CigarValus


def select_cigar_values(cigar, values):
    valueIndex = values.split(",")
    valueIndex = [int(v.strip()) for v in valueIndex]
    selectValues = []
    for i, v in enumerate(cigar):
        if i in valueIndex:
            selectValues.append(v)
    selectValues = [int(v) for v in selectValues]
    return selectValues

def unmapped_flanking_query(UnmappedQuery, QUERY, number):
    """
    number: extended number
    """
    number = int(number)
    FlankQuery = []
    for query in UnmappedQuery:
        tag, value = query.split("-")
        value = int(value)
        leftValue = value - number
        rightValue = value + number
        left = "%s-%d" % (tag, leftValue)
        right = "%s-%d" % (tag, rightValue)
        
        ### filt tag out of Query because new left or right may be out of original region
        if left in QUERY:
            FlankQuery.append(left)
        if right in QUERY:
            FlankQuery.append(right)
    return FlankQuery


def split_list(IndexList, number):
    """
    number: split distance of list
    number = 2

    Index_list = [124, 125, 126, 271, 273, 274, 276, 882, 884, 885, 887]

    return:
    [124, 125, 126], [271, 273, 274, 276], [882, 884, 885, 887]
    """
    TotalSplits = []

    number = int(number)
    listNumber = len(IndexList)
    start = IndexList[0]
    Split = []
    Split.append(start)
    for i in range(1, listNumber):
        index = IndexList[i]
        if (index - start) > number:
            TotalSplits.append(Split)
            Split = []
            start = index
            Split.append(start)
        else:
            Split.append(index)
            start = index
    TotalSplits.append(Split)

    return TotalSplits


def contig_sequence(contig_file):
    ContigSeq = {}
    for record in FastaParser(contig_file):
        desc = str(record.description)
        desc = desc.lstrip(">").split()[0]
        seq = str(record.sequence)
        ContigSeq[desc] = seq

    return ContigSeq


def bam_align_value_stats(contig_file, bam_file, out_file, query_file):
    """
    The output is a list of (operation, length) tuples, such as [(0, 30)]
    
    CIGAR:
    M   BAM_CMATCH  0
    I   BAM_CINS    1
    D   BAM_CDEL    2
    N   BAM_CREF_SKIP   3
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD    6
    =   BAM_CEQUAL  7
    X   BAM_CDIFF   8
    B   BAM_CBACK   9
    NM  NM tag  10
    
    get_cigar_stats:
    [944, 12, 73, 0, 44, 0, 0, 0, 0, 0, 127]


    """
    QueryTagValues = collections.defaultdict(list)
    QueryReference  = collections.defaultdict(list)
    UnmappedQuery = []
    QUERY = set()

    out_h = open(out_file, "w")
    query_h = open(query_file, "w")

    bamRecord = pysam.AlignmentFile(bam_file,'rb')
    for line in bamRecord:
        # try:
        #     MD = line.get_tag('MD')
        #     CIGAR = line.cigar
        # except KeyError:
        #     continue

        
        ### get flag 
        flag = line.flag
        ### flag is int
        # print(isinstance(flag, int))
        # Ture

        ### get query and reference name
        query = line.query_name
        reference = line.reference_name

        QUERY.add(query)

        if reference != None:
            ### cigar value list
            cigar = line.get_cigar_stats()[0]
            NM = line.get_tag('NM')
            AS = line.get_tag("AS")

            ### DEL,INS,SOFT,NM
            MissValues = select_cigar_values(cigar, "1,2,4,10")
            TotalMiss = sum(MissValues)

            ### get total missmatch value and align score
            SelectValue = (TotalMiss, AS)
            QueryTagValues[query].append(SelectValue)

            ### query to reference
            QueryReference[query].append(reference)
        else:
            UnmappedQuery.append(query)


    ### Query to multiple reference
    MultipleQuery = []
    for q in QueryReference:
        refs = len(QueryReference[q])
        if refs > 1:
            MultipleQuery.append(q)


    ### flanking query
    UnmappedFlank = unmapped_flanking_query(UnmappedQuery, QUERY, 1)
    MultipleFlank = unmapped_flanking_query(MultipleQuery, QUERY, 1)

    # AllTargets = UnmappedQuery + UnmappedFlank + MultipleFlank + MultipleFlank

    AllTargets = UnmappedQuery + UnmappedFlank 

    ### fille the gap if query continus
    ContigIndex = collections.defaultdict(set)
    for t in AllTargets:
        contig, index = t.split("-")
        index = int(index)
        ContigIndex[contig].add(index)

    ChrTargetContigs = {}
    for c in ContigIndex:
        index = sorted(list(ContigIndex[c]))
        TotalSplits = split_list(index, 2)

        TargetList = []
        for split in TotalSplits:
            print(split)

            Target = []
            minIndex = min(split)
            maxIndex = max(split)
            for i in range(minIndex, maxIndex+1):
                newContig = "%s-%d" % (c, i)
                Target.append(newContig)
            TargetList.append(Target)

        ChrTargetContigs[c] = TargetList


    ### get the sequence for contigs
    ContigSeq = contig_sequence(contig_file)
    for Chr in ChrTargetContigs:
        contigList = ChrTargetContigs[Chr]
        ### contigList is a two dimer list

        for clust in contigList:
            Contigs = []
            Sequences = ""
            for contig in clust:

                ### output selected contig name
                query_h.write("%s\n" % contig)


                Contigs.append(contig)

                try:
                    seq = ContigSeq[contig]
                except KeyError:
                    print("Please check whether the contig name %s is in file %s." % (contig, contig_file))
                    sys.exit(1)
                Sequences += seq

            if len(Contigs) == 1:
                newContig = Contigs[0]
            else:
                newContig = "%s_%s" % (Contigs[0], Contigs[-1])

            out_h.write(">%s\n%s\n" % (newContig, Sequences))
    out_h.close()
    query_h.close()








def main():
    parser = argparse.ArgumentParser(description="Get aligned values for seq of BAM files.")
    parser.add_argument("-b", "--bam", help="The input bam file.")
    parser.add_argument("-s", "--sequence", help="The contig sequence file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-q", "--query", help="Output selected query name of sequence.")
    args = parser.parse_args()
    bam_align_value_stats(args.sequence, args.bam, args.out, args.query)


if __name__ == "__main__":
    main()


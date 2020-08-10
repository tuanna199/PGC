#!/usr/bin/python
# from tinyfasta import FastaParser
from __future__ import division
import collections
import sys
import os
import argparse
import bisect
import operator

#usage: python ~/github/NanoHub/src/NanoHub/blastFiltRecord.py --input 115.contigs_blast_refGenome.txt --out temp.xls --column  qseqid,sseqid,evalue,qlen,length,pident,mismatch,gapopen,qstart,qend,sstart,send,stitle --lengthRatio 0.75


# def fasta_title_length(fa_file):
#     SeqTitleLen = {}
#     for record in FastaParser(fa_file):
#         desc = str(record.description)
#         desc = desc.split()[0].lstrip(">")
#         seq = str(record.sequence)
#         seqLen = len(seq)
#         SeqTitleLen[desc] = seqLen
#     return SeqTitleLen

def longest_region(regions):
    """
    regions:
    (1, 128), (1, 129), (1, 149), (1, 153), (1, 154), (1, 160), (1, 164), (1, 166), (1, 167), (2, 66), (2, 160), (3, 158), (3, 160), (4, 160), (6, 160), (6, 166), (7, 160), (28, 107)]

    return:
    [(1, 167), (2, 160), (6, 166), (28, 107)]
    """
    StartList = collections.defaultdict(list)
    for region in regions:
        start, end = region
        StartList[start].append(end)

    StartEnd = []
    for start in StartList:
        ends = StartList[start]
        maxEnd = max(ends)
        StartEnd.append((start, maxEnd))

    EndList = collections.defaultdict(list)
    for region in StartEnd:
        start, end = region
        EndList[end].append(start)

    StartEnd2 = []
    for end in EndList:
        starts = EndList[end]
        minStart = min(starts)
        StartEnd2.append((minStart, end))
    return sorted(StartEnd2)

def point_index(region, points):
    slicePoints = []
    for point in points:
        leftIndex = bisect.bisect_left(region, point)
        rightIndex = bisect.bisect_right(region, point)
        if leftIndex == 1 and rightIndex == 1:
            slicePoints.append(point)
    return slicePoints

def combine_regions(regions):
    """
    regions:
    [(1, 167), (2, 160), (6, 166), (28, 107)]

    return:
    [(1, 2), (2, 6), (6, 28), (28, 107), (107, 160), (160, 166), (166, 167)]
    """
    points = set()
    for region in regions:
        start, end = region
        points.add(start)
        points.add(end)

    SliceSet = set()
    pointList = sorted(list(points))
    for region in regions:
        start, end = region

        slicePoints = point_index(region, pointList)
        
        if slicePoints == []:
            SliceSet.add((start, end))
        else:
            pointLen = len(slicePoints)
            if pointLen == 1:
                SliceSet.add((start, slicePoints[0]))
                SliceSet.add((slicePoints[0], end))
            else:
                for i in range(pointLen):
                    point = slicePoints[i]
                    if i == 0:
                        SliceSet.add((start, slicePoints[0]))
                    else:
                        SliceSet.add((slicePoints[i-1], slicePoints[i]))
                SliceSet.add((slicePoints[i], end))
    return sorted(list(SliceSet))


def once_key(SliceList):
    """
    SliceList
    [(1, 2), (2, 6), (6, 28), (28, 107), (107, 160), (160, 166), (166, 167)]

    return:
    [1, 167]    
    """
    OnceKeys = collections.Counter()
    for region in SliceList:
        start, end = region
        OnceKeys[start] += 1
        OnceKeys[end] += 1

    OnceList = []
    for o in OnceKeys:
        v = OnceKeys[o]
        if v == 1:
            OnceList.append(o)
    return sorted(OnceList)


def filt_blast_result(blast_file, columnStr, out_file, lengthRatio):
    """
    sortedRegions:
    [(1, 128), (1, 129), (1, 149), (1, 153), (1, 154), (1, 160), (1, 164), (1, 166), (1, 167), (2, 66), (2, 160), (3, 158), (3, 160), (4, 160), (6, 160), (6, 166), (7, 160), (28, 107)]

    """
    lengthRatio = float(lengthRatio)

    columns = columnStr.split(",")
    columns = [s.strip() for s in columns]
    columnLen = len(columns)

    SeqBlastRegions = collections.defaultdict(lambda: collections.defaultdict(set))
    QueryLength = {}
    SubjectLength = {}

    out_h = open(out_file, "w")

    ### get index for columns
    blast_h = open(blast_file, "r")
    firstlines = blast_h.readline().strip().split("\t")
    if len(firstlines) != columnLen:
        print("Please check whether the length of given columns %s is identical to file %s." % (columnStr, blast_file))
        sys.exit(1)

    try:
        qseqidIndex = columns.index("qseqid")
        sseqidIndex = columns.index("sseqid")
        qlenIndex = columns.index("qlen")
        qstartIndex = columns.index("qstart")
        qendIndex = columns.index("qend")
        slenIndex = columns.index("slen")
    except ValueError:
        print("Please check whether the target column %s, %s and %s is in column string %s." % ("qseqid", "qstart", "qend", columnStr))
        sys.exit(1)


    for line in blast_h:
        lines = line.strip().split("\t")
        qseqid = lines[qseqidIndex]
        sseqid = lines[sseqidIndex]

        ### just consider the record that blasted to other sequence, not itself
        if qseqid != sseqid:
            qlen = int(lines[qlenIndex])
            QueryLength[qseqid] = qlen
            slen = int(lines[slenIndex])
            SubjectLength[sseqid] = slen

            qstart = int(lines[qstartIndex])
            qend = int(lines[qendIndex])
            alignLen = abs(qend - qstart)
            if qstart < qend:
                region = (qstart, qend)
            else:
                region = (qend, qstart)

            SeqBlastRegions[qseqid][sseqid].add(region)
        
    blast_h.close()


    for seqid in SeqBlastRegions:
        SSregions = SeqBlastRegions[seqid]

        ### different sseqid and length
        SSLength = {}
        for sseqid in SSregions:
            regions = SSregions[sseqid]

            sortRegions = sorted(list(regions))
            StartEnd = longest_region(sortRegions)
            SliceList = combine_regions(StartEnd)
            OnceList = once_key(SliceList)
            Length = oneceList_length(OnceList)

            SSLength[sseqid] = Length

        ### longest length 
        sortedLength = sorted(SSLength.items(), key=operator.itemgetter(1), reverse=True)
        selectId, maxLength = sortedLength[0]


        ### query length
        qlen = QueryLength[seqid]
        print(qlen)
        ratio = maxLength / qlen
        
        slen = SubjectLength[sseqid]
        if ratio > lengthRatio:
            ratio = "%.3f" % ratio
            out_h.write("%s\t%s\t%d\t%d\t%s\n" % (seqid, selectId,  qlen, slen, ratio))

    out_h.close()


def oneceList_length(OnceList):
    UniqRegions = []
    Lengths = 0
    ListLen = len(OnceList)
    if ListLen == 2:
        UniqRegions = OnceList
        Lengths += UniqRegions[1] - UniqRegions[0]
    else:
        for i in range(0, int(ListLen/2)):
            start = OnceList[i*2]
            end = OnceList[i*2+1]
            UniqRegions.append([start, end])
            length = end - start
            Lengths += length
    return Lengths



def main():
    parser = argparse.ArgumentParser(description="Filt the blast record based on aligned length.")
    parser.add_argument("-i", "--input", help="The input file with blast aligned format.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-c", "--column", default="qseqid,sseqid,bitscore,evalue,length,pident,mismatch,gapopen,qstart,qend,sstart,send,qlen,slen,stitle", help="The blasted column names.")
    parser.add_argument("-l", "--lengthRatio", default=0.75, help="The threshold for aligned length.")
    args = parser.parse_args()
    filt_blast_result(args.input, args.column, args.out, args.lengthRatio)



if __name__ == "__main__":
    main()

#!/usr/bin/python
# -*- coding:utf-8
from __future__ import division
import collections
import bisect
import math
import argparse
import sys


#usage: python ~/github/NanoHub/src/NanoHub/WindowSVNumber.py --genome /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38_length.txt --input M671-1_position.xls --window 1000000 --sliding 0



def Chr_length(genome_len_file):
    """
    genome_len_file:
    1       248956422       248956422       1_248956422
    10      133797422       133797422       10_133797422
    11      135086622       135086622       11_135086622
    """
    ChrLength = {}
    genome_h = open(genome_len_file, "r")
    for line in genome_h:
        lines = line.strip().split("\t")
        Chr, length = lines[:2]
        length = int(length)
        ChrLength[Chr] = length
    genome_h.close()
    return ChrLength


def window_SV_number(genome_len_file, window, sliding):
    ChrLength = Chr_length(genome_len_file)

    ChrRegions = collections.defaultdict(list)
    for Chrom in ChrLength:
        # Chr_len = max(map(int, list(ChrPosRecord[Chrom].keys())))

        # try:
        #     Chr_len = ChrLength[Chrom]
        # except KeyError:
        #     print("Please check whether the chromosome name %s is in the file %s." % (Chrom, genome_len_file))
        #     sys.exit(1)

        Chr_len = ChrLength[Chrom]
        window = int(window)
        sliding = int(sliding)

        #window should be mutiple of sliding
        if sliding == 0:
            steps = 0
        else:
            steps = int(window / sliding)
        counts = math.ceil(Chr_len / window)

        for count in range(counts):
            if steps != 0:
                for s in range(steps):
                    Start = window * count + sliding * s + 1
                    End = window * (count +1) + sliding * s
                    region = (Start, End)
            else:
                Start = window * count + 1
                End = window * (count +1) 
                region = (Start, End)
            ChrRegions[Chrom].append(region)

    ### get the chromosome and start positions
    ChrStarts = {}

    for Chr in ChrRegions:
        regions = ChrRegions[Chr]
        Starts = [r[0] for r in regions]
        ChrStarts[Chr] = Starts

    return ChrRegions, ChrStarts


def SV_Type_bin_number(genome_len_file, SV_file, window, sliding, outprefix):
    """
    input example:
    1   368827  369147  DEL
    1   598153  598365  INS
    1   610176  610253  INS
    1   610597  610727  INS

    output file:
    1       1       1000000 6
    1       1000001 2000000 11
    1       2000001 3000000 25
    1       3000001 4000000 13
    1       4000001 5000000 2
    1       5000001 6000000 3
    1       6000001 7000000 3
    """
    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    TypeList = ["DEL", "INS", "INV", "DUP"]



    ChrRegions, ChrStarts = window_SV_number(genome_len_file, window, sliding)

    TypeChrPosRecord = collections.defaultdict(lambda: collections.defaultdict(dict))

    


    in_h = open(SV_file, 'r')
    Types = set()
    for line in in_h:
        line = line.strip()
        lines = line.split('\t')
        Chr, Start, End, Type = lines[:4]
        ### get the given chromosome record
        if Chr in CHRS and Type in TypeList:
            Start = int(Start)
            TypeChrPosRecord[Type][Chr][Start] = line
            Types.add(Type)
    in_h.close()


    TypeList = sorted(list(Types))
    for t in TypeList:
        ChrPosRecord = TypeChrPosRecord[t]

        ChrRegionNumber = collections.defaultdict(lambda: collections.Counter())

        for Chr in ChrPosRecord:
            Starts = ChrStarts[Chr]
            Regions = ChrRegions[Chr]

            ### get the index for start position of record in the STARTS list
            for pos in ChrPosRecord[Chr]:
                posIndex = bisect.bisect(Starts, pos)
                if posIndex == 0:
                    startIndex = 0
                elif posIndex == len(Starts):
                    startIndex = len(Starts) - 1
                else:
                    startIndex = posIndex - 1

                region = ChrRegions[Chr][startIndex]
                ChrRegionNumber[Chr][region] += 1

        ### output the file
        outfile = outprefix + "_" + t + ".xls"
        output_chrom_region(ChrRegionNumber, outfile)



def output_chrom_region(ChrRegionNumber, outprefix):
    out_h = open(outprefix, "w")
    sortedChrs = sorted(ChrRegionNumber)
    for c in sortedChrs:
        RegionNum = ChrRegionNumber[c]
        srotedRegions = sorted(list(RegionNum.keys()))
        for r in srotedRegions:
            num = ChrRegionNumber[c][r]
            if num != 0:
                r = [str(rr) for rr in r]
                out_h.write("%s\t%s\t%d\n" % (c, "\t".join(r), num))
    out_h.close()


        





def main():
    parser = argparse.ArgumentParser(description='Get the SV type and number for compared samples with sliding window.')
    parser.add_argument("-g", "--genome", help="The file contain genome chromosome length.")
    parser.add_argument('-i', '--input', help='The input file containing alt snp ratio.')
    parser.add_argument('-w', '--window', help='The window length.')
    parser.add_argument('-s', '--sliding', help='The sliding length')
    parser.add_argument('-o', '--out', help='The output prefix name.')
    args = parser.parse_args()
    SV_Type_bin_number(args.genome, args.input, args.window, args.sliding, args.out)


if __name__ == '__main__':
    main()


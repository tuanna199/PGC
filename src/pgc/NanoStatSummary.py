#!/usr/bin/env python
from __future__ import division
import argparse

#usage: python ~/github/TrioWGS/src/TrioWGS/NanoStatSummary.py --stats NanoStats.txt --out NanoStats_summary.xls --sample M625 --clean ../../../clean/NanoPlot/M625-0/NanoStats.txt

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.04.02"



def read_nanoStat_result(nanoStat):
    """
    nanoStat:

    General summary:        
    Mean read length:                17,210.8
    Mean read quality:                    8.2
    Median read length:              15,161.0
    Median read quality:                  8.5
    Number of reads:              2,468,218.0
    Read length N50:                 22,757.0
    Total bases:             42,479,947,530.0
    """
    in_h = open(nanoStat, "r")
    for line in in_h:
        lines = line.strip().split(":")
        if lines[0].startswith("Mean read length"):
            MeanReadLen = lines[-1].strip()
        elif lines[0].startswith("Mean read quality"):
            ReadQual = lines[-1].strip()
        elif lines[0].startswith("Median read length"):
            MedianReadLen = lines[-1].strip()
        elif lines[0].startswith("Number of reads"):
            ReadNum = lines[-1].strip()
        elif lines[0].startswith("Read length N50"):
            ReadN50 = lines[-1].strip()
        elif lines[0].startswith("Total bases"):
            TotalBase = lines[-1].strip()
    in_h.close()

    Stats = [ReadNum, ReadQual, MeanReadLen, MedianReadLen, ReadN50, TotalBase]

    ### delete comma instring using replace
    Stats = [float(s.replace(",", "")) for s in Stats]
    return Stats

def clean_to_raw_ratio(raw_stat, out_file, sample="sample", clean_stat=None):
    rawStats = read_nanoStat_result(raw_stat)
    ReadNum, ReadQual, MeanReadLen, MedianReadLen, ReadN50, TotalBase = rawStats

    out_h = open(out_file, "w")
    if not clean_stat:
        out_h.write("Sample\tRead_number\tRead_quality\tMean_read_length\tMedian_read_length\tRead_length_N50\tTotal_base\n")
        ReadNum = int(ReadNum)
        ReadQual = "%.1f" % ReadQual
        MeanReadLen = "%.1f" % MeanReadLen
        MedianReadLen = "%.1f" % MedianReadLen
        ReadN50 = int(ReadN50)
        TotalBase = int(TotalBase)
        out_h.write("%s\t%d\t%s\t%s\t%s\t%d\t%d\n" % (sample, ReadNum, ReadQual, MeanReadLen, MedianReadLen, ReadN50, TotalBase))
        out_h.close()
    else:
        cleanStats = read_nanoStat_result(clean_stat)
        cleanReadNum, cleanReadQual, cleanMeanReadLen, cleanMedianReadLen, cleanReadN50, cleanTotalBase = cleanStats
        readRatio = cleanReadNum / ReadNum * 100
        baseRatio = cleanTotalBase / TotalBase * 100
        readRatio = "%.1f" % readRatio
        baseRatio = "%.1f" % baseRatio

        out_h.write("Sample\tRaw_read_number\tRaw_read_quality\tRaw_mean_read_length\tRaw_median_read_length\tRaw_read_length_N50\tRaw_total_base\tClean_read_number\tClean_read_quality\tClean_mean_read_length\tClean_median_read_length\tClean_read_length_N50\tClean_total_base\tClean_read_ratio(%)\tClean_base_ratio(%)\n")
        out_h.write("%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\n" % (sample, ReadNum, ReadQual, MeanReadLen, MedianReadLen, ReadN50, TotalBase, cleanReadNum, cleanReadQual, cleanMeanReadLen, cleanMedianReadLen, cleanReadN50, cleanTotalBase, readRatio, baseRatio))
        out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Get the summary of NanoPloy statistics.")
    parser.add_argument("-s", "--stats", help="The input read statistics file derived from NanoPlot.")
    parser.add_argument("-c", "--clean", help="The input read statistics file derived from NanoPlot after clean.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-m", "--sample", default="sample", help="The sample name.")
    args = parser.parse_args()

    if args.clean:
        clean_to_raw_ratio(args.stats, args.out, sample=args.sample, clean_stat=args.clean)
    else:
        clean_to_raw_ratio(args.stats, args.out, sample=args.sample, clean_stat=None)

if __name__ == "__main__":
    main()











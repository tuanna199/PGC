#!/usr/env/bin python
import collections
import argparse

#usage: python ~/github/TrioWGS/src/TrioWGS/SVDistribution.py --input /home/wuzhikun/Project/NanoTrio/SVStats/Sniffles/minimap2/M628-2_position.xls --summary temp

__author__ = "Zhikun Wu"
__date__ = "2019.04.15"
__email__ = "598466208@qq.com"

def get_SV_infor(sv_bed):
    """
    sv_bed:
    1   76683   76789   DEL
    1   83940   84004   DEL
    1   136275  136340  INS
    """
    ChrTypeNum = collections.defaultdict(lambda: collections.Counter())

    in_h = open(sv_bed, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        Chr, Satrt, End, Type = lines[:4]
        ChrTypeNum[Chr][Type] += 1

    in_h.close()
    return ChrTypeNum

def SV_distribution(sv_bed, summary_file, chroms, types):
    ChrTypeNum = get_SV_infor(sv_bed)

    su_h = open(summary_file, "w")
    Chroms = chroms.split(",")
    Chroms = [c.strip() for c in Chroms]
    Types = types.split(",")
    Types = [t.strip() for t in Types]
    su_h.write("Chr\t%s\tTotal\n" % "\t".join(Types))
    for Chr in Chroms:
        if Chr in ChrTypeNum:
            typeNum = []
            for t in Types:
                if t in ChrTypeNum[Chr]:
                    number = ChrTypeNum[Chr][t]
                else:
                    number = 0
                typeNum.append(number)
            totalNum = sum(typeNum)
            typeNum = [str(t) for t in typeNum]
            su_h.write("%s\t%s\t%d\n" % (Chr, "\t".join(typeNum), totalNum))


def main():
    parser = argparse.ArgumentParser(description="Get the summary of structural variants distribution in chromosome.")
    parser.add_argument("-i", "--input", help="The input file with SV position and type.")
    parser.add_argument("-s", "--summary", help="The summary of SV distribution.")
    parser.add_argument("-c", "--chromosome", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y", help="The target chromosomes which are separated with comma.")
    parser.add_argument("-t", "--type", default="DEL,INS,DUP,INV,TRA,INVDUP", help="The target SV types.")
    args = parser.parse_args()
    SV_distribution(args.input, args.summary, chroms=args.chromosome, types=args.type)

if __name__ == "__main__":
    main()


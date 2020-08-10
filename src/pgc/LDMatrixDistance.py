#!/usr/bin/python
from WindowRegionNumber import window_variant_number
import numpy
import argparse
import collections

#usage: python ~/github/NanoHub/src/NanoHub/LDMatrixDistance.py --ld /home/wuzhikun/Project/NanoTrio/population/plink/Chrs/1_matrix.ld --map /home/wuzhikun/Project/NanoTrio/population/plink/Chrs/Sample_SV_geno_1.map --out temp.txt --window 100000 --sliding 0 --winOut temp_win.txt  --lenThreshold 1000000 


__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.10.18"

def read_ld_matrix(ld_file, map_file, out_file, window_out, window, sliding, lenThreshold):
    """
    map_file:
    19  19_140239_19_140392_DEL 0   140239
    19  19_183027_19_183090_DEL 0   183027

    out_file:
    Start   End Number  Value
    1   10000   9743    0.018
    10001   20000   9284    0.016
    20001   30000   8378    0.016
    """
    ### parser the LD matrix
    ### np.loadtxt("input.txt", dtype='i', delimiter=',')
    LDMatrix = numpy.loadtxt(ld_file)

    ### parser the marker position
    map_h = open(map_file, "r")
    IndexPos = {}
    for index, line in enumerate(map_h.readlines()):
        lines = line.strip().split("\t")
        pos = lines[-1]
        pos = int(pos)
        IndexPos[index] = pos
    map_h.close()

    PosLDValue = collections.defaultdict(list)
    for i in range(len(LDMatrix)):
        for j in range(len(LDMatrix[0])):
            if j > i:
                ldValue = LDMatrix[i][j]
                pos1 = IndexPos[i]
                pos2 = IndexPos[j]
                distance = abs(pos1 - pos2)
                PosLDValue[distance] = ldValue

    maxDist = max(list(PosLDValue.keys()))

    RegionNumber, RegionValue = window_variant_number(PosLDValue, window, sliding)

    positions = sorted(RegionValue.keys())

    Records = []
    for p in positions:
        value = RegionValue[p]
        number = RegionNumber[p]
        ratio = value / number
        ratio = "%.3f" % ratio
        p = [str(r) for r in p]
        Records.append("%s\t%d\t%s" % ("\t".join(p), number, ratio))

    return PosLDValue, Records





def main():
    parser = argparse.ArgumentParser(description="Calculate LD and distance.")
    parser.add_argument("-l", "--ld", help="The input ld matrix file.")
    parser.add_argument("-m", "--map", help="The map file containing marker information.")
    parser.add_argument("-w", "--window", default=100000, help="The window length.")
    parser.add_argument("-s", "--sliding", default=0, help="The sliding length for window.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-wo", "--winOut", help="The output file with window region.")
    parser.add_argument("-t", "--lenThreshold", default=1000000, help="The threshold for cut length for distance of LD.")
    args = parser.parse_args()
    PosLDValue, Records = read_ld_matrix(args.ld, args.map, args.out, args.winOut, args.window, args.sliding, args.lenThreshold)

    win_h = open(args.winOut, "w")
    win_h.write("Start\tEnd\tNumber\tValue\n")
    for record in Records:
        win_h.write("%s\n" % record)
    win_h.close()

    lenThreshold = int(args.lenThreshold)
    out_h = open(args.out, "w")
    out_h.write("Pos\tr2\n")
    for pos in sorted(list(PosLDValue.keys())):
        if pos < lenThreshold:
            ld = PosLDValue[pos]
            out_h.write("%d\t%f\n" % (pos, ld))
    out_h.close()

if __name__ == "__main__":
    main()


#!/usr/bin/python
from __future__ import division
import collections
import argparse
import math
import os
import bisect



#usage: python ~/github/NanoHub/src/NanoHub/SVTerminalDistance.py --vcf /home/wuzhikun/Project/NanoTrio/population/Sample_common_SV_convert.vcf --out temp_dist.xls --window  100000 --sliding 0 --cytoband /home/wuzhikun/database/Annovar/hg38/hg38_cytoBand.txt

def overlap_regions(pos, regions):
    targetRegions = []
    for region in regions:
        posIndex = bisect.bisect(region, pos)
        if posIndex == 1:
            targetRegions.append(region)
    return targetRegions


def window_variant_number(ChrPosRecord, out_file, window, sliding, ChrPQLen):
    """
    input example:
    X       10072   .       G       T       33.90   .       AC=1;AF=0.500;AN=2;BaseQRankSum=-8.420e-01;ClippingRankSu
    m=0.00;DP=9;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=35.45;MQRankSum=-2.189e+00;QD=4.84;ReadPosRankSum=2.
    19;SOR=0.446  GT:AD:DP:GQ:PL  ./.:0,0:0:.:0,0,0       ./.:2,0:2:.:0,0,0       0/1:4,3:7:60:60,0,95

    out_file:
    Chr     Start   End     Number
    1       0.00000 35
    1       0.00810 65
    1       0.01621 150
    """
    ChrRegions = collections.defaultdict(list)
    ChrRegionNumber = collections.defaultdict(lambda: collections.Counter())
    ChrWinds = {}


    ### get the window regions
    for Chrom in ChrPosRecord:
        Chr_len = max(map(int, list(ChrPosRecord[Chrom].keys())))
        window = int(window)
        sliding = int(sliding)
        #window should be mutiple of sliding
        if sliding == 0:
            steps = 0
        else:
            steps = int(window / sliding)
        counts = math.ceil(Chr_len / window)
        ChrWinds[Chrom] = counts

        for count in range(counts):
            if steps != 0:
                for s in range(steps):
                    Start = window * count + sliding * s + 1
                    End = window * (count +1) + sliding * s
                    region = (Start, End)
                    ChrRegions[Chrom].append(region)
            else:
                Start = window * count + 1
                End = window * (count +1) 
                region = (Start, End)
                ChrRegions[Chrom].append(region)

    ### if sliding is zero, we just use window
    if sliding == 0:
        slidingFold = 0
    else:
        if window % sliding == 0:
            slidingFold  = int(window / sliding)
        else:
            print("Window length should be integral multiple number of sliding length.")
            sys.exit(1)   


    ### assign the variants number to the regions
    for Chr in ChrPosRecord:
        allRegions = ChrRegions[Chr]
        regionStarts = [r[0] for r in allRegions]

        for pos in ChrPosRecord[Chr]:
            ### get the start region for pos
            satrtIndex = bisect.bisect(regionStarts, pos)
            startPre = satrtIndex - slidingFold -1
            if startPre < 0:
                startPre = 0

            searchRegions = [] 
            for s in range(startPre, satrtIndex):
                seRegion = allRegions[s]
                searchRegions.append(seRegion)

            targetRegions = overlap_regions(pos, searchRegions)
            if targetRegions != []:
                for t in targetRegions:
                    # ChrRegionNumber[Chr][t] += 1

                    ### normalization
                    t1, t2 = t
                    centroLen, rightLen = ChrPQLen[Chr]
                    if t2 < centroLen:
                        normdist1 = t1 /centroLen
                        normdist2 = t2 /centroLen
                    else:
                        normdist1 = (t1 - centroLen) / rightLen + 1
                        normdist2 = (t2 - centroLen) / rightLen + 1

                    ChrRegionNumber[Chr][(normdist1, normdist2)] += 1
                    normdist1 = "%.4f" % normdist1
                    normdist2 = "%.4f" % normdist2
                    ChrRegionNumber[Chr][(normdist1, normdist2)] += 1


    ### output the result
    out_h = open(out_file, 'w')
    # out_h.write('Chr\tStart\tEnd\tNumber\n')
    out_h.write('Chr\tPos\tNumber\n')
    sortedChrs = sorted(list(ChrRegionNumber.keys()))
    for c in sortedChrs:
        sortedRegions = ChrRegionNumber[c]
        for r in sortedRegions:
            number = ChrRegionNumber[c][r]
            # out_h.write("%s\t%s\t%d\n" % (c, "\t".join(map(str, list(r))), number))
            out_h.write("%s\t%s\t%d\n" % (c, r, number))
    out_h.close()




def variant_number(ChrPosRecord, out_file, ChrPQLen):
    """
    input example:
    X       10072   .       G       T       33.90   .       AC=1;AF=0.500;AN=2;BaseQRankSum=-8.420e-01;ClippingRankSu
    m=0.00;DP=9;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=35.45;MQRankSum=-2.189e+00;QD=4.84;ReadPosRankSum=2.
    19;SOR=0.446  GT:AD:DP:GQ:PL  ./.:0,0:0:.:0,0,0       ./.:2,0:2:.:0,0,0       0/1:4,3:7:60:60,0,95

    out_file:
    Chr     Start   End     Number
    1       0.00000 35
    1       0.00810 65
    1       0.01621 150
    """
    ChrRegions = collections.defaultdict(list)
    ChrRegionNumber = collections.defaultdict(lambda: collections.Counter())


    ### assign the variants number to the regions
    for Chr in ChrPosRecord:
        centroLen, rightLen = ChrPQLen[Chr]

        for pos in ChrPosRecord[Chr]:
            if pos < centroLen:
                normdist = pos /centroLen
            else:
                normdist = (pos - centroLen) / rightLen + 1
            ChrRegionNumber[Chr][normdist] += 1


    ### output the result
    out_h = open(out_file, 'w')
    # out_h.write('Chr\tStart\tEnd\tNumber\n')
    out_h.write('Chr\tPos\tNumber\n')
    sortedChrs = sorted(list(ChrRegionNumber.keys()))
    for c in sortedChrs:
        sortedRegions = ChrRegionNumber[c]
        for r in sortedRegions:
            number = ChrRegionNumber[c][r]
            # out_h.write("%s\t%s\t%d\n" % (c, "\t".join(map(str, list(r))), number))
            out_h.write("%s\t%s\t%d\n" % (c, r, number))
    out_h.close()





def acen_position(regions):
    region1, region2 = regions
    position = 0
    if region1[0] == region2[-1]:
        position = region1[0]
    elif region1[-1] == region2[0]:
        position = region1[-1]
    return position



def centromere_position(cytoBand):
    """
    cytoband:
    1       121700000       123400000       p11.1   acen
    1       123400000       125100000       q11     acen
    10      38000000        39800000        p11.1   acen
    10      39800000        41600000        q11.1   acen

    ChrPos = {"1":123400000, ...}
    """
    ChrPos = {}
    ChrPQLen = {}

    ChrRegions = collections.defaultdict(list)
    ChrEnds = collections.defaultdict(list)
 
    cyto_h = open(cytoBand, "r")
    for line in cyto_h:
        lines = line.strip().split("\t")
        Chr, Satrt, End = lines[:3]
        Satrt = int(Satrt)
        End = int(End)
        Type = lines[-1]

        ChrEnds[Chr].append(End)

        if Type == "acen":
            ChrRegions[Chr].append([Satrt, End])
    cyto_h.close()


    for Chr in ChrRegions:
        ### get the last position
        Ends = ChrEnds[Chr]
        lastPos = max(Ends)

        ### get the center position
        regions = ChrRegions[Chr]
        regionLength = len(regions)
        if regionLength == 2:
            position = acen_position(regions)
            if position != 0:
                ### build the dict for center and last position for each chromosome
                ChrPos[Chr] = (position, lastPos)

                ### get the length of chromosome long and short arm
                rightLen = lastPos - position
                ChrPQLen[Chr] = [position, rightLen]
            else:
                print("Please check whether p and q arm have not identical acen position for chromosome %s." % Chr)
        else:
            print("Please check why the acen regions are not two for chromosome %s." % Chr)

    return ChrPos, ChrPQLen



def sv_distance_stats(cytoBand, vcf_file,  outfile, window, sliding):

    CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

    window = int(window)
    sliding = int(sliding)



    ChrPosRecord = collections.defaultdict(lambda: collections.Counter())

    ### centro potision for each chromosome
    ChrPos, ChrPQLen = centromere_position(cytoBand)
    

    vcf_h = open(vcf_file, "r")
    for line in vcf_h:
        line = line.strip()
        if not line.startswith("#"):
            lines = line.strip().split("\t")
            Chr, Pos, ID = lines[:3]
            if Chr in CHRS:
                Pos = int(Pos)
                SVType = lines[4]
                SVType = SVType.lstrip("<").rstrip(">")

                center, lastPos = ChrPos[Chr]

                # ### normalization for distance to [0, 1]
                # if Pos <= center:
                #     dist = Pos 
                # else:
                #     dist = abs(lastPos - Pos) 

                ChrPosRecord[Chr][Pos] += 1
    vcf_h.close()
    if window == 0:
        variant_number(ChrPosRecord, outfile, ChrPQLen)
    else:
        window_variant_number(ChrPosRecord, outfile, window, sliding, ChrPQLen)


def main():
    parser = argparse.ArgumentParser(description="Get the number of variants across region in the genome.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-c", "--cytoband", help="The input cytoband file.")
    parser.add_argument("-o", "--out", help="The output file")
    parser.add_argument("-w", "--window", default=100000, help="The window length of the genome.")
    parser.add_argument("-s", "--sliding", default=0, help="The sliding length in each window.")
    args = parser.parse_args()
    sv_distance_stats(args.cytoband, args.vcf, args.out, args.window, args.sliding)


if __name__ == "__main__":
    main()

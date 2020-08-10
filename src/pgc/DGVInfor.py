#!/usr/bin/env python
import collections
import argparse

def get_DGV_SV_infor(DGV_file, out_file):
    ChrRegionSV = collections.defaultdict(lambda: collections.defaultdict(list))

    in_h = open(DGV_file, "r")
    header = in_h.readline().strip()
    for line in in_h:
        lines = line.strip().split("\t")
        SVId, Chr, Start, End, SVType = lines[:5]
        supportVar = lines[11]
        ChrRegionSV[Chr][(Start, End)].append([SVId, SVType, supportVar])
    in_h.close()


    out_h = open(out_file, "w")
    out_h.write("Chr\tStart\tEnd\tDGV_infor(ID,SVTYPE,supportVar)\n")
    for c in ChrRegionSV:
        regions = list(ChrRegionSV[c].keys())
        sortRegions = sorted(regions)
        print(sortRegions)
        for region in sortRegions:
            SVs = ChrRegionSV[c][region]
            SVs = [",".join(s) for s in SVs]
            out_h.write("%s\t%s\t%s\n" % (c, "\t".join(map(str, list(region))), ";".join(SVs)))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Extract the main information of SV for dbVar database.")
    parser.add_argument("-i", "--input", help="The dbVar database with vcf format.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    get_DGV_SV_infor(args.input, args.out)

if __name__ == "__main__":
    main()
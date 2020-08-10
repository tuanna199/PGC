#!/usr/bin/env python
import argparse
import collections
from BaseFunc import Infor_target_values

#usage: python ~/github/TrioWGS/src/TrioWGS/dbVarInfor.py --vcf GRCh38.variant_region.vcf  --out GRCh38.variant_region.xls


def extract_dbVar_SV_infor(dbVar_file, out_file):
    """
    dbVar:
    1       10001   nsv482937       T       <CNV>   .       .       DBVARID;SVTYPE=CNV;IMPRECISE;END=2368561;CIPOS=.,0;CIEND=.,0
    1       10001   nsv482957       T       <CNV>   .       .       DBVARID;SVTYPE=CNV;IMPRECISE;END=15873505;CIPOS=.,0;CIEND=.,0
    1       10001   nsv483003       T       <CNV>   .       .       DBVARID;SVTYPE=CNV;IMPRECISE;END=15873505;CIPOS=.,0;CIEND=.,0

    out_file:
    1   10001   19818   nsv945697,CNV
    1   10001   22118   nsv945698,CNV
    1   10001   127330  nsv7879,CNV
    1   10001   297968  nsv870619,CNV
    """
    ChrRegionSV = collections.defaultdict(lambda: collections.defaultdict(list))
    in_h = open(dbVar_file, "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start, SVId = lines[:3]
            Start = int(Start)
            Infor = lines[7]
            SVType, End = Infor_target_values(Infor, "SVTYPE,END")
            End = int(End)

            ChrRegionSV[Chr][(Start, End)].append([SVId, SVType])
    in_h.close()

    out_h = open(out_file, "w")
    out_h.write("Chr\tStart\tEnd\tdbVar_infor(id,SVTYPE)\n")
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
    parser.add_argument("-v", "--vcf", help="The dbVar database with vcf format.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    extract_dbVar_SV_infor(args.vcf, args.out)

if __name__ == "__main__":
    main()
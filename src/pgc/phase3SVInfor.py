#!/usr/bin/env python
import sys
import argparse
import collections
from BaseFunc import Infor_target_values

#usage: python ~/github/TrioWGS/src/TrioWGS/dbVarInfor.py --vcf GRCh38.variant_region.vcf  --out GRCh38.variant_region.xls


def extract_dbVar_SV_infor(phaseSV_file, out_file):
    """
    1       710330  ALU_umary_ALU_2 A       <INS:ME:ALU>    .       .       AC=35;AF=0.00698882;AFR_AF=0;AMR_AF=0.0072;AN=5008;CS=ALU_umary;EAS_AF=0.0069;EUR_AF=0.0189;MEINFO=AluYa4_5,1,223,-;NS=2504;SAS_AF=0.0041;SITEPOST=0.9998;SVLEN=222;SVTYPE=ALU;TSD=null GT      0|0

    out_file:
    1   10001   19818   nsv945697,CNV
    1   10001   22118   nsv945698,CNV
    1   10001   127330  nsv7879,CNV
    1   10001   297968  nsv870619,CNV
    """
    ChrRegionSV = collections.defaultdict(lambda: collections.defaultdict(list))
    in_h = open(phaseSV_file, "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start, SVId = lines[:3]
            Start = int(Start)
            Infor = lines[7]
            if "END" in Infor:
                SVTYPE,END,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF = Infor_target_values(Infor, "SVTYPE,END,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF")
                End = int(End)
            elif "SVLEN" in Infor:
                SVTYPE,SVLEN,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF = Infor_target_values(Infor, "SVTYPE,SVLEN,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF")
                End = Start + abs(int(SVLEN))
            else:
                print("Please check whether the tag 'SVLEN' or 'END' is in the SV information %s." % Infor)
                sys.exit(1)


            ChrRegionSV[Chr][(Start, End)].append([SVId,SVTYPE, AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF])
    in_h.close()

    out_h = open(out_file, "w")
    out_h.write("Chr\tStart\tEnd\tphase3SV_infor(SVId,SVTYPE,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF)\n")
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
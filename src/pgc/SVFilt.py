#!/usr/bin/env python
from BaseFunc import Infor_target_values
from BaseFunc import parse_bed_regions
from BaseFunc import overlap_chrom_regions
import sys
import argparse
import collections
import bisect

#usage: python ~/github/TrioWGS/src/TrioWGS/SVFilt.py --vcf M625-0.vcf --region /home/wuzhikun/database/Annovar/hg38/hg38_cytoBand_acen_region.txt --out M625-0_filt.vcf

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.04.16"



def vcf_filt_SV(vcf_file, region_file, out_file, filtLength):
    filtLength = int(filtLength)

    ChrRegions = parse_bed_regions(region_file)

    out_h = open(out_file, "w")
    in_h = open(vcf_file, "r")
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            out_h.write("%s\n" % line)
        else:
            lines = line.split("\t")
            Chr = lines[0]
            start = lines[1]
            Infor = lines[7]
            end = Infor_target_values(Infor, "END")
            Chr2 = Infor_target_values(Infor, "CHR2")
            start = int(start)
            end = int(end)
            ### detect whether this record is overlap with the regions
            is_exist = overlap_chrom_regions(Chr, start, end, ChrRegions)
            if is_exist == True:
                continue
            else:
                SVlength = abs(end-start)
                #################################
                ########### filt TRA ############
                if Chr == Chr2 and SVlength >= filtLength:
                    print("Please check whether the SV record should be filted:\n%s" % line)
                else:
                    out_h.write("%s\n" % line)
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Filt out the SV that existed in the given regions.")
    parser.add_argument("-v", "--vcf", help="The input vcf file.")
    parser.add_argument("-r", "--region", help="The given region that should bed removed, with the format like bed.")
    parser.add_argument("-f", "--filtLength", default=50000000, help="Filt the SV with length larger than this threshold.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    vcf_filt_SV(args.vcf, args.region, args.out, args.filtLength)

if __name__ == "__main__":
    main()
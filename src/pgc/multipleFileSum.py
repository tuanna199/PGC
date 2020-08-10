#!/usr/bin/python
from __future__ import division
import argparse
import sys
import collections
import subprocess

#usage: python ~/github/NanoHub/src/NanoHub/multipleFileSum.py --vcf M668-2_filt_types.vcf,M669-1_filt_length.vcf,M669-2_filt_depth.vcf,M669-2_filt_length.vcf,M669-2_filt_types.vcf --out temp



def file_type_sv_number(vcfs, TypeFileNumber):
    files = vcfs.split(",")
    files = [f.strip() for f in files]

    for file in files:
        f = file.split("/")[-1]
        fileName = f.split(".")[0].split("_")[0]
        fileType = f.split("_")[-1].split(".")[0]

        if fileName == "" or fileType == "":
            print("Please check the file name and make sure to get name and type of file %s." % file)
            sys.exit(1)

        cmd = "grep -cv '^#' %s" % file
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        out, error = p.communicate()
        number = out.decode("utf-8").split("\n")[0]
        TypeFileNumber[fileType][fileName] = number


def multiple_files_summary(vcfs, out_file):
    TypeFileNumber = collections.defaultdict(dict)
    file_type_sv_number(vcfs, TypeFileNumber)

    types = sorted(list(TypeFileNumber.keys()))

    out_h = open(out_file, "w")
    out_h.write("Sample\t%s\n" % "\t".join(types))

    Samples = set()
    for t in types:
        sample = set(list(TypeFileNumber[t].keys()))
        Samples = Samples.union(sample)

    SAMPLES = sorted(list(Samples))

    for s in SAMPLES:
        Numbers = []
        for t in types:
            if s in TypeFileNumber[t]:
                num = TypeFileNumber[t][s]
            else:
                num = "0"
            Numbers.append(num)
        out_h.write("%s\t%s\n" % (s, "\t".join(Numbers)))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Summary SV number for multiple files of multiple samples.")
    parser.add_argument("-v", "--vcf", help="The input vcf files.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    multiple_files_summary(args.vcf, args.out)




if __name__ == "__main__":
    main()


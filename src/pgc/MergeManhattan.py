#!/usr/bin/python
import collections
import argparse
import sys
import os


def merge_multiple_traits_association(fileStr, out_file, target):
    """
    one file:
    SNP CHR BP  P
    1_136988-1_137121-148-INS   1   136988  0.4811
    1_180194-1_180382-172-DEL   1   180194  0.2447
    """
    targets = target.split(",")
    targets = [t.strip() for t in targets]

    TraitRecordPvalue = collections.defaultdict(dict)

    # NumTraits = {"2": "ASAL", "3": "AST", "4": "GCT", "5": "GLU", "11": "BA#", "12": "BA%", "19": "MCH", "21": "MCV", "24": "MPV", "31": "RDW-CV", "38": "TG", "39": "EBNA1-IgA", "40": "Zta-IgA", "42": "EA-IgA", "45": "UBACT"}
    files = fileStr.split(",")

    Records = set()
    for f in files:
        f = f.strip()
        fbase = f.split("/")[-1]
        # fnumber = f.rstrip(".assoc.txt").split("_")[-1]
        # if fnumber in targets:
        #     try:
        #         trait = NumTraits[fnumber]
        #     except KeyError:
        #         print("Please ckeck whether the numer of trait '%s' exist in the list." % fnumber)
        #         sys.exit(1)
        trait = fbase.split(".")[0]
        if trait in targets:
            in_h = open(f, "r")
            header = in_h.readline().strip()
            for line in in_h:
                line = line.strip()
                lines = line.split("\t")
                SV, Chr, Pos, Pvalue = lines
                record = "%s\t%s\t%s" % (SV, Chr, Pos)
                Records.add(record)
                TraitRecordPvalue[trait][record] = Pvalue
            in_h.close()

    out_h = open(out_file, "w")
    traits = sorted(list(TraitRecordPvalue.keys()))
    out_h.write("SNP\tCHR\tBP\t%s\n" % "\t".join(traits))

    RecordAll = sorted(list(Records))
    for r in RecordAll:
        Pvalues = []
        for t in traits:
            if r in TraitRecordPvalue[t]:
                pvalue = TraitRecordPvalue[t][r]
            else:
                pvalue = "NA"
            Pvalues.append(pvalue)
        out_h.write("%s\t%s\n" % (r, "\t".join(Pvalues)))
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Merge the pvalues of multiple GWAS traits.")
    parser.add_argument("-f", "--file", help="The input files which were separated with ','.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--target", help="The target number of traits.")
    args = parser.parse_args()
    merge_multiple_traits_association(args.file, args.out, args.target)

if __name__ == "__main__":
    main()






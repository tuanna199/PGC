#!/usr/bin/env python
import argparse

#usage: python ~/github/NanoHub/src/NanoHub/TRA2circos.py --input  M628-1_tra.xls --out temp.xls --method two


def TRA_to_circos_link(tra_file, out_file, chromosome, method):
    """
    tra_file
    1   10004   3   198172702
    1   10331   X   449317
    1   10469   11  175280
    """
    method = method.lower()

    sra_h = open(tra_file, "r")
    chroms = chromosome.split(",")
    chroms = [c.strip() for c in chroms]
    out_h = open(out_file, "w")
    for line in sra_h:
        lines = line.strip().split("\t")
        lineLength = len(lines)
        if method == "two":
            Chr1, Pos1, Chr2, Pos2 = lines[:4]
            if (Chr1 in chroms) and (Chr2 in chroms):
                record = "hs%s\t%s\t%s\ths%s\t%s\t%s" % (Chr1, Pos1, Pos1, Chr2, Pos2, Pos2)
                out_h.write("%s\n" % record)
        elif method == "one":
            Chr, Start, End, Value = lines[:4]
            if Chr in chroms:
                record = "hs%s\t%s\t%s\t%s" % (Chr, Start, End, Value)
                out_h.write("%s\n" % record)
        else:
            print("Please make sure that method is 'one' or 'two'.")
            sys.exit(1)
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Convert the link of translocation to circos format.")
    parser.add_argument("-i", "--input", help="The input file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-c", "--chromosome", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y", help="Tgaet chromosomes.")
    parser.add_argument("-m", "--method", help="The method that change the given columns,'one' or 'two' for chromosome.")
    args = parser.parse_args()
    TRA_to_circos_link(args.input, args.out, args.chromosome, args.method)

if __name__ == "__main__":
    main()

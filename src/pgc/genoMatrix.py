#!/usr/bin/python
import argparse
import collections
import sys
import os
import numpy

#usage: python ~/github/NanoHub/src/NanoHub/genoMatrix.py --genotype Sample_SV_genotype_filt_01.txt --out temp

def genotype_to_matirx(geno_file, out_file):
    """
    hr1    Pos1    Chr2    Pos2    SVlength        Type    CN001   CN002   CN003   CN004   CN005   CN006   CN007   CN008   CN009   CN010   CN011   CN012   CN013   CN014   
    1       90312   1       90388   74.3    INS     0/0     0/0     0/1     0/0     0/0     0/1     0/1     0/1     0/0     0/0     0/0     0/0     0/0     0/1     0/0     
    1       136985  1       137162  137.7   INS     1/1     0/1     0/1     0/1     0/1     0/1     0/1     0/1     0/1     0/1     0/1     0/1     0/0     1/1     0/1     
    1       181217  1       181455  236.6   DEL     0/0     0/0     0/0     0/1     0/0     0/1     0/0     0/0     0/1     0/0     0/1     0/1     0/0     0/1     0/0     
    1       257928  1       258047  97.4    INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/1     0/0     
    1       368919  1       369167  323.3   DEL     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     
    """
    geno_h = open(geno_file, "r")
    headers = geno_h.readline().strip().split("\t")
    samples = headers[6:]

    Matrix = []
    for line in geno_h:
        lines = line.strip().split("\t")
        genos = lines[6:]

        NewGenos = []
        for g in genos:
            gg = g.split("/")
            gg = [int(i) for i in gg]
            num = sum(gg) + 1
            NewGenos.append(num)
        Matrix.append(NewGenos)
    geno_h.close()

    m = numpy.array(Matrix)
    Trans = m.transpose()

    ### output the matrix
    out_h = open(out_file, "w")

    TranLen = len(Trans)
    sampleLen = len(samples)
    if TranLen == sampleLen:
        for s, v in zip(samples, Trans):
            v = [str(i) for i in v]
            out_h.write("%s\t%s\n" % (s, "\t".join(v)))
    else:
        print("Please make sure that the length of samples and matrix is identical.")
        sys.exit(1)
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Translocate the genotype to matrix.")
    parser.add_argument("-g", "--genotype", help="The input genotype file")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    genotype_to_matirx(args.genotype, args.out)



if __name__ == "__main__":
    main()

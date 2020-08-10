#!/usr/bin/python
import collections
import argparse
import os

#usage: python ~/github/NanoHub/src/NanoHub/eigenstrat2PCA.py --ind Sample_SV.ind --pca Sample_SV.pca.txt --out temp_out.txt

def parse_pca(pca_file):
    """
    10
    1.3550
    1.3100
    1.2770
    1.2580
    1.2480
    1.2460
    1.2370
    1.2250
    1.2190
    1.2120
     0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
     -0.0142 -0.0131 -0.0190  0.0136  0.0181  0.0237 -0.0057 -0.0100 -0.0376 -0.0062
      0.0573 -0.0176  0.0036  0.0108  0.0428 -0.0025  0.0113 -0.0210 -0.0371 -0.0815
     -0.0091 -0.0110 -0.0130  0.0156  0.0181  0.0110  0.0215 -0.0213 -0.0272 -0.0380
    """
    pca_h = open(pca_file, "r")
    count = 1
    IndexPCA = {}
    for line in pca_h:
        lines = line.strip().split()
        pcaLen = len(lines)
        if pcaLen > 1:
            IndexPCA[count] = lines
            count += 1
            pcaNum = len(lines)
    pca_h.close()
    return IndexPCA, pcaNum

def sample_information(ind_file):
    """
    M509-1  M   Control
    M509-2  F   Control
    M561-0  F   Case
    """
    ind_h = open(ind_file, "r")
    count = 1
    IndexSample = {}
    for line in ind_h:
        lines = line.strip().split("\t")
        sample = lines[0]
        IndexSample[count] = sample
        count += 1
    ind_h.close()
    return IndexSample

def  sample_pca(pca_file, ind_file, out_file):
    IndexSample = sample_information(ind_file)
    IndexPCA, pcaNum = parse_pca(pca_file)
    index1 = sorted(list(IndexSample.keys()))
    index2 = sorted(list(IndexPCA.keys()))

    out_h = open(out_file, "w")
    pcaNames = ["PCA" + str(i+1) for i in range(pcaNum)]
    out_h.write("Sample\t%s\n" % "\t".join(pcaNames))
    if index1 == index2:
        for i in index1:
            sample = IndexSample[i]
            pca = IndexPCA[i]
            out_h.write("%s\t%s\n" % (sample, "\t".join(pca)))
    else:
        print("Please check whether sample number in geno file %s and pca file %s is identical." % (pca_file, ind_file))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the sampel and pca information.")
    parser.add_argument("-p", "--pca", help="The pca information.")
    parser.add_argument("-g", "--ind", help="The indival file for sample.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    sample_pca(args.pca, args.ind, args.out)

if __name__ == "__main__":
    main()


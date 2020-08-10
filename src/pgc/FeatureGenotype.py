#!/usr/bin/python
import argparse
import collections
import os
import sys

#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/FeatureGenotype.py --genotype /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype.txt --annotation /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_genepred_overlap_enhancer_all.bed --outPrefix /home/wuzhikun/Project/Population/population/genotype/feature/Sample_SV_genotype > /home/wuzhikun/Project/Population/log/FeatureGenotype.log 2>&1


def feature_tag(anno_file):
    """
    anno_file:
    1   88889   88955   1_88889-1_88956-66-INS  1   88550   89550   RP11-34P13.8    ENST00000495576.1   DownStream
    1   90312   90410   1_90312-1_90313-98-INS  1   89294   91629   RP11-34P13.7    ENST00000466430.5   Exon
    1   90312   90410   1_90312-1_90313-98-INS  1   90286   91105   RP11-34P13.8    ENST00000495576.1   Exon
    """
    FeatureTags = collections.defaultdict(set)
    anno_h = open(anno_file, "r")
    for line in anno_h:
        lines = line.strip().split("\t")
        if lines != []:
            tag = lines[3]
            feature = lines[-1]

            if "utr" in feature.lower():
                feature = "UTR5_3"

            if "stream" in feature.lower():
                feature = "Up_DownStream"

            if feature == "Exon":
                feature = "NC_Exon"

            FeatureTags[feature].add(tag)
    anno_h.close()
    return FeatureTags

def divide_genotypes(anno_file, geno_file, outPrefix):
    FeatureTags = feature_tag(anno_file)

    geno_h = open(geno_file, "r")
    header = geno_h.readline().strip()

    TagRecrod = {}
    for line in geno_h:
        line = line.strip()
        lines = line.split("\t")
        infor = lines[:6]
        tag = "%s_%s-%s_%s-%s-%s" % tuple(infor)
        TagRecrod[tag] = line
    geno_h.close()

    outdir = "/".join(outPrefix.split("/")[:-1])
    outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    features = sorted(list(FeatureTags.keys()))
    for f in features:
        if "/" in f:
            ff = f.replace("/", "_")
            out_h = open("%s_%s.txt" % (outPrefix, ff), "w")
        else:
            out_h = open("%s_%s.txt" % (outPrefix, f), "w")

        out_h.write("%s\n" % header)

        tags = FeatureTags[f]
        for t in sorted(list(tags)):
            try:
                record = TagRecrod[t]
                out_h.write("%s\n" % record)
            except KeyError:
                print("Please ckeck whether the tag %s is in file %s." % (t, geno_file))
                sys.exit(1)
        out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Split the genotype file based on the feature of tag.")
    parser.add_argument("-g", "--genotype", help="The input genotype file.")
    parser.add_argument("-a", "--annotation", help="The annotation file containing tag and gene feature.")
    parser.add_argument("-o", "--outPrefix", help="The output frefix.")
    args = parser.parse_args()
    divide_genotypes(args.annotation, args.genotype, args.outPrefix)





if __name__ == "__main__":
    main()
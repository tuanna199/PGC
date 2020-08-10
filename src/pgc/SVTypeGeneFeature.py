#!/usr/bin/python
import argparse
import collections
import sys
import os


def SV_type_gene_feature(anno_file, out_file, circos_file):
    GeneSets = set()
    FeatureSets = set()
    SVFeatureGenes = collections.defaultdict(lambda: collections.defaultdict(set))
    FeatureGenes = collections.defaultdict(set)
    TypeTag = collections.defaultdict(set)

    GeneFeatureType = collections.defaultdict(lambda: collections.defaultdict(set))

    anno_h = open(anno_file, "r")
    for line in anno_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        SVType = tag.split("-")[-1]
        gene = lines[7]
        feature = lines[-1]
        if feature == "UTR3" or feature == "UTR5":
            feature = "UTR"
        TypeTag[SVType].add(tag)
        SVFeatureGenes[SVType][feature].add(gene)
        GeneSets.add(gene)
        FeatureSets.add(feature)
        FeatureGenes[feature].add(gene)
        GeneFeatureType[gene][feature].add(SVType)
    anno_h.close()

    out_h = open(out_file, "w")
    types = sorted(list(SVFeatureGenes.keys()))
    features = sorted(list(FeatureSets))
    for t in types:
        for f in features:
            if f in SVFeatureGenes[t]:
                genes = SVFeatureGenes[t][f]
                geneNum = len(list(genes))
            else:
                geneNum = 0
            out_h.write("%s\t%s\t%d\n" % (t, f, geneNum))
    out_h.close()


    ####
    circos_h = open(circos_file, "w")
    FeatureTypesCount = collections.defaultdict(lambda: collections.Counter())
    for f in features:
        for g in GeneFeatureType:
            if f in GeneFeatureType[g]:
                t = tuple(sorted(list(GeneFeatureType[g][f])))
                FeatureTypesCount[f][t] += 1


    for f in features:
        typeCount = collections.Counter()
        intial = 0
        tts = sorted(list(FeatureTypesCount[f].keys()))
        for t in tts:
            c = FeatureTypesCount[f][t]
            ts = list(t)
            end = intial + c
            for t in ts:
                tc = typeCount[t]
                tcEnd = tc + c
                circos_h.write("%s\t%d\t%d\t%s\t%d\t%d\n" % (t, tc+1, tcEnd, f, intial+1, end))
                typeCount[t] += c          
            intial += c
            




    for t in types:
        tags = TypeTag[t]
        tagNum = len(list(tags))
        print(t, tagNum)

    print("")

    for f in features:
        genes = FeatureGenes[f]
        geneNum = len(list(genes))
        print(f, geneNum)


    print("\n%s" % len(list(GeneSets)))


def main():
    parser = argparse.ArgumentParser(description='Statistics for SV type and gene feature number.')
    parser.add_argument('-a', '--annotation', help='The input annotation file.')
    parser.add_argument('-o', '--out', help='The output file.')
    parser.add_argument("-c", "--circos", help="The out file with circos format.")
    args = parser.parse_args()
    SV_type_gene_feature(args.annotation, args.out, args.circos)

    
if __name__ == '__main__':
    main()



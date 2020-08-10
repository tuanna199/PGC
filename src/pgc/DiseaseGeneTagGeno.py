#!/usr/bin/python
import collections
import sys
import os
import argparse

def disease_gene(omim_enrich):
    """
    omin_enrich
    Gene_set    Term    Overlap P-value Adjusted P-value    Old P-value Old Adjusted P-value    Odds Ratio  Combined Score  Genes
    OMIM_Disease    usher syndrome  4/11    0.02032468808335528 1.0 0   0   3.537318712415989   13.781107075056267  PCDH15;USH1C;MYO7A;USH2A
    OMIM_Disease    lymphoma    6/22    0.020580512259664605    0.9261230516849072  0   0   2.6529890343119917  10.302645890730052  CASP10;RAD54L;ATM;MALT1;MAD1L1;RAD54B
    """
    enrich_h = open(omim_enrich, "r")
    header = enrich_h.readline().strip()

    GeneDisease = collections.defaultdict(set)
    GeneSet = set()
    for line in enrich_h:
        lines = line.strip().split("\t")
        term = lines[1]
        Gene = lines[-1]
        genes = Gene.split(";")
        for g in genes:
            GeneDisease[g].add(term)
            GeneSet.add(g)
    enrich_h.close()
    return GeneSet, GeneDisease



def gene_anno_tag(anno_file):
    """
    anno_file:
    10  53740006    54054127    10_53740006-10_54054127-314121-INV  10  53806578    53807130    PCDH15  ENST00000373965.6  CDS
    10  53740006    54054127    10_53740006-10_54054127-314121-INV  10  53808691    53809546    PCDH15  ENST00000395440.5  CDS
    """
    GeneTags = collections.defaultdict(set)

    anno_h = open(anno_file, "r")
    for line in anno_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        gene = lines[7]
        frature = lines[-1]
        # if frature.lower() == "cds":
        GeneTags[gene].add(tag)
    anno_h.close()
    return GeneTags


def tag_genotype(geno_file):
    """
    geno_file:
    Chr1    Pos1    Chr2    Pos2    SVlength        Type    CN001   CN002   CN003   CN004   CN005   CN006   CN007   CN008   CN009   CN010   CN011   CN01
    1       10059   1       10060   3242    INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    1       10380   1       10381   321     INS     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0     0/0
    """
    TagValues = collections.defaultdict(list)

    TagCount = {}

    tag_h = open(geno_file, "r")
    header = tag_h.readline().strip()
    headers = header.split("\t")
    samples = headers[6:]


    for line in tag_h:
        lines = line.strip().split("\t")
        infor = lines[:6]
        tag = "%s_%s-%s_%s-%s-%s" % tuple(infor)
        genotypes = lines[6:]

        AllCount = 0
        HomoCount = 0
        for s, v in zip(samples, genotypes):
            if v != "0/0":
                value = (s, v)
                TagValues[tag].append(value)
                ### all and homo count
                AllCount += 1
                if v == "1/1":
                    HomoCount += 1
        TagCount[tag] = (AllCount, HomoCount)
    tag_h.close()
    return TagValues, TagCount


def gene_tag_genotype(omim_enrich, anno_file, geno_file, out_file):
    GeneSet, GeneDisease = disease_gene(omim_enrich)
    GeneTags = gene_anno_tag(anno_file)
    TagValues, TagCount = tag_genotype(geno_file)

    out_h = open(out_file, "w")

    genes = sorted(list(GeneSet))
    for g in genes:
        disease = GeneDisease[g]
        Disease = ",".join(sorted(list(disease)))
        if g in GeneTags:
            ### tags is a set
            tags = GeneTags[g]
            sortedTags = sorted(list(tags))
            for t in sortedTags:
                if t in TagValues:
                    ### values is a set of tuple
                    values = TagValues[t]
                    V = [",".join(list(v)) for v in values]
                    VAL = "|".join(V)

                    count = list(TagCount[t])
                    if count[1] == 0:
                        countStr = str(count[0])
                    else:
                        countStr = "%s (%s)" % tuple(count)
                    out_h.write("%s\t%s\t%s\t%s\t%s\n" % (g, Disease, t, countStr, VAL))
                else:
                    print("Please check whether the tag %s is in the genotype file %s." % (t, geno_file))
        else:
            print("Please check whether the gene %s is in annotation file %s." % (g, anno_file))
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Get the sample and genotype for genes.")
    parser.add_argument("-e", "--enrich", help="The enrich file containing term and genes.")
    parser.add_argument("-g", "--genotype", help="The input genotype file containing tag and genotypes.")
    parser.add_argument("-a", "--annotation", help="The input annotation file contain gene and tag.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    gene_tag_genotype(args.enrich, args.annotation, args.genotype, args.out)

if __name__ == "__main__":
    main()



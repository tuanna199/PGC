#!/usr/bin/python
import argparse
import collections
import sys
import os


#usage: python ~/github/NanoHub/src/NanoHub/geneDisease.py --database /home/wuzhikun/Project/Population/population_200624/Annotation/Enrichment/DEL_INS/OMIM_Expand_gene_tag_genotypes.txt,/home/wuzhikun/Project/Population/population_200624/Annotation/Enrichment/DEL_INS/GWAS_gene_tag_genotypes.txt --name OMIM,GWAS --out temp.txt

def database_genes(datas, names, out_file, feature):
    databases = datas.split(",")
    databases = [d.strip() for d in databases]

    ns = names.split(",")
    ns = [n.strip() for n in ns]

    GeneSet = set()
    
    if len(databases) == len(ns):
        DiseaseGenes = collections.defaultdict(set)
        for d, n in zip(databases, ns):
            in_h = open(d, "r")
            for line in in_h:
                lines = line.strip().split("\t")
                gene = lines[0]
                GeneSet.add(gene)
                DiseaseGenes[n].add(gene)
                number = lines[3].split()[0]
            in_h.close()

        GeneDisease = collections.defaultdict(set)
        for d in DiseaseGenes:
            genes = list(DiseaseGenes[d])
            geneNum = len(genes)

            for g in genes:
                GeneDisease[g].add(d)

        diseaseCount = collections.Counter()
        for g in GeneDisease:
            disease = tuple(sorted(list(GeneDisease[g])))
            diseaseCount[disease] += 1

        out_h = open(out_file, "w")
        inital = 0
        ddCont = collections.Counter()
        for d in sorted(list(diseaseCount.keys())):
            dds = sorted(list(d))

            for t in dds:
                c = diseaseCount[d]
                # print(t, ddCont[t]+1, ddCont[t]+c, feature, inital+1, inital+c)
                out_h.write("%s\t%d\t%d\t%s\t%d\t%d\n" % (t, ddCont[t]+1, ddCont[t]+c, feature, inital+1, inital+c))
                ddCont[t] += c
            inital += c


    else:
        print("Please check whether the length of databases and their names is identical.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="statistics of relationshrip of gene to disease.")
    parser.add_argument("-d", "--database", help="The input database of disease.")
    parser.add_argument("-v", "--name", help="The names of databasefile.")
    parser.add_argument("-f", "--feature", default="CDS", help="The feature of gene.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    database_genes(args.database, args.name, args.out, args.feature)


if __name__ == "__main__":
    main()
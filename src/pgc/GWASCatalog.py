#!/usr/bin/python
from BaseFunc import column_index
import collections
import argparse
import sys
import os

#usage: python ~/github/NanoHub/src/NanoHub/GWASCatalog.py --catalog gwas_catalog_v1.0.2-associations_e98_r2020-03-08.tsv --disease gwas_disease_gene.xls --gene gwas_gene_disease.xls

def parse_GWAS_catalog(association_file, pheno_gene, gene_pheno):
    ### get the index of column
    asso_h = open(association_file, "r")
    headers = asso_h.readline().strip().split("\t")
    diseaseIndex = column_index(headers, "DISEASE/TRAIT")
    regionIndex = column_index(headers, "REGION")
    chrIndex = column_index(headers, "CHR_ID")
    posIndex = column_index(headers, "CHR_POS")
    reportedIndex = column_index(headers, "REPORTED GENE(S)")

    DiseaseGenes = collections.defaultdict(set)
    GeneDiseases = collections.defaultdict(set)

    for line in asso_h:
        lines = line.strip().split("\t")
        Disease = lines[diseaseIndex]
        Disease = Disease.split("(")[0].strip()
        Region = lines[regionIndex]
        Chr = lines[chrIndex]
        Pos = lines[posIndex]
        ReportGene = lines[reportedIndex]
        if ReportGene == "NR":
            Genes = []
        else:
            Genes = ReportGene.split(",")
            Genes = [g.strip() for g in Genes]

        ### build the dict for diseases and genes
        if Genes != [] and Genes != [""]:
            for gene in Genes:
                DiseaseGenes[Disease].add(gene)
                GeneDiseases[gene].add(Disease)
    asso_h.close()

    pheno_h = open(pheno_gene, "w")
    pheno_h.write("Disease\tGene\n")
    DiseaseAll = sorted(list(DiseaseGenes.keys()))
    for d in DiseaseAll:
        genes = DiseaseGenes[d]
        pheno_h.write("%s\t%s\n" % (d, ",".join(list(genes))))
    pheno_h.close()

    gene_h = open(gene_pheno, "w")
    gene_h.write("Gene\tDisease\n")
    GeneAll = sorted(list(GeneDiseases.keys()))
    for g in GeneAll:
        disease = GeneDiseases[g]
        gene_h.write("%s\t%s\n" % (g, ",".join(list(disease))))
    gene_h.close()

def main():
    parser = argparse.ArgumentParser(description="Get the gene annotation for translocation.")
    parser.add_argument("-c", "--catalog", help="The input GWAS catalog file.")
    parser.add_argument("-g", "--gene", help="The output file with gene and disease.")
    parser.add_argument("-d", "--disease", help="The output file with disease and gene.")
    args = parser.parse_args()
    parse_GWAS_catalog(args.catalog, args.disease, args.gene)

if __name__ == "__main__":
    main()



#!/usr/bin/python
import collections
import os
import sys
import argparse

#usage: python ~/github/NanoHub/src/NanoHub/TagGeneFeature.py --coding /home/wuzhikun/database/GeneAnno/GRCh38/Annotation/GRCh38_protein_coding_gene.txt  --annotation /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_genepred_overlap.bed --tag /home/wuzhikun/Project/Population/population/bed/database/SV_overlap_filt_tags.xls --out temp.xls --iscoding true

def target_tag(tag_file):
    Tags = {}
    tag_h = open(tag_file, "r")
    for line in tag_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        Tags[tag] = 1
    tag_h.close()
    return Tags


def coding_gene(gene_file):
    """
    gene_file
    $ head  /home/wuzhikun/database/GeneAnno/GRCh38/Annotation/GRCh38_protein_coding_gene.txt
    TSPAN6
    TNMD
    DPM1
    SCYL3
    C1orf112
    """
    CodingGene = []
    gene_h = open(gene_file, "r")
    for line in gene_h:
        line = line.strip()
        if line != "":
            lines = line.split("\t")
            gene = lines[0]
            CodingGene.append(gene)
    return CodingGene



def SV_tag_gene_feature(coding_file, tag_file, annotation_bed, out_file, out_record, iscoding):
    """
    annotation_bed:
    1       88889   88955   1_88889-1_88956-66-INS  1       88550   89550   RP11-34P13.8    ENST00000495576.1       DownStream
    1       90047   90283   1_90047-1_90283-235-DEL 1       89294   91629   RP11-34P13.7    ENST00000466430.5       Exon
    1       90047   90283   1_90047-1_90283-235-DEL 1       89550   90050   RP11-34P13.8    ENST00000495576.1       Exon
    """
    CodingGene = coding_gene(coding_file)

    Tags = target_tag(tag_file)


    GeneRecords = collections.defaultdict(list)
    GeneFeature = collections.defaultdict(set)
    anno_h = open(annotation_bed, "r")

    for line in anno_h:
        line = line.strip()
        lines = line.split("\t")
        tag = lines[3]
        gene = lines[7]
        feature = lines[-1]
        if tag in Tags:
            GeneFeature[gene].add(feature)
            GeneRecords[gene].append(line)
    anno_h.close()

    out_h = open(out_file, "w")
    sortedGenes = sorted(list(GeneFeature.keys()))

    record_h = open(out_record, "w")

    iscoding = iscoding.lower()
    if iscoding == "true":
        for g in sortedGenes:
            if g in CodingGene:
                feature = sorted(list(GeneFeature[g]))
                out_h.write("%s\t%s\n" % (g, ",".join(feature)))

                records = GeneRecords[g]
                record_h.write("%s\n" % "\n".join(records))
    elif iscoding == "false":
        for g in sortedGenes:
            feature = sorted(list(GeneFeature[g]))
            out_h.write("%s\t%s\n" % (g, ",".join(feature)))

            records = GeneRecords[g]
            record_h.write("%s\n" % "\n".join(records))
    else:
        print("Please make sure that the parameter 'iscoding' is 'ture' or 'false'.")
        sys.exit(1)

    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the gene annotation for target tags.")
    parser.add_argument("-c", "--coding", help="The input file containing coding gene.")
    parser.add_argument("-t", "--tag", help="The input bed files which contained tags.")
    parser.add_argument("-a", "--annotation", help="The annotation file for SV tags with bed format.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-r", "--record", help="The record of target genes.")
    parser.add_argument("-s", "--iscoding", default="ture", help="Whether the gene is coding or noncoding, 'true' or 'false'.")
    args = parser.parse_args()
    SV_tag_gene_feature(args.coding, args.tag, args.annotation, args.out, args.record, args.iscoding)

if __name__ == "__main__":
    main()





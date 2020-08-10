#!/usr/bin/python
import collections
import argparse
import sys
import os
import bisect

def ensembl_gene_region(ensembl_genepred):
    """
    ensembl_genepred:
    585     ENST00000619216.1       1       -       17368   17436   17368   17368   1       17368,  17436,  0       MIR6859-1
           none    none    -1,
    585     ENST00000473358.1       1       +       29553   31097   29553   29553   3       29553,30563,30975,      30039,30667,31097,      0       MIR1302-2HG     none    none    -1,-1,-1,
    """
    GeneRegion = {}
    en_h = open(ensembl_genepred, "r")
    for line in en_h:
        lines = line.strip().split("\t")
        Chr, Strand, Start, End, CDSStart, CDSEnd = lines[2:8]
        gene = lines[12]
        Start = int(Start)
        End = int(End)
        CDSStart = int(CDSStart)
        CDSEnd = int(CDSEnd)
        GeneRegion[gene] = (Chr, CDSStart, CDSEnd)
    en_h.close()
    return GeneRegion

def split_tag(tag):
    print(tag)
    tt = tag.split("-")
    Chr, Start = tt[0].split("_")
    End = tt[1].split("_")[-1]
    Start = int(Start)
    End = int(End)
    return Chr, Start, End



def is_tag_cover_whole_gene(tStart, tEnd, rStart, rEnd):
    tRegion = [tStart, tEnd]
    leftIndex = bisect.bisect(tRegion, rStart)
    rightIndex = bisect.bisect(tRegion, rEnd)
    if leftIndex == 1 and rightIndex == 1:
        is_cover = True
    else:
        is_cover = False
    return is_cover




def overlap_gene_filter(ensembl_genepred, gene_overlap, out_file, filter_file):
    """
    gene_overlap:
    1       13102298        13115992        1_13102298-1_13115992-13694-INV 1       13115495        13115518        HNRNPCL2        ENST00000621994.2
       UTR5
    1       13102298        13115992        1_13102298-1_13115992-13694-INV 1       13115518        13116400        HNRNPCL2        ENST00000621994.2
       CDS
    """
    GeneRegion = ensembl_gene_region(ensembl_genepred)

    TagGeneRecord = collections.defaultdict(lambda: collections.defaultdict(list))
    overlap_h = open(gene_overlap, "r")
    out_h = open(out_file, "w")
    for line in overlap_h:
        line = line.strip()
        lines = line.split("\t")
        Chr, Start, End, Tag = lines[:4]
        gene = lines[7]

        SVtype = Tag.split("-")[-1]
        if SVtype == "INS" or SVtype == "DEL":
            out_h.write("%s\n" % line)
        else:
            TagGeneRecord[Tag][gene].append(line)
    overlap_h.close()



    ### filt the gene that the whole region was covered by tag
    filt_h = open(filter_file, "w")
    tags = sorted(list(TagGeneRecord.keys()))
    for t in tags:
        for g in TagGeneRecord[t]:
            tChr, tStart, tEnd = split_tag(t)
            try:
                region = GeneRegion[g]
            except KeyError:
                print("Please check whether the gene %s is in the file %s." % (g, ensembl_genepred))
                # sys.exit(1)


            rChr, rStart, rEnd = region
            print(g, region)

            if rChr == tChr:
                is_cover = is_tag_cover_whole_gene(tStart, tEnd, rStart, rEnd)
                records = TagGeneRecord[t][g]
                print(tStart, tEnd, rStart, rEnd)
                print(is_cover)
                if is_cover == True:
                    filt_h.write("%s\n" % "\n".join(records))
                else:
                    print(records)
                    out_h.write("%s\n" % "\n".join(records))
    out_h.close()
    filt_h.close()




def main():
    parser = argparse.ArgumentParser(description="Filt the record that the whole gene was covered by the INV or DUP.")
    parser.add_argument("-g", "--gene", help="The input genepred file containing gene regions.")
    parser.add_argument("-r", "--overlap", help="The input file with overlaped gene and SV tags.")
    parser.add_argument("-f", "--filt", help="The output file with filter record with INV or DUP covering the whole gene.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    overlap_gene_filter(args.gene, args.overlap, args.out, args.filt)

if __name__ == "__main__":
    main()






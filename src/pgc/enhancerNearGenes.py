#!/usr/bin/python
import argparse
import collections
import bisect

#usage: python ~/github/NanoHub/src/NanoHub/enhancerNearGenes.py --enhancer F5.hg38.enhancers_dechr.bed --genepred /home/wuzhikun/database/Annovar/hg38/hg38_ensGene_gtf95.txt --out temp.txt

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.09.20"



def get_transcript_ragion(gene_pred):
    """
    gene_pred:
    ENST00000456328 1       +       11868   14409   14409   14409   3       11868,12612,13220,      12227,12721,14409,      0       DDX11L1 none    none    -1,-1,-1,
    ENST00000450305 1       +       12009   13670   13670   13670   6       12009,12178,12612,12974,13220,13452,    12057,12227,12697,13052,13374,13670,    0       DDX11L1 none    none    -1,-1,-1,-1,-1,-1,
    ENST00000488147 1       -       14403   29570   29570   29570   11      14403,15004,15795,16606,16857,17232,17605,17914,18267,24737,29533,      14501,15038,15947,16765,17055,17368,17742,18061,18366,24891,29570,      0       WASH7P  none    none    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    """
    ChrRegionRecord = collections.defaultdict(dict)
    gene_h = open(gene_pred, "r")
    for line in gene_h:
        line = line.strip()
        lines = line.split("\t")
        Trans, Chr, Strand, GeneStart, GeneEnd, CDStart, CDEnd, ExonNum, ExonStart, ExonEnd, temp, Gene = lines[:12]
        GeneStart = int(GeneStart)
        GeneEnd = int(GeneEnd)
        # ChrRegionRecord[Chr][(GeneStart, GeneEnd)] = line
        ### get transcript start site and gene name
        if Strand == "+":
            ChrRegionRecord[Chr][(GeneStart, Gene)] = line
        elif Strand == "-":
            ChrRegionRecord[Chr][(GeneEnd, Gene)] = line
        else:
            ChrRegionRecord[Chr][(GeneStart, Gene)] = line
            print("Please check whether the strand of gene is '+' or '-'. for record %s." % line)
    gene_h.close()

    ChrRegions = {}
    for Chrom in ChrRegionRecord:
        sortRegions = sorted(list(ChrRegionRecord[Chrom].keys()))
        Starts = [i[0] for i in sortRegions]
        Genes = [i[1] for i in sortRegions]
        ChrRegions[Chrom] = [Starts, Genes]
    return ChrRegionRecord, ChrRegions


def record_gene(record):
    """
    gene_pred:
    record:
    ENST00000456328 1   +   11868   14409   14409   14409   3   11868,12612,13220,  12227,12721,14409,  0 DDX11L1   none    none    -1,-1,-1,
    ENST00000450305 1   +   12009   13670   13670   13670   6   12009,12178,12612,12974,13220,13452,    12057,12227,12697,13052,13374,13670,    0   DDX11L1 none    none    -1,-1,-1,-1,-1,-1,
    ENST00000488147 1   -   14403   29570   29570   29570   11  14403,15004,15795,16606,16857,17232,17605,17914,18267,24737,29533,  14501,15038,15947,16765,17055,17368,17742,18061,18366,24891,29570,  0   WASH7P  none    none    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,

    """
    records = record.split("\t")
    Trans, Chr, Strand, GeneStart, GeneEnd, CDStart, CDEnd, ExonNum, ExonStart, ExonEnd, temp, Gene = records[:12]
    target = "\t".join([Gene, Trans, Chr, Strand, GeneStart, GeneEnd])
    return target



def enhancer_nearest_gene(gene_pred, enhancer_file, out_file):
    """
    hypothesis: enhancer region is not overlap with gene region

    enhancer_file:
    1   905311  906011  chr1:905311-906011  63  .   905431  905432  0,0,0   2   70,528  0,172
    1   906863  906986  chr1:906863-906986  4   .   906924  906925  0,0,0   2   1,1 0,122
    1   920384  920777  chr1:920384-920777  15  .   920475  920476  0,0,0   2   1,210   0,183


    ### bisect
    a = [2, 4, 5, 8, 12, 16]
    b = bisect.bisect(a, 10) -> 4
    b = bisect.bisect(a, 8) ->4

    out_file:
    1   905311  906011  AL645608.6  ENST00000607769 1   +   904833  915976
    1   905311  906011  AL645608.2  ENST00000448179 1   +   911434  914948
    1   906863  906986  AL645608.6  ENST00000607769 1   +   904833  915976
    1   906863  906986  AL645608.2  ENST00000448179 1   +   911434  914948
    """
    ChrRegionRecord, ChrRegions = get_transcript_ragion(gene_pred)

    out_h = open(out_file, "w")
    in_h = open(enhancer_file, "r")    
    for line in in_h:
        line = line.strip()
        lines = line.split("\t")
        Chr, Start, End = lines[:3]
        Start = int(Start)
        End = int(End)
        if Chr in ChrRegions:

            targetIndex = []

            STARTS, Genes = ChrRegions[Chr]
            startIndex = bisect.bisect(STARTS, Start)
            endIndex = bisect.bisect(STARTS, End)
            if startIndex == endIndex:
                ### smaller than first one
                if startIndex == -1:
                    targetIndex.append(startIndex + 1)
                ### larger than final one
                elif startIndex == len(STARTS):
                    targetIndex.append(startIndex - 1)
                else:
                    targetIndex.append(startIndex - 1)
                    targetIndex.append(startIndex)
            elif endIndex - startIndex >= 1:
                if startIndex <= 0:
                    newStart = 0
                else:
                    newStart = startIndex - 1

                if endIndex == len(STARTS):
                    newEnd = endIndex
                else:
                    newEnd = endIndex + 1

                for i in range(newStart, newEnd):
                    targetIndex.append(i)
                print("Please check whether the record %s in file enhancer file %s overlap with multiple regions of file %s." % (line, enhancer_file, gene_pred))

            ### get the region based on index
            if len(targetIndex) >= 1:
                for t in targetIndex:
                    startGene = (STARTS[t], Genes[t])
                    record = ChrRegionRecord[Chr][startGene]
                    targetInfor = record_gene(record)
                    out_h.write("%s\t%s\t%s\t%s\n" % (Chr, Start, End, targetInfor))
    in_h.close()
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the nearest genes for enhancer.")
    parser.add_argument("-e", "--enhancer", help="The file with enhancer region.")
    parser.add_argument("-g", "--genepred", help="The input gene prediction file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    enhancer_nearest_gene(args.genepred, args.enhancer, args.out)

if __name__ == "__main__":
    main()





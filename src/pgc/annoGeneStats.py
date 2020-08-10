#!/usr/bin/python
import collections
import argparse
import sys

#usage: python ~/github/NanoHub/src/NanoHub/annoGeneStats.py --stats Samples_recessive_gene_stat.xls --location location.txt --gene genes.txt


def anno_gene_stats(anno_gene, loc_file, gene_file):
    """
    AnnotSV ID  SV chrom    SV start    SV end  SV length   SV type ID  REF ALT QUAL    FILTER  INFO    FORMAT  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-0.bam  /
    home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-1.bam  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-2.bam AnnotSV type    Gene name   NM  CDS length  tx length    location   intersectStart  intersectEnd    DGV_GAIN_IDs    DGV_GAIN_n_samples_with_SV  DGV_GAIN_n_samples_tested   DGV_GAIN_Frequency  DGV_LOSS_IDs    DGV_LOSS_n_samples_with_SV  DGV_LOSS_n_samples_tested   DGV_LOSS_Frequency  1000g_event 1000g_AF    1000g_max_AF    IMH_ID  IMH_AF  IMH_ID_others   promoters   dbVar_event dbVar_variant   dbVar_status    TADcoordinates  ENCODEexperiments   GCcontent_left  GCcontent_right Repeats_coord_left  Repeats_type_left   Repeats_coord_right Repeats_type_right  ACMG    HI_CGscore  TriS_CGscore    DDD_status  DDD_mode    DDD_consequence DDD_disease DDD_pmids   HI_DDDpercent   synZ_ExAC   misZ_ExAC   pLI_ExAC    delZ_ExAC   dupZ_ExAC   cnvZ_ExAC    morbidGenes    morbidGenesCandidates   Mim Number  Phenotypes  Inheritance AnnotSV ranking
    1_6937556_6937752_INS    1   6937556 6937752 196 INS INS00359SUR N   <INS>   .   PASS    SUPP=3;SUPP_VEC=111;AVGLEN=175;SVTYPE=INS;SVMETHOD=SURVIVORv2;CHR2=1;END=6937752;CIPOS=0,24;CIEND=-2,14;STRANDS=+-    GT:PSV:LN:DR:ST:TY:CO   1/1:NA:170:0,16:+-:INS:1_6937580-1_6937766  0/1:NA:166:11,12:+-:INS:1_6937556-1_6937752 0/1:NA:189:15,18:+-:INS:1_6937556-1_6937750  split  CAMTA1  NM_001349608    0   197 intron2-intron2 6937556 6937752     0   0   -1      0   0   -1      -1  -1      -1           probable   monoallelic loss of function    CEREBELLAR ATAXIA, NONPROGRESSIVE, WITH MENTAL RETARDATION  22693284    2.22    2.41352520512185    3.71208472645534    0.999999946187373   0.837784888692585   1.39477205262828    1.3883934369361 yes     611501  Cerebellar ataxia, nonprogressive, with mental retardation, 614756 (3)  AD  4
    11_59209875_59209943_INS   11  59209875    59209943    68  INS INS0020551SUR   N   <INS>   .   PASS    SUPP=3;SUPP_VEC=111;AVGLEN=69.6667;SVTYPE=INS;SVMETHOD=SURVIVO
    Rv2;CHR2=11;END=59209943;CIPOS=0,15;CIEND=0,28;STRANDS=+-   GT:PSV:LN:DR:ST:TY:CO   1/1:NA:69:0,11:+-:INS:11_59209877-11_59209944   0/1:NA:68:16,8:+-:INS:11_59209875-11_59209943   0/1:NA:72:27,14:+-:INS:11_59209890-11_59209971  split   MPEG1   NM_001039396    0   69  exon1-exon1 59209875    59209943        0   0   -1      0   0   -1   -1 -1      -1                                                                                       67.82  -1.31912548677862   -1.18139249142939   1.39836996659703e-11    0.26081377275072    0.410662973826508   0.514646997097845                       2
    """
    LocationNum = collections.Counter()
    LocationGene = collections.defaultdict(list)

    anno_h = open(anno_gene, "r")
    for line in anno_h:
        lines = line.strip().split("\t")
        location, number, gene = lines[1:4]
        number = int(number)
        LocationNum[location] += number
        LocationGene[location].append(gene)
    anno_h.close()

    LocAllgenes = {}
    ### output the file
    loc_h = open(loc_file, "w")
    loc_h.write("Location\tNumber\tGene\n")
    locations = sorted(list(LocationNum.keys()))
    for loc in locations:
        number = LocationNum[loc]
        try:
            genes = LocationGene[loc]
            allgenes = sorted(",".join(genes).split(","))
            newGenes = ",".join(allgenes)
            loc_h.write("%s\t%d\t%s\n" % (loc, number, newGenes))
            ### build new dict
            LocAllgenes[loc] = allgenes

        except KeyError:
            print("Please check whether the location %s is in the file %s." % (loc, anno_gene))
            sys.exit(1)
    loc_h.close()

    GeneLocNum = collections.defaultdict(lambda:collections.Counter())
    ### output genes
    gene_h = open(gene_file, "w")
    UniqueGenes = set()
    for loc in locations:
        allgenes = LocAllgenes[loc]
        if loc == "Intron":
            for ge in allgenes:
                UniqueGenes.add(ge)
                GeneLocNum["Intron"][ge] += 1 
        else:
            for ge in allgenes:
                UniqueGenes.add(ge)
                GeneLocNum["Exon"][ge] += 1

    UniqueGenesList = sorted(list(UniqueGenes))
    for ge in UniqueGenesList:
        if ge in GeneLocNum["Exon"]:
            geExon = GeneLocNum["Exon"][ge]
        else:
            geExon = 0

        if ge in GeneLocNum["Intron"]:
            geIntron = GeneLocNum["Intron"][ge]
        else:
            geIntron = 0

        num = geExon + geIntron

        gene_h.write("%s\t%d\t%d\t%d\n" % (ge, geExon, geIntron, num))

    gene_h.close()        



def main():
    parser = argparse.ArgumentParser(description="Stats for gene number and genes for location.")
    parser.add_argument("-s", "--stats", help="The input stats for annotation file.")
    parser.add_argument("-l", "--location", help="Output the location and the genes ")
    parser.add_argument("-g", "--gene", help="Output the gene and location.")
    args = parser.parse_args()
    anno_gene_stats(args.stats, args.location, args.gene)

if __name__ == "__main__":
    main()

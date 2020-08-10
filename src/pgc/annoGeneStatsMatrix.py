#!/usr/bin/python
import collections
import argparse
import sys
import operator

#usage: python ~/github/NanoHub/src/NanoHub/annoGeneStatsMatrix.py --stats Samples_recessive_gene_stat.xls --location temp_loc --gene temp_gene --convert True


def sort_dict_value(keyValue):
    sortedValue = sorted(keyValue.items(), key=operator.itemgetter(1), reverse=True)
    return sortedValue



def anno_gene_stats(anno_gene, loc_file, gene_file, isConvert):
    """
    AnnotSV ID  SV chrom    SV start    SV end  SV length   SV type ID  REF ALT QUAL    FILTER  INFO    FORMAT  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-0.bam  /
    home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-1.bam  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M671-2.bam AnnotSV type    Gene name   NM  CDS length  tx length    location   intersectStart  intersectEnd    DGV_GAIN_IDs    DGV_GAIN_n_samples_with_SV  DGV_GAIN_n_samples_tested   DGV_GAIN_Frequency  DGV_LOSS_IDs    DGV_LOSS_n_samples_with_SV  DGV_LOSS_n_samples_tested   DGV_LOSS_Frequency  1000g_event 1000g_AF    1000g_max_AF    IMH_ID  IMH_AF  IMH_ID_others   promoters   dbVar_event dbVar_variant   dbVar_status    TADcoordinates  ENCODEexperiments   GCcontent_left  GCcontent_right Repeats_coord_left  Repeats_type_left   Repeats_coord_right Repeats_type_right  ACMG    HI_CGscore  TriS_CGscore    DDD_status  DDD_mode    DDD_consequence DDD_disease DDD_pmids   HI_DDDpercent   synZ_ExAC   misZ_ExAC   pLI_ExAC    delZ_ExAC   dupZ_ExAC   cnvZ_ExAC    morbidGenes    morbidGenesCandidates   Mim Number  Phenotypes  Inheritance AnnotSV ranking
    1_6937556_6937752_INS    1   6937556 6937752 196 INS INS00359SUR N   <INS>   .   PASS    SUPP=3;SUPP_VEC=111;AVGLEN=175;SVTYPE=INS;SVMETHOD=SURVIVORv2;CHR2=1;END=6937752;CIPOS=0,24;CIEND=-2,14;STRANDS=+-    GT:PSV:LN:DR:ST:TY:CO   1/1:NA:170:0,16:+-:INS:1_6937580-1_6937766  0/1:NA:166:11,12:+-:INS:1_6937556-1_6937752 0/1:NA:189:15,18:+-:INS:1_6937556-1_6937750  split  CAMTA1  NM_001349608    0   197 intron2-intron2 6937556 6937752     0   0   -1      0   0   -1      -1  -1      -1           probable   monoallelic loss of function    CEREBELLAR ATAXIA, NONPROGRESSIVE, WITH MENTAL RETARDATION  22693284    2.22    2.41352520512185    3.71208472645534    0.999999946187373   0.837784888692585   1.39477205262828    1.3883934369361 yes     611501  Cerebellar ataxia, nonprogressive, with mental retardation, 614756 (3)  AD  4
    11_59209875_59209943_INS   11  59209875    59209943    68  INS INS0020551SUR   N   <INS>   .   PASS    SUPP=3;SUPP_VEC=111;AVGLEN=69.6667;SVTYPE=INS;SVMETHOD=SURVIVO
    Rv2;CHR2=11;END=59209943;CIPOS=0,15;CIEND=0,28;STRANDS=+-   GT:PSV:LN:DR:ST:TY:CO   1/1:NA:69:0,11:+-:INS:11_59209877-11_59209944   0/1:NA:68:16,8:+-:INS:11_59209875-11_59209943   0/1:NA:72:27,14:+-:INS:11_59209890-11_59209971  split   MPEG1   NM_001039396    0   69  exon1-exon1 59209875    59209943        0   0   -1      0   0   -1   -1 -1      -1                                                                                       67.82  -1.31912548677862   -1.18139249142939   1.39836996659703e-11    0.26081377275072    0.410662973826508   0.514646997097845                       2
    """
    LocationNum = collections.Counter()
    LocationGene = collections.defaultdict(list)


    GeneCatSample = collections.defaultdict(lambda: collections.defaultdict(list))
    CatGeneSample = collections.defaultdict(lambda: collections.defaultdict(list))

    allLocations = set()
    anno_h = open(anno_gene, "r")
    for line in anno_h:
        lines = line.strip().split("\t")
        sample, location, number, gene = lines[:4]
        number = int(number)

        ### whether convert the category to "Exon" or "Intron"
        if isConvert == "True":
            if location == "Intron":
                newLoc = "Intron"
            else:
                newLoc = "Exon"
        elif isConvert == "False":
            newLoc = location
        else:
            print("Please check whether convert the original category to 'Intron' or 'Exon' based on True of False.")
            sys.exit(1)

        allLocations.add(newLoc)
        ### get the dict of gene -> location -> sample
        genes = gene.split(",")
        for g in genes:
            GeneCatSample[g][newLoc].append(sample)

            ### get the location -> gene -> sample
            CatGeneSample[newLoc][g].append(sample)
    anno_h.close()


    ## output gene and number in samples
    ### sort all locations
    sortedAllLocation = sorted(list(allLocations))

    gene_h = open(gene_file, "w")

    headerSample = [l + "_samples" for l in sortedAllLocation]
    gene_h.write("Gene\tTotal\t%s\t%s\n" % ("\t".join(sortedAllLocation), "\t".join(headerSample)))

    GeneRecord = {}
    GeneNumber = {}

    allGenes = sorted(list(GeneCatSample.keys()))
    for ge in allGenes:
        ### get the number and samples for each location of each gene
        GeneNum = []
        GeneSample = []

        for loc in sortedAllLocation:
            if loc in GeneCatSample[ge]:
                samples = GeneCatSample[ge][loc]
                ##############################
                ####### unique for samples
                samples = sorted(list(set(samples)))
                sampleNum = len(samples)
            else:
                sampleNum = 0
                samples = ["-"]

            GeneNum.append(sampleNum)
            GeneSample.append(samples)

        GeneNumSum = sum(GeneNum)
        CatNumOut = "\t".join([str(g) for g in GeneNum])
        CatSampleOut = "\t".join([",".join(s) for s in GeneSample])

        record = "%s\t%d\t%s\t%s\t" % (ge, GeneNumSum, CatNumOut,  CatSampleOut)
        GeneNumber[ge] = GeneNumSum
        GeneRecord[ge] = record
    
    ### output
    GeneNumSorted = sort_dict_value(GeneNumber)
    for g, n in GeneNumSorted:
        r = GeneRecord[g]
        gene_h.write("%s\n" % r)

    gene_h.close()        


    ### location and genes
    loc_h = open(loc_file, "w")
    loc_h.write("Location\tGeneNumber\tGenes\tSampleNumber\tSamples\n")
    for loc in sortedAllLocation:
        geneSample = CatGeneSample[loc]
        genes = sorted(list(geneSample.keys()))
        geneNum = len(genes)
        samNum = 0
        samList = []
        for ge in geneSample:
            sam = geneSample[ge]
            samList.append(sam)
            samNum += len(sam)
        samOut = ";".join([",".join(s) for s in samList])
        loc_h.write("%s\t%d\t%s\t%d\t%s\n" % (loc, geneNum, ",".join(genes), samNum, samOut))
    loc_h.close()




def main():
    parser = argparse.ArgumentParser(description="Stats for gene number and genes for location.")
    parser.add_argument("-s", "--stats", help="The input stats for annotation file.")
    parser.add_argument("-l", "--location", help="Output the location and the genes ")
    parser.add_argument("-g", "--gene", help="Output the gene and location.")
    parser.add_argument("-c", "--convert", default="True", help="Whether convert the category.")
    args = parser.parse_args()
    anno_gene_stats(args.stats, args.location, args.gene, args.convert)

if __name__ == "__main__":
    main()

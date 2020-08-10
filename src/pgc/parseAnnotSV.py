#!/usr/bin/python
import collections
import argparse
import sys
import re

#usage: python ~/github/NanoHub/src/NanoHub/parseAnnotSV.py --annotation M416_recessive_SV.tsv --out  temp.txt


def target_index(headers, column):
    try:
        targetIndex = headers.index(column.lower())
        return targetIndex
    except ValueError:
        print("Please check whether the target column %s is in the file %s." % (column, anno_file))
        sys.exit(1)


def match_location(location):
    """
    a = "intron27"
    b = re.findall("(\D*)(\d*)", a)
    b = [('intron', '27'), ('', '')]

    a = "txStart"
    b = re.findall("(\D*)(\d*)", a)
    b = [('txStart', ''), ('', '')]
    """
    locations = location.split("-")
    locLen = len(locations)
    if locLen == 2:
        loc1, loc2 = locations
        locMatch1 = re.findall("(\D*)(\d*)", loc1)
        locMatch2 = re.findall("(\D*)(\d*)", loc2)
        locNum1 = locMatch1[0][1]
        locNum2 = locMatch2[0][1]
        if locNum1 == '' and locNum2 == '':
            newloc1 = loc1.lstrip("tx")
            newloc2 = loc2.lstrip("tx")
        elif locNum1 == '' and locNum2 != '':
            newloc1 = loc1.lstrip("tx")
            newloc2 = locMatch2[0][0]
        elif locNum1 != '' and locNum2 == '':
            newloc1 = locMatch1[0][0]
            newloc2 = loc2.lstrip("tx")
        elif locNum1 != '' and locNum2 != '':
            newloc1 = locMatch1[0][0]
            newloc2 = locMatch2[0][0]
        return newloc1, newloc2
    else:
        print("Please check whether the location of %s has start and end positions." % location)
        sys.exit(1)



def parse_annotSV(anno_file, out_file):
    """
    anno_file:
    AnnotSV ID  SV chrom    SV start    SV end  SV length   SV type ID  REF ALT QUAL    FILTER  INFO    FORMAT  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M416-0.bam  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M416-1.bam /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M416-2.bam AnnotSV type    Gene name   NM  CDS length  tx length   location    intersectStart  intersectEnd    DGV_GAIN_IDs    DGV_GAIN_n_samples_with_SV  DGV_GAIN_n_samples_tested   DGV_GAIN_Frequency  DGV_LOSS_IDs    DGV_LOSS_n_samples_with_SV  DGV_LOSS_n_samples_tested   DGV_LOSS_Frequency  1000g_event 1000g_AF    1000g_max_AF    IMH_ID  IMH_AF  IMH_ID_otherpromoters   dbVar_event dbVar_variant   dbVar_status    TADcoordinates  ENCODEexperiments   GCcontent_left  GCcontent_right Repeats_coord_left  Repeats_type_left   Repeats_coord_right Repeats_type_right  ACMG    HI_CGscore  TriS_CGscore    DDD_status  DDD_mode    DDD_consequence DDD_disease DDD_pmids   HI_DDDpercensynZ_ExAC   misZ_ExAC   pLI_ExAC    delZ_ExAC   dupZ_ExAC   cnvZ_ExAC   morbidGenes morbidGenesCandidates   Mim Number  Phenotypes  Inheritance AnnotSV ranking
    1_6937571_6937738_INS    1   6937571 6937738 167 INS INS00630SUR N   <INS>   .   PASS    SUPP=3;SUPP_VEC=111;AVGLEN=176;SVTYPE=INS;SVMETHOD=SURVIVORv2;CHR2=1;END=6937738;CIPOS=-26,0;CIEND=0,18;STRANDS=+-   GT:PSV:LN:DR:ST:TY:CO   1/1:NA:186:0,13:+-:INS:1_6937561-1_6937756  0/1:NA:163:13,14:+-:INS:1_6937571-1_6937738 0/1:NA:179:12,11:+-:INS:1_6937545-1_6937747 split   CAMTA1  NM_001349608    0   168 intron2-intron2 6937571 6937738     0   0   -1      0   0   -1      -1  -1      -1                                                          3   0   probable    monoallelic loss of function    CEREBELLAR ATAXIA, NONPROGRESSIVE, WITH MENTAL RETARDATION  22693284    2.22    2.41352520512185    3.71208472645534    0.999999946187373   0.837784888692585   1.39477205262828    1.3883934369361 yes     611501  Cerebellar ataxia, nonprogressive, with mental retardation, 614756 (3)  AD  4
    1_7115490_7115640_DEL  1   7115490 7115640 150 DEL DEL00640SUR N   <DEL>   .   PASS    SUPP=3;SUPP_VEC=111;AVGLEN=-147.667;SVTYPE=DEL;SVMETHOD=SURVIVORv2;CHR2=1;END=7115640;CIPOS=0,21;CIEND=0,16;STRANDS=+-   GT:PSV:LN:DR:ST:TY:CO   1/1:NA:145:0,11:+-:DEL:1_7115511-1_7115656  0/1:NA:150:16,11:+-:DEL:1_7115490-1_7115640 0/1:NA:148:5,6:+-:DEL:1_7115492-1_7115640   split   CAMTA1  NM_001349608    0   151 intron3-intron3 7115490 7115640     0   0   -1      0   0   -1  -1  -1      -1  

    out_file:
    Location        Number  Gene
    Exon-Intron     1       KMT5A
    Intron  40      ADAMTS2,C11orf49,CAMTA1,CAMTA1,CCSER1,CERK,CHRNA3,COL6A1,CSMD1,DIP2C,DLGAP1,DLGAP1-AS4,DPP6,EEF1DP3,FAM155A,FAM49B,GABRB3,GNS,GOLGA2P10,KHDRBS2,LINC00504,LOC727751,L
    OC729218,LRRC4C,MAP3K12,MXRA7,NRP1,PLCL2,RAPGEF5,RNF175,SEMA6D,STPG4,TKT,TMEM132B,TMEM132B,TSHZ2,TULP4,UBE2E3,WAPL,ZMYND8Start-End       7       CCNYL3,FRG2DP,LINC01566,LINC02167,LOC112268173,TP53TG3HP,UBE2MP1
    """
    ### get sample name
    sampleName = anno_file.split("/")[-1].split("_")[0]


    RegionNumber = collections.Counter()
    RegionGenes = collections.defaultdict(list)

    anno_h = open(anno_file, "r")
    headers = anno_h.readline().strip().split("\t")
    headers = [h.lower() for h in headers]

    ### get gene name, transcript, location
    geneIndex = target_index(headers, "Gene name")
    transIndex = target_index(headers, "NM")
    locationIndex = target_index(headers, "location")

    for line in anno_h:
        lines = line.strip().split("\t")
        SVID, Chr, Start, End, SVLen, SVType = lines[:6]
        Gene = lines[geneIndex]
        Trans = lines[transIndex]
        Loc = lines[locationIndex]
        loc1, loc2 = match_location(Loc)
        if loc1 == loc2:
            region = loc1.capitalize()
        else:
            region = loc1.capitalize() + "-" + loc2.capitalize()

        RegionNumber[region] += 1
        RegionGenes[region].append(Gene)
    anno_h.close()

    ### output the file
    out_h = open(out_file, "w")
    # out_h.write("Location\tNumber\tGene\n")
    for r in sorted(list(RegionNumber.keys())):
        number = RegionNumber[r]
        Genes = RegionGenes[r]
        out_h.write("%s\t%s\t%d\t%s\n" % (sampleName, r, number, ",".join(sorted(Genes))))
    out_h.close()



def main():
    parser = argparse.ArgumentParser(description="Get the annotation information for SV.")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    parse_annotSV(args.annotation, args.out)

if __name__ == "__main__":
    main()







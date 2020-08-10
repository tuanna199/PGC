#!/usr/bin/python
from BaseFunc import parse_genotype_format
from BaseFunc import Infor_target_values
import sys
import argparse

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2019.09.09"

#usage: python ~/github/NanoHub/src/NanoHub/annotsvIDModify.py --annotation /home/wuzhikun/Project/NanoTrio/Recessive/AnnotSV/M470/M470_recessive_SV.tsv  --out temp_anno.txt

#usage: python ~/github/NanoHub/src/NanoHub/annotsvIDModify.py --annotation /home/wuzhikun/Project/NanoTrio/Recessive/AnnotSV/M470/M470_recessive_SV.tsv --proband M470-0 --column 14,15,16 --out temp_anno.txt 

def modify_annotSV_ID(annotsv_file, out_file): # proband, column
    """
    annotsv_file:
    AnnotSV ID      SV chrom        SV start        SV end  SV length       SV type ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M470-0.bam     /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M470-1.bam     /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M470-2.bam     AnnotSV type    Gene name       NM      CDS length      tx length       location        intersectStart  intersectEnd    DGV_GAIN_IDs    DGV_GAIN_n_samples_with_SV      DGV_GAIN_n_samples_tested       DGV_GAIN_Frequency      DGV_LOSS_IDs    DGV_LOSS_n_samples_with_SV      DGV_LOSS_n_samples_tested       DGV_LOSS_Frequency      1000g_event     1000g_AF        1000g_max_AF    IMH_ID  IMH_AF  IMH_ID_others   promoters       dbVar_event     dbVar_variant   dbVar_status    TADcoordinates  ENCODEexperiments       GCcontent_left  GCcontent_right Repeats_coord_left      Repeats_type_left       Repeats_coord_right     Repeats_type_right      ACMG    HI_CGscore      TriS_CGscore    DDD_status      DDD_mode        DDD_consequence DDD_disease     DDD_pmids       HI_DDDpercent   synZ_ExAC       misZ_ExAC       pLI_ExAC        delZ_ExAC       dupZ_ExAC       cnvZ_ExAC       morbidGenes     morbidGenesCandidates   Mim Number      Phenotypes      Inheritance     AnnotSV ranking
    21_10458143_10458144_BND        21      10458143        10458144        1       BND     TRA0038486SUR   N       N[Y:56854863[   .       PASS    SUPP=3;SUPP_VEC=111;AVGLEN=100000;SVTYPE=BND;SVMETHOD=SURVIVORv2;CHR2=Y;END=56854863;CIPOS=0,0;CIEND=0,0;STRANDS=+-     GT:PSV:LN:DR:ST:TY:CO   1/1:NA:46396720:0,29:+-:TRA:21_10458143-Y_56854863      0/1:NA:46396720:30,15:+-:TRA:21_10458143-Y_56854863     0/1:NA:46396720:10,10:+-:TRA:21_10458143-Y_56854863     split   BAGE4   NM_181704       0       2       intron3-intron3 10458143        10458144        gssvG23344,gssvG23341,gssvG23339,gssvG23345,gssvG23342,gssvG23365,gssvG23354,gssvG23360 345     20897    0.01650955     gssvL74302,gssvL74298,gssvL74297,gssvL74305,gssvL74300  230     21917    0.01049414             -1      -1              -1
    """


    out_h = open(out_file, "w")
    anno_h = open(annotsv_file, "r")
    header = anno_h.readline().strip()
    out_h.write("%s\n" % header)

    # columns = column.split(",")
    # columns = [int(c) for c in columns]

    # headers = header.split("\t")
    # ### get target sample
    # indexColumns = [c-1 for c in columns]

    # sampleIndex = ""
    # for s in indexColumns:
    #     sample = headers[s]
    #     if proband in sample:
    #         sampleIndex = headers.index(sample)
    #         break

    # if sampleIndex == "":
    #     print("Please check whether the target proband name %s is in the header of file %s." % (proband, annotsv_file))
    #     sys.exit(1)

    ### deal with the record
    for line in anno_h:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        else:
            lines = line.split("\t")
            ID = lines[0]
            Chr, Start = lines[1:3]
            Infor = lines[11]

            # Format = lines[12]
            # ### change the ID by the Tag
            # sampleRecord = lines[sampleIndex]
            # ## TY: TRA
            # TY = parse_genotype_format(Format, sampleRecord, "TY")
            # ## CO: 21_10458143-Y_56854863
            # CO = parse_genotype_format(Format, sampleRecord, "CO")
            # if "," in TY and "," in CO:
            #     Tag = CO.split(",")[0] + "_" + TY.split(",")[0]
            # else:
            #     Tag = CO + "_" + TY

            ### change the ID based on the Infor
            Chr2 = Infor_target_values(Infor, "CHR2")
            End = Infor_target_values(Infor, "END")
            SVType = Infor_target_values(Infor, "SVTYPE")
            if "SVLEN" in Infor:
                SVLength = Infor_target_values(Infor, "SVLEN")
            elif "AVGLEN" in Infor:
                SVLength = Infor_target_values(Infor, "AVGLEN")
                
            if SVLength.startswith("-"):
                SVLength = SVLength.lstrip("-")


            if SVType == "BND":
                SVType = "TRA"

            if SVType == "TRA":
                SVLength = 0
            # Tag = "%s_%s-%s_%s_%s" % (Chr, Start, Chr2, End, SVType)
            # Tag = "%s_%s_%s_%s_%s" % (Chr, Start, Chr2, End, SVType) ### two chromosomes
            Tag = "%s_%s-%s_%s-%s-%s" % (Chr, Start, Chr2, End, SVLength, SVType) ### just one chromosomes
            out_h.write("%s\t%s\n" % (Tag, "\t".join(lines[1:])))

    anno_h.close()
    out_h.close()

def main():
    parser = argparse.ArgumentParser(description="Modify the ID for annotation file derived from annotSV.")
    parser.add_argument("-a", "--annotation", help="The annotation file from annotSV.")
    # parser.add_argument("-p", "--proband", help="The given proband name.")
    # parser.add_argument("-c", "--column", default="14,15,16", help="The given column names for samples.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    modify_annotSV_ID(args.annotation, args.out)

if __name__ == "__main__":
    main()


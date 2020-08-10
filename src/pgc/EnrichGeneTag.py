#!/usr/bin/python
import collections
import argparse
import sys
import os

def gene_record(anno_file):
    """
    anno_file:
    1       67910   68341   1_67910-1_68341-433-DEL 1       68090   69090   OR4F5   ENST00000335137.3       UpStream
    1       88684   88831   1_88684-1_88831-147-DEL 1       88294   89294   RP11-34P13.7    ENST00000466430.5       DownStream
    1       88684   88831   1_88684-1_88831-147-DEL 1       88550   89550   RP11-34P13.8    ENST00000495576.1       DownStream
    """
    GeneRecord = collections.defaultdict(list)
    gene_h = open(anno_file, "r")
    for line in gene_h:
        line = line.strip()
        lines = line.split("\t")
        gene = lines[7]
        GeneRecord[gene].append(line)
    gene_h.close()

    return GeneRecord



def gene_tags(anno_file, enrich_file, out_file, target_term):
    """
    enrich_file:
    Gene_set        Term    Overlap P-value Adjusted P-value        Old P-value     Old Adjusted P-value    Odds Ratio      Combined Score  Genes
    GWAS_Catalog_2019       Post bronchodilator FEV1/FVC ratio      132/206 1.6579505371686225e-06  0.002879860083061897    0       0       1.340677265465271
           17.84431832307067       RB1;CLSTN2;ZFAND3;ULK4;FAM13A;CEP128;GRIP1;HERC3;PPP4R4;RUVBL1;FAM19A2;KPNA4;PAMR1;MAGI1;PDGFRA;SPINK2;DAPK1;KCNH7;DNAAF1;SEMA6D;THEMIS;AFAP1;C2ORF16;MYO7A;NPNT;STON1-GTF2A1L;MSH2;KIF16B;RIN3;FBLN7;NUP205;CHRNA3;IREB2;NMD3;CACNA1C;KALRN;RLF;PKHD1;FUT8;TSPAN8;OTOG;MYL10;LRRC4C;RAB6B;EEFSEC;TGFB2;GCKR;WDR59;CADM2;ARHGEF38;ESRRG;BTBD9;GULP1;PTPRD;DLG2;GALNTL5;TCF4;CNTN4;KIAA0825;FGF12;ROBO2;DOCK5;SETD2;SERPINA1;MCTP1;SLC35F3;MAST4;KCNC2;GRIK3;CELF4;DBH;LFNG;IQCJ-SCHIP1;MTG1;ZNF385D;ERC2;DLGAP2;EDIL3;RALGPS1;ADAMTS7;PACRG;GRID2;PDE4D;TET2;COBL;KSR2;ANO4;IQCJ;SUPT3H;PARD3B;CYP2A7;CYP2A6;NPEPL1;RARB;CDH13;MAPRE3;ARHGEF5;CNTNAP5;MCPH1;LILRA6;OCA2;KANK1;SP100;DDX24;PPM1L;NRXN1;PALM2;HEATR1;ACVR1B;PPM1E;PAK1;GNG2;WIF1;PLAGL1;MGAT5;UBTD1;SPOCK1;ZPLD1;WWOX;R3HDML;EGLN2;PNPLA8;MOCOS;KCNIP4;MGAT4C;SLCO3A1;VWC2;SARDH;MSRB3;BUD13;MDGA2;HCN1
    """
    GeneRecord = gene_record(anno_file)

    enrich_h = open(enrich_file, "r")
    header = enrich_h.readline()
    out_h = open(out_file, "w")

    target_term = target_term.strip('"').lower()
    if target_term == "all":
        target = " "
    else:
        target = target_term

    for line in enrich_h:
        lines = line.strip().split("\t")
        Term = lines[1]
        Term = Term.lower()

        Gene = lines[-1]
        genes = Gene.split(";")

        if target in Term:
            for g in genes:
                records = []
                if g in GeneRecord:
                    record = GeneRecord[g]
                    records.append(record)
                else:
                    print("Please check whether the gene %s is in annotation file %s." % (g, anno_file))


                for r in records:
                    for rr in r:
                        out_h.write("%s\t%s\n" % (g, rr))
    enrich_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the target record and annotation for enriched genes.")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    parser.add_argument("-e", "--enrichment", help="The input enrichment file.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-t", "--term", help="The target term, if 'all', all terms were selected.")
    args = parser.parse_args()
    gene_tags(args.annotation, args.enrichment, args.out, args.term)

if __name__ == "__main__":
    main()


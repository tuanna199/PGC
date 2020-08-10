#!/usr/bin/python
import gseapy
import argparse
import os
import sys

#usage: python /home/wuzhikun/github/NanoHub/src/NanoHub/geneGSEA.py --annotation /home/wuzhikun/Project/NanoTrio/Denovo/Annotation/Genes_annotation_region_IDs.xls --outdir /home/wuzhikun/Project/NanoTrio/Denovo/Annotation/Enrichment --library KEGG_2019_Human,GO_Biological_Process_2018,GO_Molecular_Function_2018,GO_Cellular_Component_2018,GWAS_Catalog_2019,UK_Biobank_GWAS_v1,MGI_Mammalian_Phenotype_Level_4_2019,OMIM_Disease,OMIM_Expanded,Reactome_2016 --method website

## or

#usage: python ~/github/NanoHub/src/NanoHub/geneGSEA.py --method local --annotation /home/wuzhikun/Project/NanoTrio/Denovo/Annotation/Samples_denovo_annotation_genes.xls --outdir /home/wuzhikun/Project/NanoTrio/Denovo/Annotation/intolerant


#>>> gseapy.enrichr(gene_list=gene, description='pathway', gene_sets='KEGG_2019_Human', outdir="testGSEA")

"""
### library
['ARCHS4_Cell-lines', 'ARCHS4_IDG_Coexp', 'ARCHS4_Kinases_Coexp', 'ARCHS4_TFs_Coexp', 'ARCHS4_Tissues', 'Achilles_fitness_decrease', 'Achilles_fitness_increase', 'Aging_Perturbations_from_GEO_down', 'Aging_Perturbations_from_GEO_up', 'Allen_Brain_Atlas_down', 'Allen_Brain_Atlas_up', 'BioCarta_2013', 'BioCarta_2015', 'BioCarta_2016', 'BioPlex_2017', 'CORUM', 'Cancer_Cell_Line_Encyclopedia', 'ChEA_2013', 'ChEA_2015', 'ChEA_2016', 'Chromosome_Location', 'Chromosome_Location_hg19', 'DSigDB', 'Data_Acquisition_Method_Most_Popular_Genes', 'DepMap_WG_CRISPR_Screens_Broad_CellLines_2019', 'DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019', 'DisGeNET', 'Disease_Perturbations_from_GEO_down', 'Disease_Perturbations_from_GEO_up', 'Disease_Signatures_from_GEO_down_2014', 'Disease_Signatures_from_GEO_up_2014', 'DrugMatrix', 'Drug_Perturbations_from_GEO_2014', 'Drug_Perturbations_from_GEO_down', 'Drug_Perturbations_from_GEO_up', 'ENCODE_Histone_Modifications_2013', 'ENCODE_Histone_Modifications_2015', 'ENCODE_TF_ChIP-seq_2014', 'ENCODE_TF_ChIP-seq_2015', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'ESCAPE', 'Enrichr_Libraries_Most_Popular_Genes', 'Enrichr_Submissions_TF-Gene_Coocurrence', 'Epigenomics_Roadmap_HM_ChIP-seq', 'GO_Biological_Process_2013', 'GO_Biological_Process_2015', 'GO_Biological_Process_2017', 'GO_Biological_Process_2017b', 'GO_Biological_Process_2018', 'GO_Cellular_Component_2013', 'GO_Cellular_Component_2015', 'GO_Cellular_Component_2017', 'GO_Cellular_Component_2017b', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2013', 'GO_Molecular_Function_2015', 'GO_Molecular_Function_2017', 'GO_Molecular_Function_2017b', 'GO_Molecular_Function_2018', 'GTEx_Tissue_Sample_Gene_Expression_Profiles_down', 'GTEx_Tissue_Sample_Gene_Expression_Profiles_up', 'GWAS_Catalog_2019', 'GeneSigDB', 'Genes_Associated_with_NIH_Grants', 'Genome_Browser_PWMs', 'HMDB_Metabolites', 'HomoloGene', 'HumanCyc_2015', 'HumanCyc_2016', 'Human_Gene_Atlas', 'Human_Phenotype_Ontology', 'InterPro_Domains_2019', 'Jensen_COMPARTMENTS', 'Jensen_DISEASES', 'Jensen_TISSUES', 'KEA_2013', 'KEA_2015', 'KEGG_2013', 'KEGG_2015', 'KEGG_2016', 'KEGG_2019_Human', 'KEGG_2019_Mouse', 'Kinase_Perturbations_from_GEO_down', 'Kinase_Perturbations_from_GEO_up', 'LINCS_L1000_Chem_Pert_down', 'LINCS_L1000_Chem_Pert_up', 'LINCS_L1000_Kinase_Perturbations_down', 'LINCS_L1000_Kinase_Perturbations_up', 'LINCS_L1000_Ligand_Perturbations_down', 'LINCS_L1000_Ligand_Perturbations_up', 'Ligand_Perturbations_from_GEO_down', 'Ligand_Perturbations_from_GEO_up', 'MCF7_Perturbations_from_GEO_down', 'MCF7_Perturbations_from_GEO_up', 'MGI_Mammalian_Phenotype_2013', 'MGI_Mammalian_Phenotype_2017', 'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotype_Level_4', 'MGI_Mammalian_Phenotype_Level_4_2019', 'MSigDB_Computational', 'MSigDB_Oncogenic_Signatures', 'Microbe_Perturbations_from_GEO_down', 'Microbe_Perturbations_from_GEO_up', 'Mouse_Gene_Atlas', 'NCI-60_Cancer_Cell_Lines', 'NCI-Nature_2015', 'NCI-Nature_2016', 'NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions', 'NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions', 'NIH_Funded_PIs_2017_Human_AutoRIF', 'NIH_Funded_PIs_2017_Human_GeneRIF', 'NURSA_Human_Endogenous_Complexome', 'OMIM_Disease', 'OMIM_Expanded', 'Old_CMAP_down', 'Old_CMAP_up', 'PPI_Hub_Proteins', 'Panther_2015', 'Panther_2016', 'Pfam_Domains_2019', 'Pfam_InterPro_Domains', 'Phosphatase_Substrates_from_DEPOD', 'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO', 'Rare_Diseases_AutoRIF_ARCHS4_Predictions', 'Rare_Diseases_AutoRIF_Gene_Lists', 'Rare_Diseases_GeneRIF_ARCHS4_Predictions', 'Rare_Diseases_GeneRIF_Gene_Lists', 'Reactome_2013', 'Reactome_2015', 'Reactome_2016', 'SILAC_Phosphoproteomics', 'Single_Gene_Perturbations_from_GEO_down', 'Single_Gene_Perturbations_from_GEO_up', 'SubCell_BarCode', 'SysMyo_Muscle_Gene_Sets', 'TF-LOF_Expression_from_GEO', 'TF_Perturbations_Followed_by_Expression', 'TRANSFAC_and_JASPAR_PWMs', 'TRRUST_Transcription_Factors_2019', 'TargetScan_microRNA', 'TargetScan_microRNA_2017', 'Tissue_Protein_Expression_from_Human_Proteome_Map', 'Tissue_Protein_Expression_from_ProteomicsDB', 'Transcription_Factor_PPIs', 'UK_Biobank_GWAS_v1', 'VirusMINT', 'Virus_Perturbations_from_GEO_down', 'Virus_Perturbations_from_GEO_up', 'WikiPathways_2013', 'WikiPathways_2015', 'WikiPathways_2016', 'WikiPathways_2019_Human', 'WikiPathways_2019_Mouse', 'dbGaP', 'huMAP', 'miRTarBase_2017']

"""

# def get_target_genes(annotation_file, number):
#     """
#     annotation_file:
#     Gene    Total   Exon    Intron  Exon_samples    Intron_samples
#     PTPRN2  14      1       13      M594    M430,M452,M517,M530,M532,M534,M589,M594,M602,M609,M636,M669,M671        
#     MAD1L1  9       0       9       -       M426,M489,M507,M524,M546,M589,M605,M625,M671    
#     CDH4    8       0       8       -       M425,M426,M546,M579,M625,M633,M668,M671 
    
#     annotation_file2:
#     Gene    Enhancer        Exon    Intron  Promoter        Utr     Total
#     PTPRN2  M446,M667       M416    M416,M416,M425,M426,M426,M426,M426,M426,\
#     """
#     number = int(number)
#     anno_h = open(annotation_file, "r")
#     header = anno_h.readline().strip()

#     GeneSets = []
#     for line in anno_h:
#         lines = line.strip().split("\t")
#         gene = lines[0]
#         try:
#             sampleNum = int(lines[1])
#             if sampleNum >= number:
#                 GeneSets.append(gene)
#         except ValueError:
#             total = lines[-1]
#             totalLen = len(total)
#             if totalLen >= number:
#                 GeneSets.append(gene)
#     anno_h.close()

#     return GeneSets




def get_target_genes(gene_file):
    targetGene = set()

    gene_h = open(gene_file, "r")
    for line in gene_h:
        lines = line.strip().split("\t")
        if not lines[0].lower().startswith("gene"):
            gene = lines[0]
            targetGene.add(gene)
    gene_h.close()

    return list(targetGene)




def gene_enrichment(annotation_file, library, outdir): #number, 
    """
    :param gene_list: Flat file with list of genes, one gene id per row, or a python list object
    :param gene_sets: Enrichr Library to query. Required enrichr library name(s). Separate each name by comma.
    :param organism: Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm).
                     see here for details: https://amp.pharm.mssm.edu/modEnrichr
    :param description: name of analysis. optional.
    :param outdir: Output file directory
    :param float cutoff: Show enriched terms which Adjusted P-value < cutoff. 
                         Only affects the output figure. Default: 0.05
    :param int background: BioMart dataset name for retrieving background gene information.
                           This argument only works when gene_sets input is a gmt file or python dict.
                           You could also specify a number by yourself, e.g. total expressed genes number.
                           In this case, you will skip retrieving background infos from biomart.
    """
    # GeneSets = get_target_genes(annotation_file, number)

    targetGene = get_target_genes(annotation_file)

    if not os.path.exists(outdir):
        cmd = "mkdir -p %s" % outdir
        os.system(cmd)

    Librarys = library.split(",")
    Librarys = [ll.strip() for ll in Librarys]

    ### gene enrichment
    ### gseapy.enrichr(gene_list= GeneSets, description='pathway', gene_sets='KEGG_2019_Human', outdir= outdir)
    ### gseapy.enrichr(gene_list= GeneSets, gene_sets='KEGG_2019_Human', outdir= outdir)
    for lib in Librarys:
        gseapy.enrichr(gene_list= targetGene, organism="human", gene_sets=lib, outdir= outdir)

    return 0



def gene_enrichment_local(annotation_file, library, outdir):
    """
    :param gene_list: Flat file with list of genes, one gene id per row, or a python list object
    :param gene_sets: Enrichr Library to query. Required enrichr library name(s). Separate each name by comma.
    :param organism: Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm).
                     see here for details: https://amp.pharm.mssm.edu/modEnrichr
    :param description: name of analysis. optional.
    :param outdir: Output file directory
    :param float cutoff: Show enriched terms which Adjusted P-value < cutoff. 
                         Only affects the output figure. Default: 0.05
    :param int background: BioMart dataset name for retrieving background gene information.
                           This argument only works when gene_sets input is a gmt file or python dict.
                           You could also specify a number by yourself, e.g. total expressed genes number.
                           In this case, you will skip retrieving background infos from biomart.
    """
    # GeneSets = get_target_genes(annotation_file, number)

    targetGene = get_target_genes(annotation_file)

    if not os.path.exists(outdir):
        cmd = "mkdir -p %s" % outdir
        os.system(cmd)

    # Librarys = library.split(",")
    # Librarys = [ll.strip() for ll in Librarys]

    ### gene enrichment
    ### gseapy.enrichr(gene_list= GeneSets, description='pathway', gene_sets='KEGG_2019_Human', outdir= outdir)
    ### gseapy.enrichr(gene_list= GeneSets, gene_sets='KEGG_2019_Human', outdir= outdir)
    # geneSets = get_target_genes(library)

    gseapy.enrichr(gene_list= targetGene, gene_sets=library,  organism="human",  outdir= outdir)

    return 0






def main():
    parser = argparse.ArgumentParser(description="Gene enrichment for the annotation file.")
    parser.add_argument("-a", "--annotation", help="The input annotation file.")
    # parser.add_argument("-n", "--number", default=1, help="The number threshold of samples for given gene.")
    parser.add_argument("-l", "--library", default="GO_Biological_Process_2018,GO_Molecular_Function_2018,GO_Cellular_Component_2018,KEGG_2019_Human,GWAS_Catalog_2019,UK_Biobank_GWAS_v1,MGI_Mammalian_Phenotype_Level_4_2019,OMIM_Disease,OMIM_Expanded,Reactome_2016,GTEx_Tissue_Sample_Gene_Expression_Profiles_down,GTEx_Tissue_Sample_Gene_Expression_Profiles_up,Tissue_Protein_Expression_from_Human_Proteome_Map,Tissue_Protein_Expression_from_ProteomicsDB,WikiPathways_2019_Human,PPI_Hub_Proteins", help="The selected libraries.")
    parser.add_argument("-o", "--outdir", help="The output directory.")
    parser.add_argument("-m", "--method", help="The method used, 'website' or 'local'.")
    args = parser.parse_args()
    if args.method.lower() == "website":
        gene_enrichment(args.annotation, args.library, args.outdir) #args.number,
    elif args.method.lower() == "local":
        gene_enrichment_local(args.annotation, args.library, args.outdir)
    else:
        print("Please make sure the argument 'method' is 'website' or 'local'.")
        sys.exit(1)

if __name__ == "__main__":
    main()



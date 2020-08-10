#!/usr/bin/python
#-*- coding:utf-8 -*-
import os
import sys
import yaml
import collections
from tinyfasta import FastaParser

#usage: snakemake -p -s ~/github/NanoHub/pipeline/Population.pipeline.py --configfile ~/github/NanoHub/pipeline/Population.pipeline.yaml -j 24


### Software environment
CondaENV = config["CondaENV"]
CondaBIN = CondaENV + '/bin'

### Snakemake workflow environment
PIPENV = config["PIPENV"]
SCRIPT_DIR =  PIPENV + '/script'
PIPE_DIR = PIPENV + '/pipeline'
SRC_DIR = PIPENV + '/src/pgc'
RULE_DIR = PIPENV +  '/rule'

### General configuation
IN_PATH = config['ProjectPath']
THREADS = config["THREADS"]
ThreadFold = config["ThreadFold"]
SAMPLES = config["SAMPLES"]
CHRS = config["CHRS"]
AUTOSOMES = config["AUTOSOMES"]
SVTYPES = config["SVTYPES"]
CATEGORY = config["CATEGORY"]
TRAITS = config["TRAITS"]
PHENOS = config["PHENOS"]

### shuffle list for estimate of Fst  threshold
# SHUFFLENUM = list(range(1000))
SHUFFLENUM = list(range(100))

### randomly select of sample size, "405" indicates max number of samples.
NUMBERS = list(range(1, 405))



### Rules of snakemake  
#### Base rule
include: RULE_DIR + '/databaseFormat.rule.py'
include: RULE_DIR + '/BaseFunction.rule.py'
#### Function rule
# include: RULE_DIR + "/PopQuality.rule.py"
# include: RULE_DIR + "/NanoMappingStats.rule.py"
# include: RULE_DIR + '/NanoCallSV.rule.py'
include: RULE_DIR + '/PopMultipleSV.rule.py'
include: RULE_DIR + "/PopSVFilt.rule.py"
include: RULE_DIR + "/PopFrequency.rule.py"
include: RULE_DIR + "/populationGeno.rule.py"
include: RULE_DIR + "/PopStructure.rule.py"
include: RULE_DIR + "/PopStatistics.rule.py"
include: RULE_DIR + "/PopAnnotation.rule.py"
include: RULE_DIR + "/PopOverlap.rule.py"
include: RULE_DIR + "/PopMergeSeq.rule.py"
include: RULE_DIR + "/NanoAssembly.rule.py"



rule all:
    input:
        #################### Quality control #######################
        expand(IN_PATH + "/clean/{sample}.fastq.gz", sample=SAMPLES),
        expand(IN_PATH + "/QualityControl/{sample}_stats.xls", sample=SAMPLES),
        IN_PATH + "/QualityControl/Samples_quality_summary.xls",
        IN_PATH + "/QualityControl/Samples_quality_summary.xls",
        IN_PATH + "/QualityControl/Samples_read_length_hist.pdf",
        # ################## mapping statistics ##############
        IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.xls",
        IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.pdf",
        IN_PATH + "/mapping/Samples_mapping_summary.xls",
        IN_PATH + "/mapping/Samples_mapping_summary.pdf",
        ##################### Call SV ###############################
        expand(IN_PATH + "/SVCall/Sniffles/minimap2/{sample}.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/NanoSV/minimap2/{sample}.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/NanoSV/minimap2/{sample}_type.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/nanovar/{sample}/{sample}.nanovar.pass.vcf", sample=SAMPLES),
        ####################  Filt SV ################################
        expand(IN_PATH + "/SVCall/Final/{sample}_common_original.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/Final/{sample}_common_original_read3.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_length.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/Final/bed/{sample}_SV.bed", sample=SAMPLES),
        expand(IN_PATH + "/mapping/minimap2/bed/{sample}_SV_overlap.bed", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/FiltStats/{sample}_SV_chr_length.xls", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/FiltStats/{sample}_summary.xls", sample=SAMPLES),
        IN_PATH + "/SVCall/FiltStats/Samples_SV_type_number.xls",
        IN_PATH + "/SVCall/FiltStats/Samples_SV_number_hist.pdf",
        IN_PATH + "/SVCall/FiltStats/Samples_SV_length_summary.xls",
        IN_PATH + "/SVCall/FiltStats/Samples_types_summary.xls",
        IN_PATH + "/SVCall/FiltStats/Samples_types_filt_steps.pdf",
        #################### Merge SV ###########################
        IN_PATH + "/population/Merge/Sample_SV_common.vcf",
        IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
        IN_PATH + "/population/Merge/Sample_SV_common_CN329.vcf",
        IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
        IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL.txt",
        IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL_length_1000.pdf",
        # # # ####################### SV Statistics ################################
        IN_PATH + "/population/Stats/Sample_merge_summary.xls",
        IN_PATH + "/population/Stats/Sample_merge_record_length_density.pdf",
        IN_PATH + "/population/Stats/Sample_merge_dist_summary.xls",
        IN_PATH + "/population/Stats/Sample_merge_chroms_distribution.pdf",
        IN_PATH + "/population/Stats/Sample_merge_SV_number_ChrLength_cor.pdf",
        IN_PATH + "/population/Stats/Sample_merge_distribution.pdf",
        IN_PATH + "/population/Stats/Sample_merge_dist_summary.pdf",
        IN_PATH + "/population/Stats/Sample_common_SV_chr_length.xls",
        IN_PATH + "/population/Stats/Sample_common_SV_chr_length.pdf",
        IN_PATH + "/population/Stats/Sample_DEL.xls",
        IN_PATH + "/population/Stats/Sample_DEL_circos.xls",
        IN_PATH + "/population/Stats/Sample_merge_tra_circos.xls",
        IN_PATH + "/population/Stats/Sample_SV_type_tag_individulas.pdf",
        IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        IN_PATH + "/population/Stats/Sample_SV_distance_category.pdf",
        # ################################## Density ####################################
        IN_PATH + "/population/SVDensity/Sample_common_to_terminal_distance.xls",
        IN_PATH + "/population/SVDensity/Sample_common_to_terminal_distance.pdf",
        expand(IN_PATH + "/population/SVDensity/Sample_common_SV_{SVtype}.vcf", SVtype=SVTYPES),
        expand(IN_PATH + "/population/SVDensity/Sample_{SVtype}_to_terminal_distance.xls", SVtype=SVTYPES),
        expand(IN_PATH + "/population/SVDensity/Sample_{SVtype}_to_terminal_distance.pdf", SVtype=SVTYPES),
        ################################## SV genotype ###############################
        IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Singleton.txt",
        IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Singleton.txt",
        IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_category_summary.txt",
        IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_category_summary.pdf",
        IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
        IN_PATH + "/population/genotype/Sample_SV_genotype_frequency_stat.xls",
        IN_PATH + "/population/genotype/Sample_SV_genotype_frequency_stat.pdf",
        IN_PATH + "/population/plink/Sample_SV_geno.ped",
        IN_PATH + "/population/plink/Sample_SV_geno.bed",
        ############################ sub-group genotype  ######################
        IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        IN_PATH + "/population/subgroup/Sample_SV_geno.ped",
        IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
        # ############################## population structure #####################
        IN_PATH + "/population/admixture/log1.out",
        IN_PATH + "/population/admixture/Sample_assign_population.xls",
        IN_PATH + "/population/admixture/Sample_assign_population.pdf",
        IN_PATH + "/population/eigenstrat/Sample_SV.ind",
        IN_PATH + "/population/eigenstrat/Sample_SV.pca.txt",
        IN_PATH + "/population/eigenstrat/Sample_PCA.txt",
        IN_PATH + "/population/eigenstrat/Sample_PCA-province_region.txt",
        IN_PATH + "/population/eigenstrat/Sample_PCA12-province.pdf",
        # # ####################### Fst ###################################
        expand(IN_PATH + "/population/ShuffleFst/group/Sample_SV_genotype_fst_{number}.txt", number=SHUFFLENUM),
        expand(IN_PATH + "/population/ShuffleFst/Sample_SV_genotype_fst_{number}.txt", number=SHUFFLENUM),
        IN_PATH + "/population/Fst/Fst_permutation_values_sorted.txt",
        IN_PATH + "/population/Fst/Sample_SV_genotype_fst.txt",
        IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig.txt",
        IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig_overlap.bed",
        IN_PATH + "/population/Fst/Sample_SV_genotype_fst_table.txt",
        IN_PATH + "/population/Fst/Sample_SV_genotype_fst.pdf",
        IN_PATH + "/population/Fst/Sample_SV_genotype_fst_manhattan.pdf",
        IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst.txt",
        IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst.pdf",
        IN_PATH + "/population/ShuffleFstAge/Fst_permutation_values_sorted.txt",
        IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_30_60.txt",
        IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_sig_overlap.bed",
        # # # ############################### Annotation ###############################
        IN_PATH + "/population/Annotation/Sample_common_SV.tsv",
        IN_PATH + "/population/Annotation/Sample_common_SV_modify.tsv",
        IN_PATH + "/population/bed/Sample_common_SV_DEL.bed",
        IN_PATH + "/population/bed/Sample_common_SV_INS_point.bed",
        IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP_point.bed",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
        IN_PATH + "/population/Annotation/Sample_common_SV_annotation_stats.xls",
        IN_PATH + "/population/Annotation/Sample_common_SV_annotation_promoter_stats.xls",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof.bed",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_stats.txt",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover_inv_cds_stats.txt",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion.txt",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion_IGV.batch",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding.bed",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq.bed",
        IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Expand_gene_tag_genotypes.txt",
        IN_PATH + "/population/Annotation/Enrichment/Promoter/COSMIC_Disease_gene_tag_genotypes.txt",
        IN_PATH + "/population/Annotation/Enrichment/Promoter/Annotation_database_category.pdf",
        IN_PATH + "/population/Annotation/Sample_common_SV_anno_odds_ratio.xls",
        IN_PATH + "/population/Annotation/Sample_common_SV_anno_odds_ratio.pdf",
        IN_PATH + "/population/Annotation/Sample_common_SV_overlap_circos.txt",
        IN_PATH + "/population/Annotation/Enrichment/All/KEGG_2019_Human..enrichr.reports.txt",
        IN_PATH + "/population/Annotation/Enrichment/All/Plots/KEGG_2019_Human..enrichr.reports.pdf",
        IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Disease_gene_tag_genotypes.txt",
        IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_Disease_gene_tag_genotypes.txt",
        IN_PATH + "/population/Annotation/Enrichment/All/Annotation_database_category.txt",
        IN_PATH + "/population/Annotation/Enrichment/All/Annotation_database_category.pdf",
        IN_PATH + "/population/Annotation/Enrichment/All/disease_no_catalog_gene.txt",
        IN_PATH + "/population/Annotation/Enrichment/All/isease_no_catalog_gene_malacards.txt",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist_category.txt",
        IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist_category.pdf",
        IN_PATH + "/population/bed/frequency/Sample_common_SV_Rare.bed",
        IN_PATH + "/population/bed/frequency/Sample_common_SV_rare_genepred_overlap.bed",
        # # # # # # ############################# Overlap #############################
        IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap_filt.bed",
        IN_PATH + "/population/bed/LRS15/SV_overlap_filt_chr_length_summary.xls",
        IN_PATH + "/population/bed/LRS15/Sample_common_SV_number.pdf",
        expand(IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap.bed", SVtype=SVTYPES),
        expand(IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        IN_PATH + "/population/bed/gnomAD/Sample_overlap_filt_chr_length_summary.xls",
        IN_PATH + "/population/bed/gnomAD/Sample_common_overlap_SV_number.pdf",
        expand(IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap.bed", SVtype=SVTYPES),
        IN_PATH + "/population/bed/DGV/Sample_overlap_filt_chr_length_summary.xls",
        IN_PATH + "/population/bed/DGV/Sample_common_overlap_SV_length.pdf",
        expand(IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap.bed", SVtype=SVTYPES),
        expand(IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tags.xls",
        IN_PATH + "/population/bed/WGS17795/Sample_common_overlap_SV_number.pdf",
        expand(IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap.bed", SVtype=SVTYPES),
        expand(IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        IN_PATH + "/population/bed/WGS911/Sample_overlap_filt_chr_length_summary.xls",
        IN_PATH + "/population/bed/WGS911/Sample_common_overlap_SV_number.pdf",
        expand(IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt_cds.bed", SVtype=SVTYPES),
        IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap.xls",
        IN_PATH + "/population/bed/All/Sample_overlap_summary.txt",
        IN_PATH + "/population/bed/All/Sample_overlap_summary.pdf",
        IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tag_freq.xls",
        IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tag_freq_stat.xls",
        IN_PATH + "/population/bed/SV_overlap_filt_tag_freq_stat_summary.xls",
        IN_PATH + "/population/bed/SV_overlap_filt_tag_freq_stat_summary.pdf",
        IN_PATH + "/population/bed/database/SV_overlap_filt_tags.xls",
        IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq.xls",
        IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq_stat.xls",
        IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq_stat.pdf",
        IN_PATH + "/population/bed/novel/SV_type_known_and_novwl_stats.txt",
        IN_PATH + "/population/bed/novel/SV_type_known_and_novwl_stats.pdf",
        # ########################### Category #########################
        IN_PATH + "/population/Category/Sample_type_freq_percentage.pdf",
        IN_PATH + "/population/Category/Sample_SV_poly_hist.pdf",
        IN_PATH + "/population/Category/Sample_SV_types_hist.pdf",
        IN_PATH + "/population/Category/Kown_and_novel_tags_stats.txt",
        IN_PATH + "/population/Category/Kown_and_novel_tags_stats.pdf",
        IN_PATH + "/population/bed/database/Sample_SV_overlap_category_stats.pdf",
        IN_PATH + "/population/bed/LRS15/Sample_SV_overlap_category_stats.xls",
        IN_PATH + "/population/bed/LRS15/Sample_SV_overlap_category_stats.pdf",
        IN_PATH + "/population/bed/gnomAD/Sample_SV_overlap_category_stats.pdf",
        IN_PATH + "/population/bed/DGV/Sample_SV_overlap_category_stats.pdf",
        IN_PATH + "/population/bed/WGS911/Sample_SV_overlap_category_stats.pdf",
        ######################## novel ####################
        expand(IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag_seq.fasta", SVtype=SVTYPES),
        IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common.xls",
        IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_gene.xls",
        # # ############################ Heter to homo ####################################
        IN_PATH + "/population/genotype/Sample_SV_heter2homo.txt",
        IN_PATH + "/population/genotype/Sample_SV_heter2homo.pdf",
        IN_PATH + "/population/genotype/Sample_SV_type_heter2homo.txt",
        IN_PATH + "/population/genotype/Sample_SV_type_heter2homo.pdf",
        IN_PATH + "/population/genotype/homo/Sample_SV_genotype_homo.txt",
        IN_PATH + "/population/genotype/homo/Sample_SV_tag_homo_category.txt",
        IN_PATH + "/population/genotype/homo/Sample_homo_SV_genepred_overlap.bed",
        IN_PATH + "/population/genotype/homo/Sample_homo_SV_annotation_stats.xls",
        expand(IN_PATH + "/population/genotype/shuffle/Sample_geno_{number}_stat.txt", number=NUMBERS),
        IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary.txt",
        IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary_reshape.txt",
        IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary_reshape.pdf",
        # ############################ sequence ############################
        # # IN_PATH + "/population/Merge/Sample_common_SV_tag_seq.fasta",
        # # IN_PATH + "/population/Repeat/Sample_common_SV_tag_DEL.fasta",
        # # IN_PATH + "/population/Repeat/Sample_common_SV_tag_DUP.fasta",
        # # expand(IN_PATH + "/population/Repeat/{SVtype}/Sample_common_SV_tag_{SVtype}.fasta.out", SVtype=SVTYPES),
        # # IN_PATH + "/population/Dfam/SV_seq_tblout.txt",
        # # IN_PATH + "/population/Dfam/SV_seq_dfam_type.xls",
        # # IN_PATH + "/population/Dfam/SV_seq_dfam_INS_pie.pdf",
        # ################### IGV #############################
        # expand(IN_PATH + "/IGV/{sample}_SV.bed", sample=SAMPLES),
        # expand(IN_PATH + "/IGV/{sample}_SV.batch", sample=SAMPLES),
        # expand(IN_PATH + "/IGV/{sample}_SV_run.sh", sample=SAMPLES),
        ################## Assembly ############################
        expand(IN_PATH + "/Assembly/{sample}/{sample}_assembly.ctg.lay.gz", sample=SAMPLES),
        expand(IN_PATH + "/Assembly/{sample}/{sample}_assembly.fasta", sample=SAMPLES),
        expand(IN_PATH + "/Assembly/mapping/{sample}_assembly.bam", sample=SAMPLES),
        expand(IN_PATH + "/Assembly/{sample}/{sample}_assembly_polish.fasta", sample=SAMPLES),
        expand(IN_PATH + '/Assembly/Quast/{sample}/polish/report.txt', sample=SAMPLES),
        IN_PATH + "/Assembly/Assembly_stats.pdf",






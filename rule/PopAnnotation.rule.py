
############################## annotation based on genepred ###############
rule GenePredOverlap:
    input:
        # bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
        bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP_point.bed",

    output:
        genepred = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
    params:
        genePred = config["genePred"],
    log:
        IN_PATH + "/log/GenePredOverlap.log",
    run:
        shell("bedtools intersect -a {input.bed} -b {params.genePred} -wa -wb > {output.genepred} 2>{log}")



#################################################################################################



###################################################################################
####################################################################################
###################################################################################
################################### tag ferquency #################################
rule TagCalss:
    input:
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
        bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
    output:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        sig = IN_PATH + "/population/bed/frequency/Sample_common_SV_Singleton.bed",
        rare = IN_PATH + "/population/bed/frequency/Sample_common_SV_Rare.bed",
        low = IN_PATH + "/population/bed/frequency/Sample_common_SV_Low.bed",
        common = IN_PATH + "/population/bed/frequency/Sample_common_SV_Common.bed",
    params:
        TagFreqClass = SRC_DIR + "/TagFreqClass.py",
        outPrefix = IN_PATH + "/population/bed/frequency/Sample_common_SV",
    log:
        IN_PATH + "/log/TagCalss.log",
    run:
        shell("python {params.TagFreqClass} --frequency {input.freq} --bed {input.bed} --tag {output.category} --outPrefix {params.outPrefix} > {log} 2>&1")





rule FreqOverlap:
    input:
        sig = IN_PATH + "/population/bed/frequency/Sample_common_SV_Singleton.bed",
        rare = IN_PATH + "/population/bed/frequency/Sample_common_SV_Rare.bed",
        low = IN_PATH + "/population/bed/frequency/Sample_common_SV_Low.bed",
        common = IN_PATH + "/population/bed/frequency/Sample_common_SV_Common.bed",
    output:
        sig = IN_PATH + "/population/bed/frequency/Sample_common_SV_singleton_genepred_overlap.bed",
        rare = IN_PATH + "/population/bed/frequency/Sample_common_SV_rare_genepred_overlap.bed",
        low = IN_PATH + "/population/bed/frequency/Sample_common_SV_low_genepred_overlap.bed",
        common = IN_PATH + "/population/bed/frequency/Sample_common_SV_common_genepred_overlap.bed",
    params:
        genePred = config["genePred"],
    log:
        IN_PATH + "/log/FreqOverlap.log",
    run:
        shell("bedtools intersect -a {input.sig} -b {params.genePred} -wa -wb > {output.sig} 2>{log}")
        shell("bedtools intersect -a {input.rare} -b {params.genePred} -wa -wb > {output.rare} 2>>{log}")
        shell("bedtools intersect -a {input.low} -b {params.genePred} -wa -wb > {output.low} 2>>{log}")
        shell("bedtools intersect -a {input.common} -b {params.genePred} -wa -wb > {output.common} 2>>{log}")





rule AnnoStats:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
        # category = IN_PATH + "/population/bed/Sample_common_SV_frequency_tag.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        stat = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_stats.xls",
        ratio = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_stats_ratio.xls",
    params:
        CategoryAnnoStats = SRC_DIR + "/CategoryAnnoStats.py",
    log:
        IN_PATH + "/log/AnnoStats.log",
    run:
        shell("python {params.CategoryAnnoStats} --category {input.category} --annotation {input.anno} --out {output.ratio} --stats {output.stat} > {log} 2>&1")




################################### Stats add enhancer #########################

rule AnnoEnhancerStats:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
        # regulator = IN_PATH + "/population/Annotation/Sample_common_SV_regulator_all.bed",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        # anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all.bed",
        select = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_promoter_stats_select.xls",
        stat = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_promoter_stats.xls",
        ratio = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_promoter_stats_ratio.xls",
    params:
        # CategoryAnnoStats = SRC_DIR + "/CategoryAnnoStats.py",
        hierarchyAnnotation = SRC_DIR + "/hierarchyAnnotation.py",
    log:
        IN_PATH + "/log/AnnoEnhancerStats.log",
    run:
        # shell("cat {input.anno}  {input.regulator} > {output.anno}")
        # shell("python {params.CategoryAnnoStats} --category {input.category} --annotation {output.anno} --out {output.ratio} --stats {output.stat} > {log} 2>&1")
        shell("python {params.hierarchyAnnotation} --annotation {input.anno} --category {input.category} --out {output.select} --stats {output.stat} --ratio {output.ratio} > {log} 2>&1")





rule overlapFilt:
    input:
        overlap = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
    output:
        out = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof.bed",
        filt = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover.bed",
    params:
        geneOverlapFilt = SRC_DIR + "/geneOverlapFilt.py",
        ensGene = config["ensGene"],
    log:
        IN_PATH + "/log/overlapFilt.log",
    run:
        shell("python {params.geneOverlapFilt} --gene {params.ensGene} --overlap {input.overlap} --out {output.out} --filt {output.filt} > {log} 2>&1")


rule lofCDS:
    input:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof.bed",
    output:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
    run:
        shell("grep CDS {input.lof} > {output.lof}")












rule LofStats:
    input:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
        matrix = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    output:
        stat = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_stats.txt",
        sample = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion.txt",
    params:
        LoFStats = SRC_DIR + "/LoFStats.py",
    log:
        IN_PATH + "/log/LofStats.log",
    run:
        shell("python {params.LoFStats} --tag {input.lof} --matrix {input.matrix} --out {output.stat} --name LoF --sample {output.sample} > {log} 2>&1")



rule INVDUPStat:
    input:
        filt = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover.bed",
        matrix = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    output:
        inv = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover_inv_cds.bed",
        dup = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover_dup_cds.bed",
        invStat = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover_inv_cds_stats.txt",
        dupStat = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_cover_dup_cds_stats.txt",
        temp = IN_PATH + "/population/Annotation/all_cover_dup_inv_cds_temp.txt",
    params:
        LoFStats = SRC_DIR + "/LoFStats.py",
    log:
        IN_PATH + "/log/INVDUPStat.log",
    run:
        shell("grep INV {input.filt} | grep CDS > {output.inv}")
        shell("grep DUP {input.filt} | grep CDS > {output.dup}")
        shell("python {params.LoFStats} --tag {output.inv} --matrix {input.matrix} --out {output.invStat} --name INV --sample {output.temp} > {log} 2>&1")
        shell("python {params.LoFStats} --tag {output.dup} --matrix {input.matrix} --out {output.dupStat} --name DUP --sample {output.temp} >>{log}")



rule longDELIGV:
    input:
        tag = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion.txt",
    output:
        IGVbatch  = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_long_deletion_IGV.batch",
    params:
        longDeletionIGV = SRC_DIR + "/longDeletionIGV.py",
        IGVDir = "/home/wuzhikun/Project/Population/IGVRegion/LongDEL",
    log:
        IN_PATH + "/log/longDELIGV.log",
    run:
        shell("python {params.longDeletionIGV} --tag {input.tag} --out {output.IGVbatch} --IGVDir {params.IGVDir} > {log} 2>&1")



################################################################################




##################################################################################################
rule OddsRatio:
    input:
        # stat = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_stats.xls",
        stat = IN_PATH + "/population/Annotation/Sample_common_SV_annotation_promoter_stats.xls",
    output:
        odds = IN_PATH + "/population/Annotation/Sample_common_SV_anno_odds_ratio.xls",
        pvalue = IN_PATH + "/population/Annotation/Sample_common_SV_anno_odds_pvalue.xls",
    params:
        OddsRatio = SCRIPT_DIR + "/OddsRatio.R",
    log:
        IN_PATH + "/log/OddsRatio.log",
    run:
        shell("Rscript {params.OddsRatio} --input {input.stat} --odds {output.odds} --pvalue {output.pvalue} > {log} 2>&1")




rule OddsRatioPlot:
    input:
        odds = IN_PATH + "/population/Annotation/Sample_common_SV_anno_odds_ratio.xls",
    output:
        pdf = IN_PATH + "/population/Annotation/Sample_common_SV_anno_odds_ratio.pdf",
    params:
        OddsRatioPlot = SCRIPT_DIR + "/OddsRatioPlot.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/OddsRatioPlot.log",
    run:
        shell("Rscript {params.OddsRatioPlot} --input {input.odds} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")




rule DELOverlap:
    input:
        # bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS.bed",
        genepred = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
    output:
        genepred = IN_PATH + "/population/Annotation/Sample_DEL_INS_genepred_overlap.bed",
    # params:
    #     genePred = config["genePred"],
    log:
        IN_PATH + "/log/DELOverlap.log",
    run:
        # shell("bedtools intersect -a {input.bed} -b {params.genePred} -wa -wb > {output.genepred} 2>{log}")
        shell("grep -v 'INV' {input.genepred} | grep -v 'DUP' > {output.genepred}")

rule DELStats:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_DEL_INS_genepred_overlap.bed",
        # category = IN_PATH + "/population/bed/Sample_common_SV_frequency_tag.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        stat = IN_PATH + "/population/Annotation/Sample_DEL_INS_annotation_stats.xls",
        ratio = IN_PATH + "/population/Annotation/Sample_DEL_INS_annotation_stats_ratio.xls",
    params:
        CategoryAnnoStats = SRC_DIR + "/CategoryAnnoStats.py",
    log:
        IN_PATH + "/log/AnnoStats.log",
    run:
        shell("python {params.CategoryAnnoStats} --category {input.category} --annotation {input.anno} --out {output.ratio} --stat {output.stat} > {log} 2>&1")





rule overlapCircos:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
    output:
        region = IN_PATH + "/population/Annotation/Sample_common_SV_overlap_region.txt",
        circos = IN_PATH + "/population/Annotation/Sample_common_SV_overlap_circos.txt",
    params:
        SVTypeGeneFeature = SRC_DIR + "/SVTypeGeneFeature.py",
    log:
        IN_PATH + "/log/overlapCircos.log",
    run:
        shell("python {params.SVTypeGeneFeature} --annotation {input.anno} --out {output.region} --circos {output.circos} > {log} 2>&1")





###################################################################################








###################################### category enrichment #####################################




# rule CategoryGene:
#     input:
#         sig = IN_PATH + "/population/bed/frequency/Sample_common_SV_singleton_genepred_overlap.bed",
#         rare = IN_PATH + "/population/bed/frequency/Sample_common_SV_rare_genepred_overlap.bed",
#         low = IN_PATH + "/population/bed/frequency/Sample_common_SV_low_genepred_overlap.bed",
#         common = IN_PATH + "/population/bed/frequency/Sample_common_SV_common_genepred_overlap.bed",
#     output:
#         rare = IN_PATH + "/population/bed/frequency/Sample_common_SV_overlap_rare.bed",
#         rareGene = IN_PATH + "/population/bed/frequency/Sample_common_SV_overlap_rare_gene.txt",
#         commonGene = IN_PATH + "/population/bed/frequency/Sample_common_SV_overlap_common_gene.txt",
#     run:
#         shell("cat {input.sig} {input.rare} > {output.rare}")
#         shell("grep CDS {output.rare} | cut -f 8 | sort | uniq > {output.rareGene}")
#         shell("grep CDS {input.common} | cut -f 8 | sort | uniq > {output.commonGene}")






# rule categoryRareEnrich:
#     input:
#         rareGene = IN_PATH + "/population/bed/frequency/Sample_common_SV_overlap_rare_gene.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/frequency/Enrichment/Rare/KEGG_2019_Human..enrichr.reports.txt", 
#     params:
#         outdir = IN_PATH + "/population/bed/frequency/Enrichment/Rare/",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/DELEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("python {params.geneGSEA} --annotation {input.rareGene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")




# rule categoryRareEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/bed/frequency/Enrichment/Rare/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/frequency/Enrichment/Rare/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/bed/frequency/Enrichment/Rare/",
#         outdir = IN_PATH + "/population/bed/frequency/Enrichment/Rare/Plots/",
#         width = 8,
#         height = 3,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)






# rule categoryCommonEnrich:
#     input:
#         commonGene = IN_PATH + "/population/bed/frequency/Sample_common_SV_overlap_common_gene.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/frequency/Enrichment/Common/KEGG_2019_Human..enrichr.reports.txt", 
#     params:
#         outdir = IN_PATH + "/population/bed/frequency/Enrichment/Common/",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/DELEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("python {params.geneGSEA} --annotation {input.commonGene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")




# rule categoryCommonEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/bed/frequency/Enrichment/Common/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/frequency/Enrichment/Common/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/bed/frequency/Enrichment/Common/",
#         outdir = IN_PATH + "/population/bed/frequency/Enrichment/Common/Plots/",
#         width = 8,
#         height = 3,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)

###################################################################################################




########################################## ALL,  INS and DEL enrichment ###############################
rule AllCDSEnrich:
    input:
        bed = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
    output:
        gene = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_gene.txt",
        kegg = IN_PATH + "/population/Annotation/Enrichment/All/KEGG_2019_Human..enrichr.reports.txt", 
        omim = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Disease..enrichr.reports.txt",
        omim_expand = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Expanded..enrichr.reports.txt",
        gwas = IN_PATH + "/population/Annotation/Enrichment/All/GWAS_Catalog_2019..enrichr.reports.txt",
        clinvar = IN_PATH + "/population/Annotation/Enrichment/All/ClinVar_2019..enrichr.reports.txt",
        biobank = IN_PATH + "/population/Annotation/Enrichment/All/UK_Biobank_GWAS_v1..enrichr.reports.txt",
    params:
        outdir = IN_PATH + "/population/Annotation/Enrichment/All/",
        geneGSEA = SRC_DIR + "/geneGSEA.py",
        EnrichLibrary = config["EnrichLibrary"],
        number = 1,
    log:
        IN_PATH + "/log/DELEnrich.log",
    run:
        ### mu01
        ### --number {params.number} 
        # shell("grep CDS {input.bed} | cut -f 8 | sort | uniq > {output.gene}")
        shell("cut -f 8  {input.bed} | sort | uniq > {output.gene}")
        shell("python {params.geneGSEA} --annotation {output.gene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


rule AllCDSEnrichPlot:
    input:
        kegg = IN_PATH + "/population/Annotation/Enrichment/All/KEGG_2019_Human..enrichr.reports.txt",
    output:
        kegg = IN_PATH + "/population/Annotation/Enrichment/All/Plots/KEGG_2019_Human..enrichr.reports.pdf",
    params:
        geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
        indir = IN_PATH + "/population/Annotation/Enrichment/All/",
        outdir = IN_PATH + "/population/Annotation/Enrichment/All/Plots/",
        width = 8,
        height = 3,
        selectNum = 12,
    log:
        IN_PATH + "/log/CommonEnrichPlot.log",
    run:
        # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
        # print(cmd)
        # os.system(cmd)
        files = os.listdir(params.indir)
        for f in files:
            if f.endswith("enrichr.reports.txt"):
                fname = params.indir + f
                pdf = params.outdir + f.rstrip("txt") + "pdf"
                cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
                os.system(cmd)
########################################################################################################



#########################################################################################################

# rule DELEnrich:
#     input:
#         DEL = IN_PATH + "/population/Annotation/Sample_DEL_INS_genepred_overlap.bed",
#     output:
#         DEL = IN_PATH + "/population/Annotation/Sample_DEL_INS_genepred_overlap_gene.xls",
#         kegg = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/KEGG_2019_Human..enrichr.reports.txt", 
#         omim = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/OMIM_Disease..enrichr.reports.txt",
#         omim_expand = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/OMIM_Expanded..enrichr.reports.txt",
#         gwas = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/GWAS_Catalog_2019..enrichr.reports.txt",
#         clinvar = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/ClinVar_2019..enrichr.reports.txt",
#         biobank = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/UK_Biobank_GWAS_v1..enrichr.reports.txt",
#     params:
#         outdir = IN_PATH + "/population/Annotation/Enrichment/DEL_INS",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/DELEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("grep CDS {input.DEL} | cut -f 8 | sort | uniq > {output.DEL}")
#         shell("python {params.geneGSEA} --annotation {output.DEL} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")




# rule DELEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/",
#         outdir = IN_PATH + "/population/Annotation/Enrichment/DEL_INS/Plots/",
#         width = 8,
#         height = 3,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)




rule enrichGeneGeno:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
        # anno = IN_PATH + "/population/Annotation/Sample_DEL_INS_genepred_overlap.bed",
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        omim = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Disease..enrichr.reports.txt", 
        omim_expand = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Expanded..enrichr.reports.txt",
        clinvar = IN_PATH + "/population/Annotation/Enrichment/All/ClinVar_2019..enrichr.reports.txt",
        gwas = IN_PATH + "/population/Annotation/Enrichment/All/GWAS_Catalog_2019..enrichr.reports.txt",
        biobank = IN_PATH + "/population/Annotation/Enrichment/All/UK_Biobank_GWAS_v1..enrichr.reports.txt",
    output:
        gene = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Disease_gene_tag_genotypes.txt", 
        gene_expand = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Expand_gene_tag_genotypes.txt", 
        clinvar = IN_PATH + "/population/Annotation/Enrichment/All/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/All/GWAS_gene_tag_genotypes.txt", 
        biobank = IN_PATH + "/population/Annotation/Enrichment/All/Biobank_gene_tag_genotypes.txt", 
    log:
        IN_PATH + "/log/enrichGeneGeno.log",
    params:
        DiseaseGeneTagGeno = SRC_DIR + "/DiseaseGeneTagGeno.py",
    run:
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.omim} --annotation {input.anno} --genotype {input.geno} --out {output.gene} > {log} 2>&1")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.omim_expand} --annotation {input.anno} --genotype {input.geno} --out {output.gene_expand} 2>>{log}")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.clinvar} --annotation {input.anno} --genotype {input.geno} --out {output.clinvar} 2>>{log}")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.gwas} --annotation {input.anno} --genotype {input.geno} --out {output.gwas} 2>>{log}")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.biobank} --annotation {input.anno} --genotype {input.geno} --out {output.biobank} 2>>{log}")
        

rule cosmicOverlap:
    input:
        cosmic = "/home/wuzhikun/database/COSMIC/cancer_gene_census.gene",
        bed = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
    output:
        gene = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_gene.txt",
        cosmic = IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_gene_overlap_report.txt", 
    run:
        shell("cut -f 8 {input.bed} | sort | uniq > {output.gene}")
        shell("""echo "COSMIC\tCancer\tGene" > {output.cosmic}""")
        cmd = """cat %s %s | sort | uniq -c | sort -k 1nr | sed 's/      //g' | grep '^2' | cut -f 2 -d ' ' | awk '{print "COSMIC\tCancer\t"$0}' >> %s""" % (output.gene, input.cosmic, output.cosmic)
        os.system(cmd)





rule lofSVDist:
    input:
        bed = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
    output:
        dist = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist.txt",
    run:
        target_SV_distribution(input.freq, input.bed, output.dist)




rule CategoryLength:
    input:
        dist = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        dist = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist_category.txt",
    run:
        catagory_SV_distance(input.category, input.dist, output.dist)



rule CategoryLengthPlot:
    input:
        dist = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist_category.txt",
    output:
        pdf = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist_category.pdf",
    params:
        lengthBox = SCRIPT_DIR + "/lengthBox.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/CategoryLengthPlot.log",
    run:
        shell("Rscript {params.lengthBox} --input {input.dist} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")








rule cosmicGeno:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        cosmic = IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_gene_overlap_report.txt", 
    output:
        cosmic = IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_Disease_gene_tag_genotypes.txt",
    log:
        IN_PATH + "/log/cosmicGeno.log",
    params:
        DiseaseGeneTagGeno = SRC_DIR + "/DiseaseGeneTagGeno.py",
    run:
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.cosmic} --annotation {input.anno} --genotype {input.geno} --out {output.cosmic} > {log} 2>&1")




rule diseaseGene:
    input:
        gene_expand = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Expand_gene_tag_genotypes.txt", 
        # clinvar = IN_PATH + "/population/Annotation/Enrichment/All/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/All/GWAS_gene_tag_genotypes.txt", 
        # biobank = IN_PATH + "/population/Annotation/Enrichment/All/Biobank_gene_tag_genotypes.txt", 
        cosmic = IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_Disease_gene_tag_genotypes.txt",
    output:
        gene = IN_PATH + "/population/Annotation/Enrichment/All/disease_gene_circos.txt",
    params:
        geneDisease = SRC_DIR + "/geneDisease.py",
    log:
        IN_PATH + "/log/diseaseGene.log",
    run:
        Databases = ",".join([input.gene_expand, input.gwas, input.cosmic])
        Names = "OMIM,GWAS,COSMIC" 
        shell("python {params.geneDisease} --database {Databases} --name {Names} --out {output.gene} --feature CDS > {log} 2>&1")
    



rule diseaseCategoory:
    input:
        ommin = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Expand_gene_tag_genotypes.txt", 
        # clinvar = IN_PATH + "/population/Annotation/Enrichment/All/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/All/GWAS_gene_tag_genotypes.txt", 
        # biobank = IN_PATH + "/population/Annotation/Enrichment/All/Biobank_gene_tag_genotypes.txt", 
        cosmic = IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_Disease_gene_tag_genotypes.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        cat = IN_PATH + "/population/Annotation/Enrichment/All/Annotation_database_category.txt", 
    params:
        diseaseGeneCatrgory = SRC_DIR + "/diseaseGeneCatrgory.py",
        number = 405,
        name = "OMIM,COSMIC,GWAS",
    log:
        IN_PATH + "/log/diseaseCategoory.log",
    run:
        Database = ",".join([input.ommin, input.cosmic, input.gwas])
        shell("python {params.diseaseGeneCatrgory} --annotation {Database} --name {params.name} --category {input.category} --number {params.number} --out {output.cat} > {log} 2>&1")



rule diseaseCategooryPlot:
    input:
        cat = IN_PATH + "/population/Annotation/Enrichment/All/Annotation_database_category.txt", 
    output:
        pdf = IN_PATH + "/population/Annotation/Enrichment/All/Annotation_database_category.pdf", 
    params:
        diseaseCategory = SCRIPT_DIR + "/diseaseCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/diseaseCategooryPlot.log",
    run:
        shell("Rscript {params.diseaseCategory} --input {input.cat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")




rule nocatalogGene:
    input:
        all_gene = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_gene.txt",
        omim = IN_PATH + "/population/Annotation/Enrichment/All/OMIM_Expand_gene_tag_genotypes.txt", 
        # clinvar = IN_PATH + "/population/Annotation/Enrichment/All/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/All/GWAS_gene_tag_genotypes.txt", 
        # biobank = IN_PATH + "/population/Annotation/Enrichment/All/Biobank_gene_tag_genotypes.txt", 
        cosmic = IN_PATH + "/population/Annotation/Enrichment/All/COSMIC_Disease_gene_tag_genotypes.txt",
    output:
        combine = temp(IN_PATH + "/population/Annotation/Enrichment/All/disease_combine_gene_temp.txt"),
        nocat = IN_PATH + "/population/Annotation/Enrichment/All/disease_no_catalog_gene.txt",
    run:
        shell("cat {input.omim} {input.cosmic} {input.gwas} | cut -f 1 > {output.combine}")
        shell("cat {input.all_gene} {output.combine} | sort | uniq -c | sort -k 1n | cut -f 7- -d ' ' | grep '^1' | cut -f 2 -d ' ' > {output.nocat}")


rule GeneDisease:
    input:
        gene = IN_PATH + "/population/Annotation/Enrichment/All/disease_no_catalog_gene.txt",
    output:
        disease = IN_PATH + "/population/Annotation/Enrichment/All/isease_no_catalog_gene_malacards.txt",
    params:
        geneCardDisease = config["geneCardDisease"],
    run:
        cmd = """cat %s  | while read marker ; do sed -n '/'"$marker"'\t.*/p' %s ; done > %s """ % (input.gene, params.geneCardDisease, output.disease)
        os.system(cmd)


#####################################################################################################################################################





################################################# regulation region #######################################

rule lofnocoding:
    input:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof.bed",
    output:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding.bed",
    run:
        shell("grep Promoter {input.lof} > {output.lof}")
        shell("grep UTR {input.lof} >>{output.lof}")


def exclude_tag(target_tag, all_tag, out_file):
    TagSets = set()
    ta_h = open(target_tag, "r")
    for line in ta_h:
        lines = line.strip().split("\t")
        tag = lines[3]
        TagSets.add(tag)
    ta_h.close()

    all_h = open(all_tag, "r")
    out_h = open(out_file, "w")
    for line in all_h:
        line = line.strip()
        lines = line.split("\t")
        tag = lines[3]
        if tag not in TagSets:
            out_h.write("%s\n" % (line))
    all_h.close()
    out_h.close()



rule lofnocodingUniq:
    input:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding.bed",
        target = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
    output:
        lof = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq.bed",
    run:
        exclude_tag(input.target, input.lof, output.lof)





rule PromoterEnrich:
    input:
        bed = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq.bed",
    output:
        gene = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq_gene.txt",
        kegg = IN_PATH + "/population/Annotation/Enrichment/Promoter/KEGG_2019_Human..enrichr.reports.txt", 
        omim = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Disease..enrichr.reports.txt",
        omim_expand = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Expanded..enrichr.reports.txt",
        gwas = IN_PATH + "/population/Annotation/Enrichment/Promoter/GWAS_Catalog_2019..enrichr.reports.txt",
        clinvar = IN_PATH + "/population/Annotation/Enrichment/Promoter/ClinVar_2019..enrichr.reports.txt",
        biobank = IN_PATH + "/population/Annotation/Enrichment/Promoter/UK_Biobank_GWAS_v1..enrichr.reports.txt",
    params:
        outdir = IN_PATH + "/population/Annotation/Enrichment/Promoter/",
        geneGSEA = SRC_DIR + "/geneGSEA.py",
        EnrichLibrary = config["EnrichLibrary"],
        number = 1,
    log:
        IN_PATH + "/log/DELEnrich.log",
    run:
        ### mu01
        ### --number {params.number} 
        # shell("grep CDS {input.bed} | cut -f 8 | sort | uniq > {output.gene}")
        shell("cut -f 8  {input.bed} | sort | uniq > {output.gene}")
        shell("python {params.geneGSEA} --annotation {output.gene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


rule PromoterEnrichPlot:
    input:
        kegg = IN_PATH + "/population/Annotation/Enrichment/Promoter/KEGG_2019_Human..enrichr.reports.txt",
    output:
        kegg = IN_PATH + "/population/Annotation/Enrichment/Promoter/Plots/KEGG_2019_Human..enrichr.reports.pdf",
    params:
        geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
        indir = IN_PATH + "/population/Annotation/Enrichment/Promoter/",
        outdir = IN_PATH + "/population/Annotation/Enrichment/Promoter/Plots/",
        width = 8,
        height = 3,
        selectNum = 12,
    log:
        IN_PATH + "/log/CommonEnrichPlot.log",
    run:
        # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
        # print(cmd)
        # os.system(cmd)
        files = os.listdir(params.indir)
        for f in files:
            if f.endswith("enrichr.reports.txt"):
                fname = params.indir + f
                pdf = params.outdir + f.rstrip("txt") + "pdf"
                cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
                os.system(cmd)


rule PromoterenrichGeneGeno:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq.bed",
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        omim = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Disease..enrichr.reports.txt", 
        omim_expand = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Expanded..enrichr.reports.txt",
        clinvar = IN_PATH + "/population/Annotation/Enrichment/Promoter/ClinVar_2019..enrichr.reports.txt",
        gwas = IN_PATH + "/population/Annotation/Enrichment/Promoter/GWAS_Catalog_2019..enrichr.reports.txt",
        biobank = IN_PATH + "/population/Annotation/Enrichment/Promoter/UK_Biobank_GWAS_v1..enrichr.reports.txt",
    output:
        gene = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Disease_gene_tag_genotypes.txt", 
        gene_expand = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Expand_gene_tag_genotypes.txt", 
        clinvar = IN_PATH + "/population/Annotation/Enrichment/Promoter/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/Promoter/GWAS_gene_tag_genotypes.txt", 
        biobank = IN_PATH + "/population/Annotation/Enrichment/Promoter/Biobank_gene_tag_genotypes.txt", 
    log:
        IN_PATH + "/log/PromoterenrichGeneGeno.log",
    params:
        DiseaseGeneTagGeno = SRC_DIR + "/DiseaseGeneTagGeno.py",
    run:
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.omim} --annotation {input.anno} --genotype {input.geno} --out {output.gene} > {log} 2>&1")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.omim_expand} --annotation {input.anno} --genotype {input.geno} --out {output.gene_expand} 2>>{log}")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.clinvar} --annotation {input.anno} --genotype {input.geno} --out {output.clinvar} 2>>{log}")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.gwas} --annotation {input.anno} --genotype {input.geno} --out {output.gwas} 2>>{log}")
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.biobank} --annotation {input.anno} --genotype {input.geno} --out {output.biobank} 2>>{log}")


rule cosmicOverlapPromoter:
    input:
        cosmic = "/home/wuzhikun/database/COSMIC/cancer_gene_census.gene",
        gene = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq_gene.txt",
    output:
        cosmic = IN_PATH + "/population/Annotation/Enrichment/Promoter/COSMIC_gene_overlap_report.txt", 
    run:
        shell("""echo "COSMIC\tCancer\tGene" > {output.cosmic}""")
        cmd = """cat %s %s | sort | uniq -c | sort -k 1nr | sed 's/      //g' | grep '^2' | cut -f 2 -d ' ' | awk '{print "COSMIC\tCancer\t"$0}' >> %s""" % (input.gene, input.cosmic, output.cosmic)
        os.system(cmd)


rule cosmicGenoPromoter:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_nocoding_uniq.bed",
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        cosmic = IN_PATH + "/population/Annotation/Enrichment/Promoter/COSMIC_gene_overlap_report.txt",  
    output:
        cosmic = IN_PATH + "/population/Annotation/Enrichment/Promoter/COSMIC_Disease_gene_tag_genotypes.txt",
    log:
        IN_PATH + "/log/cosmicGenoPromoter.log",
    params:
        DiseaseGeneTagGeno = SRC_DIR + "/DiseaseGeneTagGeno.py",
    run:
        shell("python {params.DiseaseGeneTagGeno} --enrich {input.cosmic} --annotation {input.anno} --genotype {input.geno} --out {output.cosmic} > {log} 2>&1")





rule diseaseGenePromoter:
    input:
        gene_expand = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Expand_gene_tag_genotypes.txt", 
        # clinvar = IN_PATH + "/population/Annotation/Enrichment/Promoter/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/Promoter/GWAS_gene_tag_genotypes.txt", 
        # biobank = IN_PATH + "/population/Annotation/Enrichment/Promoter/Biobank_gene_tag_genotypes.txt", 
        cosmic = IN_PATH + "/population/Annotation/Enrichment/Promoter/COSMIC_Disease_gene_tag_genotypes.txt",
    output:
        gene = IN_PATH + "/population/Annotation/Enrichment/Promoter/disease_gene_circos.txt",
    params:
        geneDisease = SRC_DIR + "/geneDisease.py",
    log:
        IN_PATH + "/log/diseaseGenePromoter.log",
    run:
        Databases = ",".join([input.gene_expand, input.gwas, input.cosmic])
        Names = "OMIM,GWAS,COSMIC" 
        shell("python {params.geneDisease} --database {Databases} --name {Names} --out {output.gene} --feature CDS > {log} 2>&1")
    



rule diseaseCategooryPromoter:
    input:
        ommin = IN_PATH + "/population/Annotation/Enrichment/Promoter/OMIM_Expand_gene_tag_genotypes.txt", 
        # clinvar = IN_PATH + "/population/Annotation/Enrichment/Promoter/Clinvar_gene_tag_genotypes.txt", 
        gwas = IN_PATH + "/population/Annotation/Enrichment/Promoter/GWAS_gene_tag_genotypes.txt", 
        # biobank = IN_PATH + "/population/Annotation/Enrichment/Promoter/Biobank_gene_tag_genotypes.txt", 
        cosmic = IN_PATH + "/population/Annotation/Enrichment/Promoter/COSMIC_Disease_gene_tag_genotypes.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        cat = IN_PATH + "/population/Annotation/Enrichment/Promoter/Annotation_database_category.txt", 
    params:
        diseaseGeneCatrgory = SRC_DIR + "/diseaseGeneCatrgory.py",
        number = 405,
        name = "OMIM,COSMIC,GWAS",
    log:
        IN_PATH + "/log/diseaseCategooryPromoter.log",
    run:
        Database = ",".join([input.ommin, input.cosmic, input.gwas])
        shell("python {params.diseaseGeneCatrgory} --annotation {Database} --name {params.name} --number {params.number} --category {input.category} --out {output.cat} > {log} 2>&1")



rule diseaseCategooryPromoterPlot:
    input:
        cat = IN_PATH + "/population/Annotation/Enrichment/Promoter/Annotation_database_category.txt", 
    output:
        pdf = IN_PATH + "/population/Annotation/Enrichment/Promoter/Annotation_database_category.pdf", 
    params:
        diseaseCategory = SCRIPT_DIR + "/diseaseCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/diseaseCategooryPromoterPlot.log",
    run:
        shell("Rscript {params.diseaseCategory} --input {input.cat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")







#########################################################################################################









##########################################################################################################################################
# rule DELAllEnrich:
#     input:
#         DEL = IN_PATH + "/population/Annotation/Sample_DEL_INS_genepred_overlap.bed",
#     output:
#         DEL = IN_PATH + "/population/Annotation/Sample_DEL_INS_ALL_genepred_overlap_gene.xls",
#         kegg = IN_PATH + "/population/Annotation/Enrichment/DEL_INS_ALL/KEGG_2019_Human..enrichr.reports.txt", 
#     params:
#         outdir = IN_PATH + "/population/Annotation/Enrichment/DEL_INS_ALL",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/DELEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("cut -f 8 {input.DEL} | sort | uniq > {output.DEL}")
#         shell("python {params.geneGSEA} --annotation {output.DEL} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")



# rule DELAllEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/DEL_INS_ALL/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/DEL_INS_ALL/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/Annotation/Enrichment/DEL_INS_ALL/",
#         outdir = IN_PATH + "/population/Annotation/Enrichment/DEL_INS_ALL/Plots/",
#         width = 8,
#         height = 3,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)






##################################################################################################



########################################## AnnotSV annotation ######################################

rule SVAnnoAll:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        temp = IN_PATH + "/population/Merge/Sample_SV_common_temp.vcf",
        vcf = IN_PATH + "/population/Annotation/Sample_common_SV.tsv",
    params:
        AnnotSV = config["AnnotSV"],
        bedtools = config["bedtools"],
        annotsv_overlap = config["annotsv_overlap"],
        promoterSize = config["promoterSize"],
        outdir = IN_PATH + "/population/Annotation",
    log:
        IN_PATH + "/log/SVAnnoAll.log",
    run:
        ### export ANNOTSV=/home/wuzhikun/anaconda3/envs/NanoSV/bin/AnnotSV_2.1
        shell("cut -f 1-10 {input.vcf} > {output.temp}")
        shell("{params.AnnotSV} -bedtools {params.bedtools}   -SVinputFile  {output.temp} -genomeBuild GRCh38 -outputDir {params.outdir} -outputFile Sample_common_SV  -overlap {params.annotsv_overlap} -promoterSize {params.promoterSize} -reciprocal yes -typeOfAnnotation split > {log} 2>&1")

# rule AnnoStat:
#     input:
#         sv = IN_PATH + "/population/Annotation/Sample_common_SV.tsv",
#     output:
#         sv = IN_PATH + "/population/Annotation/Sample_SV_class.xls",
#     params:
#         AnnoSVClass = SRC_DIR + "/AnnoSVClass.py",
#     run:
#         shell("python {params.AnnoSVClass} --input {input.sv} --out {output.sv}")


rule PopAnnoModify:
    input:
        vcf = IN_PATH + "/population/Annotation/Sample_common_SV.tsv",
    output:
        vcf = IN_PATH + "/population/Annotation/Sample_common_SV_modify.tsv",
    params:
        annotsvIDModify = SRC_DIR + "/annotsvIDModify.py",
    threads:
        THREADS
    log:
        IN_PATH + "/log/PopAnnoModify.log",
    run:
        shell("python {params.annotsvIDModify} --annotation {input.vcf} --out {output.vcf} > {log} 2>&1")
#######################################################################################




# ################################ non-coding annotation #######################################

# def uniq_region(bed1, bed2, out_file):
#     Record = {}
#     bed_h = open(bed2, "r")
#     for line in bed_h:
#         line = line.strip()
#         lines = line.split("\t")
#         # record = "\t".join(lines[:3])
#         record = lines[3]
#         Record[record] = 1
#     bed_h.close()
#     bed_h = open(bed1, "r")
#     out_h = open(out_file, "w")
#     for line in bed_h:
#         line = line.strip()
#         lines = line.split("\t")
#         # record = "\t".join(lines[:3])
#         record = lines[3]
#         if record not in Record:
#             out_h.write("%s\n" % line)
#     bed_h.close()
#     out_h.close()


# rule NoncodingRegion:
#     input:
#         all = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
#         anno = IN_PATH + "/population/Annotation/Sample_common_SV_modify.tsv",
#     output:
#         gene = IN_PATH + "/population/Annotation/Sample_common_SV_gene.bed",
#         intergenic = IN_PATH + "/population/Annotation/Sample_common_SV_intergenic.bed",
#     log:
#         IN_PATH + "/log/NoncodingRegion.log",
#     run:
#         cmd = '''sed '1d' %s | awk '{print $2"\t"$3"\t"$4"\t"$1}' | uniq | sort -k 1,1n -k 2,2n | uniq > %s 2>%s''' % (input.anno, output.gene, log)
#         os.system(cmd)
#         # shell("cat {input.all} {output.gene} | sort -k 1,1n -k 2,2n | uniq -u > {output.intergenic} 2>>{log}")
#         uniq_region(input.all, output.gene, output.intergenic)




# rule Popenhancer:
#     input:
#         bed = IN_PATH + "/population/Annotation/Sample_common_SV_intergenic.bed",
#     output:
#         enhancer = IN_PATH + "/population/NonCoding/Sample_common_SV_enhancer.bed",
#         # promoter = IN_PATH + "/population/NonCoding/Sample_common_SV_promoter.bed",
#         regulator = IN_PATH + "/population/NonCoding/Sample_common_SV_regulator.bed",
#     params:
#         enhancerRegionGene = config["enhancerRegion"],
#         promoterRegion = config["promoterRegion"],
#         regulatorRegion = config["regulatorRegion"],
#     log:
#         IN_PATH + "/log/Popenhancer.log",
#     run:
#         shell("bedtools intersect -a {input.bed}  -b {params.enhancerRegionGene}  -wa -wb > {output.enhancer} 2>{log}")
#         # shell("bedtools intersect -a {input.bed}  -b {params.promoterRegion}  -wa -wb > {output.promoter} 2>>{log}")
#         shell("bedtools intersect -a {input.bed}  -b {params.regulatorRegion}  -wa -wb > {output.regulator} 2>>{log}")


# rule enhancerTag:
#     input:
#         enhancer = IN_PATH + "/population/NonCoding/Sample_common_SV_enhancer.bed",
#         regulator = IN_PATH + "/population/NonCoding/Sample_common_SV_regulator.bed",
#     output:
#         enhancer = IN_PATH + "/population/NonCoding/Sample_common_SV_enhancer_tag.txt",
#         regulator = IN_PATH + "/population/NonCoding/Sample_common_SV_regulator_tag.txt",
#     run:
#         shell("cut -f 4 {input.enhancer} | sort | uniq > {output.enhancer}")
#         shell("cut -f 4 {input.regulator} | sort | uniq > {output.regulator}")


# ##################################################################################################




# ################################# all enhancer and promoter  ########################
# rule Popenhancer2:
#     input:
#         # bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
#         bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP_point.bed",
#     output:
#         enhancer = IN_PATH + "/population/Annotation/Sample_common_SV_enhancer_all.bed",
#         regulator = IN_PATH + "/population/Annotation/Sample_common_SV_regulator_all.bed",
#     params:
#         enhancerRegionGene = config["enhancerRegion"],
#         regulatorRegion = config["regulatorRegionGene"],
#     log:
#         IN_PATH + "/log/Popenhancer2.log",
#     run:
#         shell("bedtools intersect -a {input.bed}  -b {params.enhancerRegionGene}  -wa -wb > {output.enhancer} 2>{log}")
#         shell("bedtools intersect -a {input.bed}  -b {params.regulatorRegion}  -wa -wb > {output.regulator} 2>>{log}")






# rule enhancerTag2:
#     input:
#         enhancer = IN_PATH + "/population/Annotation/Sample_common_SV_enhancer_all.bed",
#         regulator = IN_PATH + "/population/Annotation/Sample_common_SV_regulator_all.bed",
#     output:
#         enhancer = IN_PATH + "/population/Annotation/Sample_common_SV_enhancer_tag_all.txt",
#         regulator = IN_PATH + "/population/Annotation/Sample_common_SV_regulator_tag_all.txt",
#     run:
#         shell("cut -f 4 {input.enhancer} | sort | uniq > {output.enhancer}")
#         shell("cut -f 4 {input.regulator} | sort | uniq > {output.regulator}")




# rule Popenhancer3:
#     input:
#         enhancer = IN_PATH + "/population/Annotation/Sample_common_SV_enhancer_all.bed",
#         regulator = IN_PATH + "/population/Annotation/Sample_common_SV_regulator_all.bed",
#     output:
#         enhancer = IN_PATH + "/population/Annotation/Sample_common_SV_enhancer_all-2.bed",
#         regulator_temp = temp(IN_PATH + "/population/Annotation/Sample_common_SV_regulator_all_temp.bed"),
#         regulator = IN_PATH + "/population/Annotation/Sample_common_SV_regulator_all-2.bed",
#     log:
#         IN_PATH + "/log/Popenhancer3.log",
#     run:
#         cmd = """ awk '{print $0"\tEnhancer\tEnhancer"}' %s > %s """ % (input.enhancer, output.enhancer)
#         os.system(cmd)
#         shell("cut -f 1-8 {input.regulator} > {output.regulator_temp}")
#         cmd = """ awk '{print $0"\tPromoter\tPromoter"}' %s > %s """ % (output.regulator_temp, output.regulator)
#         os.system(cmd)


# ###################################################################################



############################### Annotation for special tags ##########################

rule SpecialEnrich:
    input:
        gene = IN_PATH + "/population/bed/CellExclude/Sample_common_DEL_INS_INV_special_anno_gene.txt",
    output:
        kegg = IN_PATH + "/population/bed/CellExclude/Enrichment/KEGG_2019_Human..enrichr.reports.txt", 
    params:
        outdir = IN_PATH + "/population/bed/CellExclude/Enrichment",
        geneGSEA = SRC_DIR + "/geneGSEA.py",
        EnrichLibrary = config["EnrichLibrary"],
        number = 1,
    log:
        IN_PATH + "/log/SpecialEnrich.log",
    run:
        ### mu01
        ### --number {params.number} 
        shell("python {params.geneGSEA} --annotation {input.gene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


rule SpecialEnrichPlot:
    input:
        kegg = IN_PATH + "/population/bed/CellExclude/Enrichment/KEGG_2019_Human..enrichr.reports.txt", 
    output:
        kegg = IN_PATH + "/population/bed/CellExclude/Enrichment/Plots/KEGG_2019_Human..enrichr.reports.pdf",
    params:
        geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
        indir = IN_PATH + "/population/bed/CellExclude/Enrichment/",
        outdir = IN_PATH + "/population/bed/CellExclude/Enrichment/Plots/",
        width = 8,
        height = 4,
        selectNum = 20,
    log:
        IN_PATH + "/log/SpecialEnrichPlot.log",
    run:
        # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
        # print(cmd)
        # os.system(cmd)
        files = os.listdir(params.indir)
        for f in files:
            if f.endswith("enrichr.reports.txt"):
                fname = params.indir + f
                pdf = params.outdir + f.rstrip("txt") + "pdf"
                cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
                os.system(cmd)
##########################################################################################




################################## categories of annotations ##############################
# rule annoCat:
#     input:
#         cat = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
#         anno = IN_PATH + "/population/Annotation/Sample_common_SV_modify.tsv",
#     output:
#         common = IN_PATH + "/population/Annotation/Sample_common_SV_anno_Common.xls",
#         major = IN_PATH + "/population/Annotation/Sample_common_SV_anno_Major.xls",
#         poly = IN_PATH + "/population/Annotation/Sample_common_SV_anno_Poly.xls",
#         single = IN_PATH + "/population/Annotation/Sample_common_SV_anno_Single.xls",
#     params:
#         outPrefix = IN_PATH + "/population/Annotation/Sample_common_SV_anno",
#         AnnoCategory = SRC_DIR + "/AnnoCategory.py",
#     log:
#         IN_PATH + "/log/annoCat.log"
#     run:
#         shell("python {params.AnnoCategory} --category {input.cat} --annotation {input.anno} --out {params.outPrefix} > {log} 2>&1")


############################################################################################




###################################### Annotation stats  #######################################
# rule GeneAnnoStat:
#     input:
#         anno = IN_PATH + "/population/Annotation/Sample_common_SV_anno_{category}.xls",
#     output:
#         stat = IN_PATH + "/population/Annotation/Sample_common_SV_gene_stat_{category}.xls",
#     params:
#         parseAnnotSV = SRC_DIR + "/parseAnnotSV.py",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/GeneAnnoStat_{category}.log",
#     run:
#         shell("python {params.parseAnnotSV} --annotation {input.anno} --out {output.stat} >{log} 2>&1")


# rule GeneSumStat:
#     input:
#         stat = IN_PATH + "/population/Annotation/Sample_common_SV_gene_stat_{category}.xls",
#     output:
#         location = IN_PATH + "/population/Annotation/Sample_common_SV_gene_stat_location_{category}.xls",
#         gene = IN_PATH + "/population/Annotation/Sample_common_SV_gene_stat_gene_{category}.xls",
#     params:
#         # annoGeneStats = SRC_DIR + "/annoGeneStats.py",
#         annoGeneStatsMatrix = SRC_DIR + "/annoGeneStatsMatrix.py",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/GeneSumStat_{category}.log",
#     run:
#         shell("python {params.annoGeneStatsMatrix} --stats {input.stat} --location {output.location} --gene {output.gene} >{log} 2>&1")



# rule GeneEnrichment:
#     input:
#         gene = IN_PATH + "/population/Annotation/Sample_common_SV_gene_stat_gene_{category}.xls",
#     output:
#         gene = IN_PATH + "/population/Annotation/Enrichment/Sample_common_SV_exon_gene_{category}.xls",
#     run:
#         cmd = "awk '{if($4==1){print $0}}' %s > %s" % (input.gene, output.gene)
#         os.system(cmd)


# rule GeneEnrich:
#     input:
#         gene = IN_PATH + "/population/Annotation/Enrichment/Sample_common_SV_exon_gene_{category}.xls",
#     output:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/{category}/KEGG_2019_Human..enrichr.reports.txt", 
#     params:
#         outdir = IN_PATH + "/population/Annotation/Enrichment/{category}",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/GeneEnrichment_{category}.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("python {params.geneGSEA} --annotation {input.gene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


# rule GeneEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/{category}/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/Annotation/Enrichment/{category}/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/Annotation/Enrichment/",
#         outdir = IN_PATH + "/population/Annotation/Enrichment/Plots/",
#         width = 8,
#         height = 4,
#         selectNum = 20,
#     log:
#         IN_PATH + "/log/GeneEnrichPlot_{category}.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)
##########################################################################################################




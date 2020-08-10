########################################## WGS 17795 #################################
rule ExtractBed:
    input:
        bed = "/home/wuzhikun/database/WGS17795/Build38.public.v2.filter.bedpe",
    output:
        DEL = "/home/wuzhikun/database/WGS17795/Build38_DEL.bed",
        INS = "/home/wuzhikun/database/WGS17795/Build38_INS.bed",
        DUP = "/home/wuzhikun/database/WGS17795/Build38_DUP.bed",
        INV = "/home/wuzhikun/database/WGS17795/Build38_INV.bed",
    params:
        ExtractBed = SRC_DIR + "/ExtractBed.py",
        outPrefix = "/home/wuzhikun/database/WGS17795/Build38",
    log:
        IN_PATH + "/log/ExtractBed.log"
    run:
        shell("python {params.ExtractBed} --bed {input.bed} --outPrefix {params.outPrefix} > {log} 2>&1")


######################################################################################




############################################################################################
rule SVExclude:
    input:
        allList = IN_PATH + "/preStudies/Cell2019_Table_S1.txt",
        Del = IN_PATH + "/preStudies/Cell2019_SV_dechr_DEL.bed",
        Ins = IN_PATH + "/preStudies/Cell2019_SV_dechr_INS_newEnd.bed",
        Inv = IN_PATH + "/preStudies/Cell2019_SV_dechr_INV.bed",
    output:
        Del = IN_PATH + "/preStudies/Cell2019_SV_dechr_DEL_exclude_AKHX.bed",
        Ins = IN_PATH + "/preStudies/Cell2019_SV_dechr_INS_exclude_AKHX.bed",
        Inv = IN_PATH + "/preStudies/Cell2019_SV_dechr_INV_exclude_AKHX.bed",
    params:
        SVListExclude = SRC_DIR + "/SVListExclude.py",
        exclude = "AK1,HX1",
    log:
        IN_PATH + "/log/SVExclude.log"
    run:
        shell("python {params.SVListExclude} --list {input.allList} --bed {input.Del} --out {output.Del} --exclude {params.exclude} > {log} 2>&1")
        shell("python {params.SVListExclude} --list {input.allList} --bed {input.Ins} --out {output.Ins} --exclude {params.exclude} >> {log} 2>&1")
        shell("python {params.SVListExclude} --list {input.allList} --bed {input.Inv} --out {output.Inv} --exclude {params.exclude} >> {log} 2>&1")




# rule dbVarInfor:
#     input:
#         vcf = config["dbVar"],
#     output:
#         xls = IN_PATH + "/database/dbVar/GRCh38.variant_call.xls",
#     params:
#         dbVarInfor = SRC_DIR + "/dbVarInfor.py",
#     log:
#         IN_PATH + "/log/dbVarInfor.log"
#     run:
#         shell("python {params.dbVarInfor} --vcf {input.vcf} --out {output.xls} >{log} 2>&1")


# rule phase3SVInfor:
#     input:
#         vcf = config["phase3SV"],
#     output:
#         xls = IN_PATH + "/database/phase3SV/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.xls",
#     params:
#         phase3SVInfor = SRC_DIR + "/phase3SVInfor.py",
#     log:
#         IN_PATH + "/log/phase3SVInfor.log"
#     run:
#         shell("python {params.phase3SVInfor} --vcf {input.vcf} --out {output.xls} >{log} 2>&1")


# rule DGVInfor:
#     input:
#         DGV = config["DGV"],
#     output:
#         xls = IN_PATH + "/database/DGV/GRCh38_hg38_variants_2016-08-31.xls",
#     params:
#         DGVInfor = SRC_DIR + "/DGVInfor.py",
#     log:
#         IN_PATH + "/log/DGVInfor.log"
#     run:
#         shell("python {params.DGVInfor} --input {input.DGV} --out {output.xls} >{log} 2>&1")



rule DGVType:
    input:
        DGV = config["DGV"],
    output:
        DEL = IN_PATH + "/database/DGV/type/GRCh38_hg38_variants_DEL.bed",
        INS = IN_PATH + "/database/DGV/type/GRCh38_hg38_variants_INS.bed",
        DUP = IN_PATH + "/database/DGV/type/GRCh38_hg38_variants_DUP.bed",
        INV = IN_PATH + "/database/DGV/type/GRCh38_hg38_variants_INV.bed",
    params:
        DGVType = SRC_DIR + "/DGVType.py",
        outPrefix = IN_PATH + "/database/DGV/type/GRCh38_hg38_variants",
    log:
        IN_PATH + "/log/DGVType.log"
    run:
        shell("python {params.DGVType} --input {input.DGV} --outPrefix {params.outPrefix} > {log} 2>&1")




# rule enhancerGene:
#     input:
#         enhancer = config["enhancerRegion"],
#         ensGene_gtf95 = config["ensGene_gtf95"],
#     output:
#         enhancer = config["enhancerRegionGene"],
#     params:
#         enhancerNearGenes = SRC_DIR + "/enhancerNearGenes.py",
#     log:
#         IN_PATH + "/log/enhancerGene.log"
#     run:
#         shell("python {params.enhancerNearGenes} --enhancer {input.enhancer} --genepred {input.ensGene_gtf95} --out {output.enhancer} > {log} 2>&1")

# rule regulatorGene:
#     input:
#         regulator = config["regulatorRegion"],
#         ensGene_gtf95 = config["ensGene_gtf95"],
#     output:
#         regulator = config["regulatorRegionGene"],
#     params:
#         enhancerNearGenes = SRC_DIR + "/enhancerNearGenes.py",
#     log:
#         IN_PATH + "/log/regulatorGene.log"
#     run:
#         shell("python {params.enhancerNearGenes} --enhancer {input.regulator} --genepred {input.ensGene_gtf95} --out {output.regulator} > {log} 2>&1")

###################################################################################################





##################################### GWAS catalog ################################################
rule GWASCatalog:
    input:
        association = "/home/wuzhikun/database/GWAS/gwas_catalog_v1.0.2-associations_e98_r2020-03-08.tsv",
    output:
        disease = "/home/wuzhikun/database/GWAS/gwas_disease_gene.xls",
        gene = "/home/wuzhikun/database/GWAS/gwas_gene_disease.xls",
    params:
        GWASCatalog = SRC_DIR + "/GWASCatalog.py",
    log:
        IN_PATH + "/log/GWASCatalog.log"
    run:
        shell("python {params.GWASCatalog} --catalog {input.association} --disease {output.disease} --gene {output.gene} > {log} 2>&1")

####################################################################################################



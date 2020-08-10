

#################################### Merge all SVs ##########################
rule MergeSVALL:
    input:
        vcf = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=SAMPLES),
    output:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common.vcf",
    params:
        metafile = config["metafile"],
        clique_maxflow_SV = SRC_DIR + "/clique_maxflow_SV.py",
        allele_freq = 0.2,
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/MergeSV.log",
    run:
        Files = ",".join(sorted(input.vcf))
        shell("python {params.clique_maxflow_SV} --workers {threads} -v {Files} -o {output.vcf} --allele_freq {params.allele_freq} > {log} 2>&1")


rule MergeSVFilt:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common.vcf",
    output:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    params:
        commonSVFilt = SRC_DIR + "/commonSVFilt.py",
    log:
        IN_PATH + "/log/MergeSVFilt.log",
    run:
        shell("python {params.commonSVFilt} --vcf {input.vcf} --out {output.vcf} > {log} 2>&1")



rule commonSVTag:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        tag = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    params:
        CombineSVTag = SRC_DIR + "/CombineSVTag.py",
    log:
        IN_PATH + "/log/SVTag.log"
    run:
        shell("python {params.CombineSVTag} --vcf {input.vcf} --out {output.tag} > {log} 2>&1")





rule SubPopulation:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_CN329.vcf",
    params:
        SubPopSVGeno = SRC_DIR + "/SubPopSVGeno.py",
    log:
        IN_PATH + "/log/SubPopulation.log",
    run:
        shell("python {params.SubPopSVGeno} --vcf {input.vcf} --out {output.vcf} --sample CN > {log} 2>&1")




########################################################################################################





################################################ SV genotype ########################################


rule SVGenotype:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_numeric.txt",
    params:
        SVGenotype = SRC_DIR + "/SVGenotype.py",
    log:
        IN_PATH + "/log/SVGenotype.log", 
    run:
        shell("python {params.SVGenotype} --vcf {input.vcf} --out {output.geno} --frequency {output.freq} >{log} 2>&1")




rule SVFrequencyFilt:
    input:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype_filt.txt",
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
    params:
        SVGenoFreqFilt = SRC_DIR + "/SVGenoFreqFilt.py",
        freqThreshold = 0.05,
    log:
        IN_PATH + "/log/SVFrequencyFilt.log", 
    run:
        shell("python {params.SVGenoFreqFilt} --genotype {input.geno} --freqThreshold {params.freqThreshold} --out {output.geno} --frequency {output.freq} > {log} 2>&1")





rule AlleleFrequencyStat:
    input:
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
    output:
        stat = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency_stat.xls",
    params:
        TypeAlleleFrequency = SRC_DIR + "/TypeAlleleFrequency.py",
    log:
        IN_PATH + "/log/AlleleFrequencyStat.log", 
    run:
        shell("python {params.TypeAlleleFrequency} --input {input.freq} --out {output.stat} > {log} 2>&1")


rule AlleleFrequencyPlot:
    input:
        stat = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency_stat.xls",
    output:
        pdf = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency_stat.pdf",
    params:
        TypeFreqBar = SCRIPT_DIR + "/TypeFreqBar.R",
        width = 6, 
        height = 4,
    log:
        IN_PATH + "/log/AlleleFrequencyPlot.log", 
    run:
        shell("Rscript {params.TypeFreqBar} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
##############################################################################################




##################################### change format #################################
rule Meta2PED:
    input:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype_filt.txt",
    output:
        Ped = IN_PATH + "/population/plink/Sample_SV_geno.ped",
        Map = IN_PATH + "/population/plink/Sample_SV_geno.map",
    params:
        # Triometa2ped = SRC_DIR + "/Triometa2ped.py",
        meta2ped = SRC_DIR + "/meta2ped.py",
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/Meta2PED.log", 
    run:
        ### just get autosomal SV
        shell("python {params.meta2ped} --meta {params.metafile} --genotype {input.geno} --map {output.Map} --ped {output.Ped} >{log} 2>&1")



rule ped2bed:
    input:
        Ped = IN_PATH + "/population/plink/Sample_SV_geno.ped",
        Map = IN_PATH + "/population/plink/Sample_SV_geno.map",
    output:
        bed = IN_PATH + "/population/plink/Sample_SV_geno.bed",
        bim = IN_PATH + "/population/plink/Sample_SV_geno.bim",
        fam = IN_PATH + "/population/plink/Sample_SV_geno.fam",
    # params:
    #     plink = config["plink"],
    log:
        IN_PATH + "/log/ped2bed.log", 
    run:
        ### plink run in mu01
        pedPrefix = input.Ped.rstrip(".ped")
        shell("plink --file {pedPrefix}  --make-bed --out  {pedPrefix} >{log} 2>&1")


#######################################################################################################



###################################### category genotype ###########################
    

rule genoCat:
    input:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        singleton = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Singleton.txt",
        rare = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Rare.txt",
        low = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Low.txt",
        common = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Common.txt",
    params:
        outPrefix = IN_PATH + "/population/genotype/Category/Sample_SV_genotype",
    run:
        genotype_category(input.category, input.geno, params.outPrefix)



rule CatHeterHomoRatio:
    input:
        singleton = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Singleton.txt",
        rare = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Rare.txt",
        low = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Low.txt",
        common = IN_PATH + "/population/genotype/Category/Sample_SV_genotype_Common.txt",
    output:
        singleton = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Singleton.txt",
        rare = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Rare.txt",
        low = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Low.txt",
        common = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Common.txt",
    params:
        HeterHomoRatio = SRC_DIR + "/HeterHomoRatioAll.py",
    log:
        IN_PATH + "/log/CatHeterHomoRatio.log", 
    run:
        shell("python {params.HeterHomoRatio} --genotype {input.singleton} --out {output.singleton} > {log} 2>&1")
        shell("python {params.HeterHomoRatio} --genotype {input.rare} --out {output.rare} >> {log} 2>&1")
        shell("python {params.HeterHomoRatio} --genotype {input.low} --out {output.low} >> {log} 2>&1")
        shell("python {params.HeterHomoRatio} --genotype {input.common} --out {output.common} >> {log} 2>&1")


rule CatHeterHomoRatioSum:
    input:
        singleton = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Singleton.txt",
        rare = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Rare.txt",
        low = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Low.txt",
        common = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_Common.txt",
    output:
        category = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_category_summary.txt",
    run:
        shell("head -n 1 {input.singleton} | sed 's/Chr/Category/' > {output.category}")
        shell("tail -n 1 {input.singleton} | sed 's/Total/Singleton/' >> {output.category}")
        shell("tail -n 1 {input.rare} | sed 's/Total/Rare/' >> {output.category}")
        shell("tail -n 1 {input.low} | sed 's/Total/Low/' >> {output.category}")
        shell("tail -n 1 {input.common} | sed 's/Total/Common/' >> {output.category}")


rule CatHeterHomoRatioSumPlot:
    input:
        category = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_category_summary.txt",
    output:
        pdf = IN_PATH + "/population/genotype/Category/Sample_SV_heter2homo_category_summary.pdf",
    params:
        CategoryBox = SCRIPT_DIR + "/CategoryBox.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/CatHeterHomoRatioSumPlot.log", 
    run:
        shell("Rscript {params.CategoryBox} --input {input.category} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
##################################################################################  





################################## Heter to Homo Ratio #######################################
rule ChrHeterHomoRatio:
    input:
        genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        ratio = IN_PATH + "/population/genotype/Sample_SV_heter2homo.txt",
    params:
        HeterHomoRatio = SRC_DIR + "/HeterHomoRatio.py",
    log:
        IN_PATH + "/log/HeterHomoRatio.log", 
    run:
        shell("python {params.HeterHomoRatio} --genotype {input.genotype} --out {output.ratio} > {log} 2>&1")


rule HeterHomoRatioPlot:
    input:
        ratio = IN_PATH + "/population/genotype/Sample_SV_heter2homo.txt",
    output:
        pdf = IN_PATH + "/population/genotype/Sample_SV_heter2homo.pdf",   
    params:
        HeterRatioScatter = SCRIPT_DIR + "/HeterRatioScatter.R",
        width = 7,
        height = 4,
    log:
        IN_PATH + "/log/HeterHomoRatioPlot.log", 
    run:
        shell("Rscript {params.HeterRatioScatter} --input {input.ratio} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")



rule TypeHeterHomoRatio:
    input:
        genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        ratio = IN_PATH + "/population/genotype/Sample_SV_type_heter2homo.txt",
    params:
        CategoryHeterHomo = SRC_DIR + "/CategoryHeterHomo.py",
    log:
        IN_PATH + "/log/TypeHeterHomoRatio.log", 
    run:
        shell("python {params.CategoryHeterHomo} --genotype {input.genotype} --category type --out {output.ratio} > {log} 2>&1")


rule TypeHeterHomoRatioPlot:
    input:
        ratio = IN_PATH + "/population/genotype/Sample_SV_type_heter2homo.txt",
    output:
        pdf = IN_PATH + "/population/genotype/Sample_SV_type_heter2homo.pdf",
    params:
        CategoryBox = SCRIPT_DIR + "/CategoryBox.R",
        width = 4, # 2.6,
        height = 4, #2.6,
    log:
        IN_PATH + "/log/TypeHeterHomoRatioPlot.log", 
    run:
        shell("Rscript {params.CategoryBox} --input {input.ratio} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")



# rule FeatureGenotype:
#     input:
#         genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
#         feature = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_enhancer_all.bed",
#     output:
#         geno = IN_PATH + "/population/genotype/feature/Sample_SV_genotype_{feature}.txt",
#     params:
#         FeatureGenotype = SRC_DIR + "/FeatureGenotype.py",
#         outPrefix = IN_PATH + "/population/genotype/feature/Sample_SV_genotype",
#     log:
#         IN_PATH + "/log/FeatureGenotype_{feature}.log", 
#     run:
#         shell("python {params.FeatureGenotype} --genotype {input.genotype} --annotation {input.feature} --outPrefix {params.outPrefix} > {log} 2>&1")


# rule featureHeter:
#     input:
#         genotype = IN_PATH + "/population/genotype/feature/Sample_SV_genotype_{feature}.txt",
#     output:
#         ratio = IN_PATH + "/population/genotype/feature/Sample_SV_heter_{feature}.txt",
#     params:
#         CategoryHeterHomo = SRC_DIR + "/CategoryHeterHomo.py",
#     log:
#         IN_PATH + "/log/featureHeter_{feature}.log", 
#     run:
#         shell("python {params.CategoryHeterHomo} --genotype {input.genotype} --category type --out {output.ratio} > {log} 2>&1")

# rule heterSum:
#     input:
#         ratio = expand(IN_PATH + "/population/genotype/feature/Sample_SV_heter_{feature}.txt", feature=FEATURES),
#     output:
#         ratio = IN_PATH + "/population/genotype/feature/All_feature_heter_retio.txt",
#     run:
#         RATIOS = sorted(input.ratio)
#         for i in range(len(RATIOS)):
#             ratio = RATIOS[i]
#             if i == 0:
#                 cmd1 = "sed -n '1p' %s > %s" % (ratio, output.ratio)
#                 os.system(cmd1)
#                 cmd = "sed -n '$p' %s >> %s" % (ratio, output.ratio)
#             else:
#                 cmd = "sed -n '$p' %s >> %s" % (ratio, output.ratio)
#             print(cmd)
#             os.system(cmd)



# rule FeatureHomoRatioPlot:
#     input:
#         ratio = IN_PATH + "/population/genotype/feature/All_feature_heter_retio.txt",
#     output:
#         pdf = IN_PATH + "/population/genotype/feature/All_feature_heter_retio.pdf",
#     params:
#         CategoryBox = SCRIPT_DIR + "/CategoryBox.R",
#         width = 2.6,
#         height = 2.6,
#     log:
#         IN_PATH + "/log/FeatureHomoRatioPlot.log", 
#     run:
#         shell("Rscript {params.CategoryBox} --input {input.ratio} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")



# rule LengthHeter:
#     input:
#         genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
#     output:
#         genotype = IN_PATH + "/population/genotype/Sample_SV_genotype_{SVtype}.txt",
#         ratio = IN_PATH + "/population/genotype/Sample_{SVtype}_heter2homo.txt",
#         pdf = IN_PATH + "/population/genotype/Sample_{SVtype}_heter2homo.pdf",
#     params:
#         CategoryHeterHomo = SRC_DIR + "/CategoryHeterHomo.py",
#         LengthHeterBox = SCRIPT_DIR + "/LengthHeterBox.R",
#         width = 5,
#         height = 4,
#     log:
#         IN_PATH + "/log/LengthHeter_{SVtype}.log", 
#     run:
#         shell("sed -n '1p' {input.genotype} > {output.genotype}")
#         shell("grep {wildcards.SVtype} {input.genotype}  >> {output.genotype}")
#         shell("python {params.CategoryHeterHomo} --genotype {output.genotype} --category length --out {output.ratio} > {log} 2>&1")
#         shell("Rscript {params.LengthHeterBox} --input {output.ratio} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")

        
# rule SVBurden:
#     input:
#         genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
#     output:
#         sex = IN_PATH + "/population/genotype/Sample_SV_burden_sex.txt",
#         age = IN_PATH + "/population/genotype/Sample_SV_burden_age.txt",
#     params:
#         SVBurden = SRC_DIR + "/SVBurden.py",
#         metafile = config["metafile"],
#     log:
#         IN_PATH + "/log/SVBurden.log", 
#     run:
#         shell("python {params.SVBurden} --genotype {input.genotype} --out {output.sex} --metafile {params.metafile} --category sex > {log} 2>&1")
#         shell("python {params.SVBurden} --genotype {input.genotype} --out {output.age} --metafile {params.metafile} --category age >> {log} 2>&1")


# rule SVBurdenPlot:
#     input:
#         burden = IN_PATH + "/population/genotype/Sample_SV_burden_sex.txt",
#     output:
#         burden = IN_PATH + "/population/genotype/Sample_SV_burden_sex.pdf",
#     params:
#         SVBurdenBox = SCRIPT_DIR + "/SVBurdenBox.R",
#         width = 6,
#         height = 3,
#     log:
#         IN_PATH + "/log/SVBurdenPlot.log", 
#     run:
#         shell("Rscript {params.SVBurdenBox} --input {input.burden} --pdf {output.burden} --width {params.width} --height {params.width} > {log} 2>&1")
###############################################################################################





###################################  genotype with MAF 0.1 ###################################
rule SVFrequencyFilt01:
    input:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype_filt_01.txt",
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency_01.txt",
    params:
        SVGenoFreqFilt = SRC_DIR + "/SVGenoFreqFilt.py",
        freqThreshold = 0.1,
    log:
        IN_PATH + "/log/SVFrequencyFilt01.log", 
    run:
        shell("python {params.SVGenoFreqFilt} --genotype {input.geno} --freqThreshold {params.freqThreshold} --out {output.geno} --frequency {output.freq} > {log} 2>&1")


rule Meta2PED01:
    input:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype_filt_01.txt",
    output:
        Ped = IN_PATH + "/population/plink/Sample_SV_geno_01.ped",
        Map = IN_PATH + "/population/plink/Sample_SV_geno_01.map",
    params:
        # Triometa2ped = SRC_DIR + "/Triometa2ped.py",
        meta2ped = SRC_DIR + "/meta2ped.py",
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/Meta2PED01.log", 
    run:
        ### just get autosomal SV
        shell("python {params.meta2ped} --meta {params.metafile} --genotype {input.geno} --map {output.Map} --ped {output.Ped} >{log} 2>&1")



rule ped2bed01:
    input:
        Ped = IN_PATH + "/population/plink/Sample_SV_geno_01.ped",
        Map = IN_PATH + "/population/plink/Sample_SV_geno_01.map",
    output:
        bed = IN_PATH + "/population/plink/Sample_SV_geno_01.bed",
        bim = IN_PATH + "/population/plink/Sample_SV_geno_01.bim",
        fam = IN_PATH + "/population/plink/Sample_SV_geno_01.fam",
    # params:
    #     plink = config["plink"],
    log:
        IN_PATH + "/log/ped2bed01.log", 
    run:
        ### plink run in mu01
        pedPrefix = input.Ped.rstrip(".ped")
        shell("plink --file {pedPrefix}  --make-bed --out  {pedPrefix} >{log} 2>&1")

###############################################################################################


########################### SV LD #####################################
rule SVLD:
    input:
        bed = IN_PATH + "/population/plink/Sample_SV_geno_01.bed",
        bim = IN_PATH + "/population/plink/Sample_SV_geno_01.bim",
        fam = IN_PATH + "/population/plink/Sample_SV_geno_01.fam",
    output:
        ld = IN_PATH + "/population/LD/Sample_SV.ld",
    params:
        inputPrefix = IN_PATH + "/population/plink/Sample_SV_geno",
        outprefix = IN_PATH + "/population/LD/Sample_SV",
        outdir = IN_PATH + "/population/LD/",
        ld_window_kb = 300000,
        ld_window_r2 = 0.2,
    log:
        IN_PATH + "/log/SVLD.log",
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        shell("plink --bfile {params.inputPrefix} --r2 --ld-window-kb {params.ld_window_kb}  --ld-window-r2 {params.ld_window_r2} --out {params.outprefix} > {log} 2>&1")


rule SVLDMatrix:
    input:
        Ped = IN_PATH + "/population/plink/Sample_SV_geno_01.ped",
        Map = IN_PATH + "/population/plink/Sample_SV_geno_01.map",
    output:
        bed = IN_PATH + "/population/plink/Chrs/Sample_SV_geno_{chr}.bed",
        bim = IN_PATH + "/population/plink/Chrs/Sample_SV_geno_{chr}.bim",
        fam = IN_PATH + "/population/plink/Chrs/Sample_SV_geno_{chr}.fam",
        Map = IN_PATH + "/population/plink/Chrs/Sample_SV_geno_{chr}.map",
        ld = IN_PATH + "/population/plink/Chrs/{chr}_matrix.ld",
    params:
        inputPrefix = IN_PATH + "/population/plink/Sample_SV_geno",
        outprefix = IN_PATH + "/population/plink/Chrs/Sample_SV_geno_{chr}",
        ldPrefix = IN_PATH + "/population/plink/Chrs/{chr}_matrix",
        outdir = IN_PATH + "/population/plink/Chrs/",
    log:
        IN_PATH + "/log/SVLDMatrix_{chr}.log",
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        shell("plink --file {params.inputPrefix} --chr {wildcards.chr}  --out {params.outprefix} > {log} 2>&1")
        cmd = "grep '^%s\t' %s > %s" % (wildcards.chr, input.Map, output.Map)
        os.system(cmd)
        shell("plink --bfile {params.outprefix}  --matrix --r2 --out {params.ldPrefix} >> {log} 2>&1")


rule RegionLD:
    input:
        Map = IN_PATH + "/population/plink/Chrs/Sample_SV_geno_{chr}.map",
        ld = IN_PATH + "/population/plink/Chrs/{chr}_matrix.ld",
    output:
        ld = IN_PATH + "/population/LD/{chr}_ld_value.xls",
        ld_win = IN_PATH + "/population/LD/{chr}_ld_value_win.xls",
    params:
        LDMatrixDistance = SRC_DIR + "/LDMatrixDistance.py",
        window = 1000,
        sliding = 0,
        lenThreshold = 1000000,
    log:
        IN_PATH + "/log/RegionLD_{chr}.log",
    run:
        shell("python {params.LDMatrixDistance} --ld {input.ld} --map {input.Map} --out {output.ld} --winOut {output.ld_win} --window {params.window} --sliding {params.sliding} --lenThreshold {params.lenThreshold} > {log} 2>&1")
        
########################################################################



######################### tree ########################################
rule Distance:
    input:
        geno = IN_PATH + "/population/genotype/Sample_SV_genotype_filt_01.txt",
    output:
        geno = IN_PATH + "/population/tree/Sample_SV_genotype_matrix.txt",
    params:
        genoMatrix = SRC_DIR + "/genoMatrix.py",
    log:
        IN_PATH + "/log/Distance.log",
    run:
        shell("python {params.genoMatrix} --genotype {input.geno} --out {output.geno} > {log} 2>&1")



rule Distance2:
    input:
        geno = IN_PATH + "/population/tree/Sample_SV_genotype_matrix.txt",
    output:
        distance = IN_PATH + "/population/tree/Sample_SV_genotype_distance.txt",
    params:
        simpleDistance = SCRIPT_DIR + "/simpleDistance.R",
    log:
        IN_PATH + "/log/Distance2.log",
    run:
        shell("Rscript {params.simpleDistance} --input {input.geno} --out {output.distance} > {log} 2>&1")
        
#######################################################################



################################### heter to homo ###############
# rule FunctionRegion:
#     input:
#         anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
#     output:
#         cds = IN_PATH + "/population/Annotation/Region/SV_cds_tag.txt",
#         exon = IN_PATH + "/population/Annotation/Region/SV_exon_tag.txt",
#         intron = IN_PATH + "/population/Annotation/Region/SV_intron_tag.txt",
#         utr = IN_PATH + "/population/Annotation/Region/SV_utr_tag.txt",
#         updown = IN_PATH + "/population/Annotation/Region/SV_updown_tag.txt",
#     run:
#         shell("grep CDS {input.anno} | cut -f 4 | sort | uniq > {output.cds}")
#         shell("grep Exon {input.anno} | cut -f 4 | sort | uniq > {output.exon}")
#         shell("grep Intron {input.anno} | cut -f 4 | sort | uniq > {output.intron}")
#         shell("grep UTR3 {input.anno} | cut -f 4 | sort | uniq > {output.utr}")
#         shell("grep UTR5 {input.anno} | cut -f 4 | sort | uniq >> {output.utr}")
#         shell("grep DownStream {input.anno} | cut -f 4 | sort | uniq > {output.updown}")
#         shell("grep UpStream {input.anno} | cut -f 4 | sort | uniq >> {output.updown}")


# rule CatAnno:
#     input:
#         cds = IN_PATH + "/population/Annotation/Region/SV_cds_tag.txt",
#         exon = IN_PATH + "/population/Annotation/Region/SV_exon_tag.txt",
#         intron = IN_PATH + "/population/Annotation/Region/SV_intron_tag.txt",
#         utr = IN_PATH + "/population/Annotation/Region/SV_utr_tag.txt",
#         updown = IN_PATH + "/population/Annotation/Region/SV_updown_tag.txt",
#         category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
#     output:
#         cds = IN_PATH + "/population/Annotation/Region/SV_cds_tag.txt",
#         exon = IN_PATH + "/population/Annotation/Region/SV_exon_tag.txt",
#         intron = IN_PATH + "/population/Annotation/Region/SV_intron_tag.txt",
#         utr = IN_PATH + "/population/Annotation/Region/SV_utr_tag.txt",
#         updown = IN_PATH + "/population/Annotation/Region/SV_updown_tag.txt",


rule HomoTags:
    input:
        genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        homo = IN_PATH + "/population/genotype/homo/Sample_SV_genotype_homo.txt",
        tag = IN_PATH + "/population/genotype/homo/Sample_SV_tag_homo.txt",
    params:
        homoGeno = SRC_DIR + "/homoGeno.py",
    log:
        IN_PATH + "/log/HomoTags.log",
    run:
        shell("python {params.homoGeno} --genotype {input.genotype} --tag {output.tag} --filt {output.homo} > {log} 2>&1")




rule HomoTagCat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        tag = IN_PATH + "/population/genotype/homo/Sample_SV_tag_homo.txt",
    output:
        category = IN_PATH + "/population/genotype/homo/Sample_SV_tag_homo_category.txt",
    run:
        target_category(input.category, input.tag, output.category)
        # homo_tag_freq(input.tag, output.category, 405)



rule HomoAnno:
    input:
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
        tag = IN_PATH + "/population/genotype/homo/Sample_SV_tag_homo.txt",
    output:
        anno = IN_PATH + "/population/genotype/homo/Sample_homo_SV_genepred_overlap.bed",
    run:
        target_tag_annotation(input.anno, input.tag, output.anno)



rule HomoAnnoStats:
    input:
        anno = IN_PATH + "/population/genotype/homo/Sample_homo_SV_genepred_overlap.bed",
        category = IN_PATH + "/population/genotype/homo/Sample_SV_tag_homo_category.txt",
    output:
        stat = IN_PATH + "/population/genotype/homo/Sample_homo_SV_annotation_stats.xls",
        ratio = IN_PATH + "/population/genotype/homo/Sample_homo_SV_annotation_stats_ratio.xls",
    params:
        CategoryAnnoStats = SRC_DIR + "/CategoryAnnoStats.py",
    log:
        IN_PATH + "/log/HomoAnnoStats.log",
    run:
        shell("python {params.CategoryAnnoStats} --category {input.category} --annotation {input.anno} --out {output.ratio} --stat {output.stat} > {log} 2>&1")

##################################################################



######################### Shuffle genotypes ##########################


rule shuffleGeno:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        genotype = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        stat = IN_PATH + "/population/genotype/shuffle/Sample_geno_{number}_stat.txt",
    params:
        randomGeno = SRC_DIR + "/randomGeno.py",
        repeat = 10,
    threads:
        THREADS
    log:
        IN_PATH + "/log/genoShuffle/shuffleGeno_{number}.log",
    run:
        shell("python {params.randomGeno} --genotype {input.genotype} --category {input.category} --out {output.stat} --number {wildcards.number} --repeat {params.repeat} > {log} 2>&1")



def reshape_catgegory_value(sum_file, out_file):
    in_h = open(sum_file, "r")
    headers = in_h.readline().strip().split("\t")
    cats = headers[1:]
    out_h = open(out_file, "w")
    out_h.write("Sample\tCategory\tNumber\n")
    for line in in_h:
        lines = line.strip().split("\t")
        number = lines[0]
        values = lines[1:]
        for i, j in zip(cats, values):
            out_h.write("%s\t%s\t%s\n" % (number, i, j))
    in_h.close()
    out_h.close()



rule shuffleGenoSum:
    input:
        stat = expand(IN_PATH + "/population/genotype/shuffle/Sample_geno_{number}_stat.txt", number=NUMBERS),
    output:
        stat = IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary.txt",
        reshape = IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary_reshape.txt",
    run:
        STATS = sorted(input.stat)
        for i in range(len(STATS)):
            stat = STATS[i]
            if i == 0:
                cmd = "cat %s > %s" % (stat, output.stat)
            else:
                cmd = "sed -n '2p' %s >> %s" % (stat, output.stat)
            os.system(cmd)

        reshape_catgegory_value(output.stat, output.reshape)


rule shuffleGenoSumPlot:
    input:
        reshape = IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary_reshape.txt",
    output:
        reshape = IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary_reshape-1.txt",
        pdf = IN_PATH + "/population/genotype/All_Sample_geno_shuffle_summary_reshape.pdf",
    params:
        shuffleSampleLine = SCRIPT_DIR + "/shuffleSampleLine.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/shuffleGenoSumPlot.log",
    run:
        shell("grep -v All {input.reshape} > {output.reshape}")
        shell("Rscript {params.shuffleSampleLine} --input {output.reshape} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")

######################################################################
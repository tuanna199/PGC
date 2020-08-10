
######################################## Statistics ###############################################
rule SVSatats:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        Type = IN_PATH + "/population/Category/Sample_type_freq.xls",
        poly = IN_PATH + "/population/Category/Sample_SV_poly.xls",
        sample = IN_PATH + "/population/Category/Sample_SV_types.xls",
        Common = IN_PATH + "/population/Category/Sample_SV_record_Common.vcf",
        Major = IN_PATH + "/population/Category/Sample_SV_record_Low.vcf",
        Poly = IN_PATH + "/population/Category/Sample_SV_record_Rare.vcf",
        Single = IN_PATH + "/population/Category/Sample_SV_record_Singleton.vcf",
    log:
        IN_PATH + "/log/SVSatatsAll.log"
    params:
        SVAlleleCalculateMod = SRC_DIR + "/SVAlleleCalculateMod.py",
        outPrefix = IN_PATH + "/population/Category/Sample_SV_record",
    run:
        shell("python {params.SVAlleleCalculateMod} --category {input.category} --vcf {input.vcf}  --type {output.Type} --poly {output.poly} --sample {output.sample} --outPrefix {params.outPrefix} >{log} 2>&1") 





rule typeFreqPlot:
    input:
        Type = IN_PATH + "/population/Category/Sample_type_freq.xls",
    output:
        pdf =  IN_PATH + "/population/Category/Sample_type_freq.pdf",
        pdf2 =  IN_PATH + "/population/Category/Sample_type_freq_percentage.pdf",
    params:
        SVFreqClass = SCRIPT_DIR + "/SVFreqClass.R",
        height = 4,
        width = 5,
    log:
        IN_PATH + "/log/typeFreqPlot.log"
    run:
        # shell("source activate Rmeta && Rscript {params.SVFreqClass} --input {input.Type} --pdf {output.pdf} --height {params.height} --width {params.width} --pdf2  {output.pdf2} >{log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --height %s --width %s --pdf2  %s >  %s 2>&1" % (params.SVFreqClass, input.Type, output.pdf, params.height, params.width, output.pdf2, log)
        print(cmd)
        os.system(cmd)


rule SampleFreqPlot:
    input:
        poly = IN_PATH + "/population/Category/Sample_SV_poly.xls",
    output:
        hist = IN_PATH + "/population/Category/Sample_SV_poly_hist.pdf",
        box = IN_PATH + "/population/Category/Sample_SV_poly_box.pdf",
    params:
        SVFreqBarMultiple = SCRIPT_DIR + "/SVFreqBarMultiple.R",
        width = 20,
        height = 4,
    log:
        IN_PATH + "/log/SampleFreqPlot.log"
    run:
        # shell("source activate Rmeta && Rscript {params.SVFreqBarMultiple} --input {input.poly} --hist {output.hist} --box {output.box} --width {params.width} --height {params.height} >{log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --hist %s --box %s --width %s --height %s > %s 2>&1" % (params.SVFreqBarMultiple, input.poly, output.hist, output.box, params.width, params.height, log)
        print(cmd)
        os.system(cmd)



rule SampleTypePlot:
    input:
        summary = IN_PATH + "/population/Category/Sample_SV_types.xls",
    output:
        hist = IN_PATH + "/population/Category/Sample_SV_types_hist.pdf",
        box = IN_PATH + "/population/Category/Sample_SV_types_box.pdf",
    params:
        SVTypeBarMultiple = SCRIPT_DIR + "/SVTypeBarMultiple.R",
        ### height and width of hist plot
        width = 8,
        height = 4,
    log:
        IN_PATH + "/log/SummaryPlot.log", 
    run:
        # shell("source activate Rmeta && Rscript {params.SVTypeBarMultiple} --input {input.summary} --hist {output.hist} --box {output.box} --width {params.width} --height {params.height} >{log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --hist %s --box %s --width %s --height %s > %s 2>&1" % (params.SVTypeBarMultiple, input.summary, output.hist, output.box, params.width, params.height, log)
        print(cmd)
        os.system(cmd)
#######################################################################################################


# ##################################### SV Tags ###############################################

# rule commonSVTag:
#     input:
#         vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
#     output:
#         tag = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
#     params:
#         CombineSVTag = SRC_DIR + "/CombineSVTag.py",
#     log:
#         IN_PATH + "/log/SVTag.log"
#     run:
#         shell("python {params.CombineSVTag} --vcf {input.vcf} --out {output.tag} > {log} 2>&1")



# rule categoryTag:
#     input:
#         vcf = IN_PATH + "/population/Category/Sample_SV_record_{category}.vcf",
#     output:
#         tag = IN_PATH + "/population/Category/Sample_SV_record_{category}_SV_tag.xls",
#     params:
#         CombineSVTag = SRC_DIR + "/CombineSVTag.py",
#     log:
#         IN_PATH + "/log/categoryTag_{category}.log"
#     run:
#         shell("python {params.CombineSVTag} --vcf {input.vcf} --out {output.tag} > {log} 2>&1")


# rule categoryTag2:
#     input:
#         tag = expand(IN_PATH + "/population/Category/Sample_SV_record_{category}_SV_tag.xls", category=CATEGORY),
#     output:
#         tag = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
#     params:
#         categoryTag = SRC_DIR + "/categoryTag.py",
#     log:
#         IN_PATH + "/log/categoryTag2.log"
#     run:
#         TAGS = ",".join(input.tag)
#         shell("python {params.categoryTag} --files {TAGS} --output {output.tag} > {log} 2>&1")
# ###############################################################################################












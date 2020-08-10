

#####################################  SV merge and filt #################################################################
rule ChangeFormat:
    input:
        nanosv =  IN_PATH + "/SVCall/NanoSV/minimap2/{sample}.vcf",
        nanovar = IN_PATH + "/SVCall/nanovar/{sample}/{sample}.nanovar.pass.vcf",
    output:
        nanosv =  IN_PATH + "/SVCall/Format/{sample}_nanosv.vcf",
        nanovar = IN_PATH + "/SVCall/Format/{sample}_nanovar.vcf",
    params:
        nanovar_change_format = SRC_DIR + "/nanovar_change_format.py",
        nanosv_change_format = SRC_DIR + "/nanosv_change_format.py",
    log:
        IN_PATH + "/log/{sample}.ChangeFormat.log", 
    run:
        shell("python {params.nanosv_change_format} --raw_vcf {input.nanosv}  --new_vcf {output.nanosv} --score 20 > {log} 2>&1")
        shell("python {params.nanovar_change_format} --raw_vcf {input.nanovar} --new_vcf {output.nanovar} >> {log} 2>&1")



rule MultipleFilt:
    input:
        sniffles = IN_PATH + "/SVCall/Sniffles/minimap2/{sample}.vcf",
        nanosv =  IN_PATH + "/SVCall/Format/{sample}_nanosv.vcf",
        nanovar = IN_PATH + "/SVCall/Format/{sample}_nanovar.vcf",
    output:
        sniffles = IN_PATH + "/SVCall/SVFilt/{sample}_sniffles.vcf",
        nanosv = IN_PATH + "/SVCall/SVFilt/{sample}_nanosv.vcf",
        nanovar = IN_PATH + "/SVCall/SVFilt/{sample}_nanovar.vcf",
    threads:
        THREADS
    params:
        SVFiltReadsMultiple = SRC_DIR + "/SVFiltReadsMultiple.py",
        ratioThreshold = 0.1,
        column = "Clean_total_base",
    log:
        IN_PATH + "/log/{sample}.MultipleFilt.log", 
    run:
        ### TypeList = ["INS", "DEL", "INV", "DUP"]
        ### CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
        shell("python {params.SVFiltReadsMultiple} --support 2 --vcf {input.sniffles} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.sniffles}  >{log} 2>&1")
        shell("python {params.SVFiltReadsMultiple} --support 2 --vcf {input.nanosv} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.nanosv}   >>{log} 2>&1")
        shell("python {params.SVFiltReadsMultiple} --support 2 --vcf {input.nanovar} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.nanovar}  >>{log} 2>&1")


###################################################################################################################################





################################################  merge SV ########################################################
rule MergeMultipleSV:
    input:
        sniffles = IN_PATH + "/SVCall/SVFilt/{sample}_sniffles.vcf",
        nanovar = IN_PATH + "/SVCall/SVFilt/{sample}_nanovar.vcf",
        nanosv = IN_PATH + "/SVCall/SVFilt/{sample}_nanosv.vcf",
    output:
        vcf = IN_PATH + "/SVCall/Multiple/{sample}/{sample}_common_SV.vcf",
    params:
        metafile = config["metafile"],
        clique_maxflow_SV = SRC_DIR + "/clique_maxflow_SV.py",
        allele_freq = 0.2,
    threads:
        THREADS
    log:
        IN_PATH + "/log/MergeMultipleSV_{sample}.log", 
    run:
        Files = ",".join([input.sniffles, input.nanovar, input.nanosv])
        shell("python {params.clique_maxflow_SV} -v {Files} -o {output.vcf} --allele_freq {params.allele_freq} > {log} 2>&1")


rule OriginalRecord:
    input:
        vcf = IN_PATH + "/SVCall/Multiple/{sample}/{sample}_common_SV.vcf",
        sniffles = IN_PATH + "/SVCall/SVFilt/{sample}_sniffles.vcf",
        nanovar = IN_PATH + "/SVCall/SVFilt/{sample}_nanovar.vcf",
    output:
        vcf = IN_PATH + "/SVCall/Final/{sample}_common_original.vcf",
    params:
        MergeSVOriginal = SRC_DIR + "/MergeSVOriginal.py",
    threads:
        THREADS
    log:
        IN_PATH + "/log/OriginalRecord_{sample}.log", 
    run:
        ### dominant
        shell("python {params.MergeSVOriginal} --merge {input.vcf} --sniffle {input.sniffles} --nanovar {input.nanovar} --out {output.vcf} --method dominant > {log} 2>&1")




################################################################################





# #####################################  SV merge and filt #################################################################
# rule ChangeFormat:
#     input:
#         nanosv =  IN_PATH + "/SVCall/NanoSV/{sample}.vcf",
#         nanovar = IN_PATH + "/SVCall/NanoVar/{sample}.nanovar.pass.vcf",
#     output:
#         nanosv =  IN_PATH + "/SVCall/Format/{sample}_nanosv.vcf",
#         nanovar = IN_PATH + "/SVCall/Format/{sample}_nanovar.vcf",
#     params:
#         nanovar_change_format = SRC_DIR + "/nanovar_change_format.py",
#         nanosv_change_format = SRC_DIR + "/nanosv_change_format.py",
#     log:
#         IN_PATH + "/log/{sample}.ChangeFormat.log", 
#     run:
#         shell("python {params.nanosv_change_format} --raw_vcf {input.nanosv}  --new_vcf {output.nanosv} --score 20 > {log} 2>&1")
#         shell("python {params.nanovar_change_format} --raw_vcf {input.nanovar} --new_vcf {output.nanovar} >> {log} 2>&1")



# rule MultipleFilt:
#     input:
#         sniffles = IN_PATH + "/SVCall/Sniffles/{sample}.vcf",
#         nanosv =  IN_PATH + "/SVCall/Format/{sample}_nanosv.vcf",
#         nanovar = IN_PATH + "/SVCall/Format/{sample}_nanovar.vcf",
#         # quality = IN_PATH + "/QualityControl/Samples_quality_summary.xls",
#     output:
#         sniffles = IN_PATH + "/SVCall/SVFilt/{sample}_sniffles.vcf",
#         nanosv = IN_PATH + "/SVCall/SVFilt/{sample}_nanosv.vcf",
#         nanovar = IN_PATH + "/SVCall/SVFilt/{sample}_nanovar.vcf",
#     threads:
#         THREADS
#     params:
#         SVFiltReadsMultiple = SRC_DIR + "/SVFiltReadsMultiple.py",
#         ratioThreshold = 0.1,
#         column = "Clean_total_base",
#     log:
#         IN_PATH + "/log/{sample}.MultipleFilt.log", 
#     run:
#         ### TypeList = ["INS", "DEL", "INV", "DUP"]
#         ### CHRS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
#         shell("python {params.SVFiltReadsMultiple} --support 2 --vcf {input.sniffles} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.sniffles}  >{log} 2>&1")
#         shell("python {params.SVFiltReadsMultiple} --support 2 --vcf {input.nanosv} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.nanosv}   >>{log} 2>&1")
#         shell("python {params.SVFiltReadsMultiple} --support 2 --vcf {input.nanovar} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.nanovar}  >>{log} 2>&1")


# ###################################################################################################################################









########################################## filt based on SV reads number #########################





rule MultipleFilt2:
    input:
        vcf = IN_PATH + "/SVCall/Final/{sample}_common_original.vcf",
        # quality = IN_PATH + "/QualityControl/Samples_quality_summary.xls",
    output:
        vcf = IN_PATH + "/SVCall/Final/{sample}_common_original_read3.vcf",
    threads:
        THREADS
    params:
        SVFiltReadsMultiple = SRC_DIR + "/SVFiltReadsMultiple.py",
        ratioThreshold = 0.1,
        column = "Clean_total_base",
    log:
        IN_PATH + "/log/{sample}.MultipleFilt2.log", 
    run:
        ### support read = 3 
        shell("python {params.SVFiltReadsMultiple}  --support 3  --vcf {input.vcf} --quality . --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.vcf} >{log} 2>&1")


# rule SVFiltTypes:
#     input:
#         vcf = IN_PATH + "/SVCall/Final/{sample}_common_original_read3.vcf",
#     output:
#         vcf = IN_PATH + "/SVCall/SVFilt/{sample}_filt_types.vcf",
#     params:
#         SVFiltTypeChr = SRC_DIR + "/SVFiltTypeChr.py",
#         types = "DEL,INS,DUP,INV",
#         chrs = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{sample}.SVFiltTypes.log", 
#     run:
#         shell("python {params.SVFiltTypeChr} --vcf {input.vcf} --out {output.vcf} --types {params.types} --chrs {params.chrs} > {log} 2>&1")



# rule SVFiltReads:
#     input:
#         vcf = IN_PATH + "/SVCall/SVFilt/{sample}_filt_types.vcf",
#         quality = IN_PATH + "/Samples_quality_summary.xls",
#     output:
#         vcf = IN_PATH + "/SVCall/SVFilt/{sample}_filt_depth.vcf",
#     threads:
#         THREADS
#     params:
#         SVFiltReads = SRC_DIR + "/SVFiltReads.py",
#         ratioThreshold = 0.1,
#         column = "Clean_total_base",
#     log:
#         IN_PATH + "/log/{sample}.SVFiltReads.log", 
#     run:
#         # shell("python {params.SVFiltReads} --input {input.vcf} --out {output.vcf} --recordNumber {params.recordNumber} >{log} 2>&1")
#         shell("python {params.SVFiltReads} --vcf {input.vcf} --quality {input.quality} --ratioThreshold {params.ratioThreshold}  --column {params.column} --out {output.vcf} >{log} 2>&1")



rule FiltLength:
    input:
        vcf = IN_PATH + "/SVCall/Final/{sample}_common_original_read3.vcf",
    output:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_length.vcf",
        length = IN_PATH + "/SVCall/FiltStats/{sample}_length_length_stat.xls",
        number = IN_PATH + "/SVCall/FiltStats/{sample}_length_number_stat.xls",
    threads:
        THREADS
    params:
        SVFiltLength = SRC_DIR + "/SVFiltLength.py",
        delThreshold = 2000000,
        invThreshold = 20000000,
    log:
        IN_PATH + "/log/FiltLength_{sample}.log"
    run:
        ### 
        shell("python {params.SVFiltLength} --vcf {input.vcf} --out {output.vcf} --number {output.number} --length {output.length} --delThreshold {params.delThreshold}  --invThreshold {params.invThreshold} > {log} 2>&1")




rule FiltCentro:
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_length.vcf",
    output:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere-0.vcf",
    threads:
        THREADS
    params:
        SVFilt = SRC_DIR + "/SVFilt.py",
        centromere = config["centromere"],
        filtLength = 20000000,
    log:
        IN_PATH + "/log/FiltCentro_{sample}.log",
    run:
        shell("python {params.SVFilt} --vcf {input.vcf} --region {params.centromere} --out {output.vcf} --filtLength {params.filtLength} >{log} 2>&1")



rule DepthBed:
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere-0.vcf",
    output:
        bed_extend = IN_PATH + "/SVCall/Final/bed/{sample}_SV.bed",
    threads:
        THREADS
    params:
        SVvcf2bed = SRC_DIR + "/SVvcf2bed.py",
    log:
        IN_PATH + "/log/{sample}_SVBed.log",
    run:
        shell("python {params.SVvcf2bed} --vcf {input.vcf} --out {output.bed_extend} --distance 0 --select 'chrs' --method region 2>>{log}")



rule dpethOverlap:
    input:
        bed_extend = IN_PATH + "/SVCall/Final/bed/{sample}_SV.bed",
        bed_depth = IN_PATH + "/mapping/minimap2/bed/{sample}.regions.bed.gz",
    output:
        bed = IN_PATH + "/mapping/minimap2/bed/{sample}.regions.bed",
        bed_gap = IN_PATH + "/mapping/minimap2/bed/{sample}.regions_gap.bed",
        overlap = IN_PATH + "/mapping/minimap2/bed/{sample}_SV_overlap.bed",
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}_dpethOverlap.log",
    run:
        shell("gzip -dc {input.bed_depth} > {output.bed}")
        # cmd = "awk '{if ($4==0 || $4>100){print $0}}' %s > %s" % (output.bed, output.bed_gap)
        cmd = "awk '{if ($4>500){print $0}}' %s > %s" % (output.bed, output.bed_gap)
        print(cmd)
        os.system(cmd)
        shell("bedtools intersect -a {input.bed_extend} -b {output.bed_gap} -wa -wb > {output.overlap} 2>{log}")



rule excludeDepth:
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere-0.vcf",
        overlap = IN_PATH + "/mapping/minimap2/bed/{sample}_SV_overlap.bed",
    output:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf",
    threads:
        THREADS
    params:
        SVVCFExclude = SRC_DIR + "/SVVCFExclude.py",
    log:
        IN_PATH + "/log/{sample}_excludeDepth.log",
    run:
        shell("python {params.SVVCFExclude} --vcf {input.vcf} --bed {input.overlap} --out {output.vcf} --method exclude > {log} 2>&1")


#########################################################################################



################################# SV statistics summary ###################################
rule MultipleSummary:
    input:
        original = expand(IN_PATH + "/SVCall/Final/{sample}_common_original.vcf", sample=SAMPLES),
        types = expand(IN_PATH + "/SVCall/Final/{sample}_common_original_read3.vcf", sample=SAMPLES),
        length = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_length.vcf", sample=SAMPLES),
        centro = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=SAMPLES),
    output:
        summary = IN_PATH + "/SVCall/FiltStats/Samples_types_summary.xls",
    params:
        multipleFileSum = SRC_DIR + "/multipleFileSum.py",
    log:
        IN_PATH + "/log/MultipleSummary.log"
    run:
        VCFS = ",".join(input.original + input.types + input.length + input.centro)
        shell("python {params.multipleFileSum} --vcf {VCFS} --out {output.summary} > {log} 2>&1")



rule SummaryPlot:
    input:
        summary = IN_PATH + "/SVCall/FiltStats/Samples_types_summary.xls",
    output:
        pdf = IN_PATH + "/SVCall/FiltStats/Samples_SV_number_hist.pdf",
    params:
        SVDistHist = SCRIPT_DIR + "/SVDistHist.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/SummaryPlot.log"
    run:
        shell("Rscript {params.SVDistHist} --input {input.summary} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")


rule filtStep:
    input:
        summary = IN_PATH + "/SVCall/FiltStats/Samples_types_summary.xls",
    output:
        pdf = IN_PATH + "/SVCall/FiltStats/Samples_types_filt_steps.pdf",
    params:
        SVFilterBar = SCRIPT_DIR + "/SVFilterBar.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/filtStep.log"
    run:
        shell("Rscript {params.SVFilterBar} --input {input.summary} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")

#########################################################################################






############################## SV length summary #####################################
rule SVLengthSample:
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf",
    output:
        length = IN_PATH + "/SVCall/FiltStats/{sample}_SV_chr_length.xls",
    params:
        ChromSVLength = SRC_DIR + "/ChromSVLength.py",
        cutLength = 20000000,
    log:
        IN_PATH + "/log/SVLength_{sample}.log",
    run:
         shell("python {params.ChromSVLength} --vcf {input.vcf} --out {output.length} --cutLength {params.cutLength} --LenTag 'SVLEN' > {log} 2>&1")


rule SVLengthSum:
    input:
        vcf = expand(IN_PATH + "/SVCall/FiltStats/{sample}_SV_chr_length.xls", sample=SAMPLES),
    output:
        length = IN_PATH + "/SVCall/FiltStats/Samples_SV_length_summary.xls",
    log:
        IN_PATH + "/log/SVLengthSum.log",
    run:
        Files = sorted(input.vcf)
        for i in range(len(Files)):
            file = Files[i]
            if i == 0:
                cmd1 = "sed -n '1p' %s > %s" % (file, output.length)
                os.system(cmd1)
            cmd = "sed -n '$p' %s >> %s" % (file, output.length)
            os.system(cmd)


############################################################################################





############################################# IGV  #########################################

rule SVBed:
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf",
    output:
        bed_extend = IN_PATH + "/IGV/{sample}_SV.bed",
    threads:
        THREADS
    params:
        SVvcf2bed = SRC_DIR + "/SVvcf2bed.py",
        distance = config["IGV_distance"],
    log:
        IN_PATH + "/log/{sample}_SVBed.log",
    run:
        shell("python {params.SVvcf2bed} --vcf {input.vcf} --out {output.bed_extend} --distance {params.distance} --select 'chrs' 2>>{log}")


rule SVSnapshot:
    input:
        bed_extend = IN_PATH + "/IGV/{sample}_SV.bed",
        bam = "/home/wuzhikun/Project/NewChinese/mapping/minimap2/{sample}.bam",
    output:
        batch = IN_PATH + "/IGV/{sample}_SV.batch",
    threads:
        THREADS
    params:
        IGVsnapshot = SRC_DIR + "/IGVsnapshot.py",
        RefGenome = config["RefGenome"],
        outdir = IN_PATH + "/IGV/{sample}_SV",
        height = 500,
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/SVSnapshot_{sample}.log", 
    run:
        shell("python {params.IGVsnapshot} --bed {input.bed_extend} --out {output.batch} --outdir {params.outdir} --genome {params.RefGenome} --bam {input.bam} --sample {wildcards.sample} --height {params.height} >{log} 2>&1")


rule SVIGVScript:
    input:
        batch =  IN_PATH + "/IGV/{sample}_SV.batch",
    output:
        script = IN_PATH + "/IGV/{sample}_SV_run.sh",
    params:
        IGVBatchSplit = SRC_DIR + "/IGVBatchSplit.py",
        IGV = config["IGV"],
        IGVmemory = config["IGVmemory"],
        outdir = IN_PATH + "/IGV/IGV_batch/{sample}",
        number = 200,
    log:
        IN_PATH + "/log/{sample}_SVIGVScript.log",
    run:
        # shell("echo 'java {params.IGVmemory} -jar {params.IGV} -b  {input.batch}' > {output.script} 2>&1")
        # shell("echo 'python {params.IGVBatchSplit} --input {input.batch} --outdir {params.outdir} --number {params.number} >{log} 2>&1 ' ")
        shell("python {params.IGVBatchSplit} --input {input.batch} --outdir {params.outdir} --number {params.number} >{log} 2>&1")
        if os.path.exists(output.script):
            cmd = "rm %s" % output.script
            os.system(cmd)
        shell("""for i in $(ls {params.outdir}/*) ; do echo "java {params.IGVmemory} -jar {params.IGV} -b $i" >> {output.script} ; done""")
        shell("chmod +x {output.script}")

#############################################################################################
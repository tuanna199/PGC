

##################################### Quality control for nanopore #############################

rule NanoQCRaw:
    input:
        fastq = IN_PATH + "/raw/{sample}.fastq.gz",
        # fastq = IN_PATH + "/raw/{sample}.fastq",
    output:
        temp1 = temp(IN_PATH + "/raw/{sample}_sub.fastq"),
        qc = IN_PATH + "/QualityControl/raw/NanoQC/{sample}/nanoQC.html",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/raw/NanoQC/{sample}",
        sub_freq = config["sub_freq"],
    log:
        IN_PATH + "/log/nanoQC_raw_{sample}.log"
    run:
        shell("seqtk sample -s100 {input.fastq} {params.sub_freq}  > {output.temp1} 2>>{log}")
        shell("nanoQC --outdir {params.outdir} {output.temp1} 2>>{log}")


rule NanoPlotRaw:
    ### It contain the result of NanoStats
    input:
        fastq = IN_PATH + "/raw/{sample}.fastq.gz",
    output:
        stat = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}/NanoStats.txt",
        hist = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}/HistogramReadlength.pdf",
        qua = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}/LengthvsQualityScatterPlot_hex.pdf",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}",
    log:
        IN_PATH + "/log/NanoPlot_raw_{sample}.log"
    run:
        # shell("NanoPlot  --outdir {params.outdir}  --fastq {input.fastq} --loglength --format pdf --plots hex dot >{log} 2>&1")
        cmd = "source activate nanovar && NanoPlot  --outdir %s  --fastq %s --loglength --format pdf --plots hex dot > %s 2>&1" % (params.outdir, input.fastq, log)
        print(cmd)
        os.system(cmd)



rule NanoFilt:
    input:
        fastq = IN_PATH + "/raw/{sample}.fastq.gz",
        # fastq = IN_PATH + "/raw/{sample}.fastq",
    output:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    threads:
        THREADS
    log:
        IN_PATH + "/log/NanoFilt_{sample}.log"
    params:
        readtype = config["readtype"],
        minQuality = config["minQuality"],
        minLength = config["minLength"],
        headcrop = config["headcrop"],
        tailcrop = config["tailcrop"],
        # phage_lambda = config["phage_lambda"],
    run:
        # ### Filt out reads of phage lambda
        # shell("gunzip -c {input.fastq} | NanoLyse --reference {params.phage_lambda} | NanoFilt --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} | gzip > {output.fastq} 2>{log}")
        # shell("gunzip -c {input.fastq} | NanoFilt --readtype {params.readtype} --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} | gzip > {output.fastq} 2>{log}")
        cmd = "source activate nanovar &&  gunzip -c %s | NanoFilt --readtype %s --quality  %s --length %s --headcrop %s --tailcrop %s | gzip > %s 2>%s" % (input.fastq, params.readtype, params.minQuality, params.minLength, params.headcrop, params.tailcrop, output.fastq, log)
        print(cmd)
        os.system(cmd)
        # shell("NanoFilt --readtype {params.readtype} --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} < {input.fastq} | gzip > {output.fastq} 2>{log}")



rule md5:
    input:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        md5 = IN_PATH + "/QualityControl/md5/{sample}.md5.txt",
    run:
        shell("md5sum {input.fastq} > {output.md5}")



rule NanoQCClean:
    input:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        temp1 = temp(IN_PATH + "/clean/{sample}_sub.fastq"),
        qc = IN_PATH + "/QualityControl/clean/NanoQC/{sample}/nanoQC.html",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/clean/NanoQC/{sample}",
        sub_freq = config["sub_freq"],
    log:
        IN_PATH + "/log/nanoQC_clean_{sample}.log"
    run:
        shell("seqtk sample -s100 {input.fastq} {params.sub_freq}  > {output.temp1} 2>>{log}")
        shell("nanoQC --outdir {params.outdir} {output.temp1} 2>>{log}")



rule NanoPlotClean:
    ### It contain the result of NanoStats
    input:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        stat = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}/NanoStats.txt",
        hist = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}/HistogramReadlength.pdf",
        qua = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}/LengthvsQualityScatterPlot_hex.pdf",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}",
    log:
        IN_PATH + "/log/NanoPlot_clean_{sample}.log"
    run:
        # shell("NanoPlot  --outdir {params.outdir}  --fastq {input.fastq} --loglength --format pdf --plots hex dot >{log} 2>&1")
        cmd = "source activate nanovar &&  NanoPlot  --outdir %s  --fastq %s --loglength --format pdf --plots hex dot > %s 2>&1" % (params.outdir, input.fastq, log)
        print(cmd)
        os.system(cmd)


rule QualityStats:
    input:
        raw = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}/NanoStats.txt",
        clean = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}/NanoStats.txt",
    output:
        stats = IN_PATH + "/QualityControl/{sample}_stats.xls",
    threads:
        THREADS
    params:
        NanoStatSummary = SRC_DIR + "/NanoStatSummary.py",
    log:
        IN_PATH + "/log/QualityStats_{sample}.log"
    run:
        shell("python {params.NanoStatSummary} --stats {input.raw} --clean {input.clean} --sample {wildcards.sample} --out {output.stats} >{log} 2>&1")



rule pdf2png:
    input:
        raw = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}/LengthvsQualityScatterPlot_hex.pdf",
        clean = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}/LengthvsQualityScatterPlot_hex.pdf",
    output:
        raw = IN_PATH + "/QualityControl/raw/NanoPlot/{sample}/LengthvsQualityScatterPlot_hex.png",
        clean = IN_PATH + "/QualityControl/clean/NanoPlot/{sample}/LengthvsQualityScatterPlot_hex.png",
    threads:
        THREADS
    params:
        pdf2png = SRC_DIR + "/pdf2png.py",
    run:
        shell("python {params.pdf2png} --input {input.raw} --threads 1")
        shell("python {params.pdf2png} --input {input.clean} --threads 1")



rule QualityStatsMerge:
    input:
        stats = expand(IN_PATH + "/QualityControl/{sample}_stats.xls", sample=SAMPLES),
    output:
        stats = IN_PATH + "/QualityControl/Samples_quality_summary.xls",
    run:
        files = input.stats
        for i in range(len(files)):
            if i == 0:
                cmd = "cat %s > %s" % (files[i], output.stats)
            else:
                cmd = "sed '1d' %s >> %s" % (files[i], output.stats)
            os.system(cmd)


rule ReadBaseHist:
    input:
        stats = IN_PATH + "/QualityControl/Samples_quality_summary.xls",
    output:
        hist = IN_PATH + "/QualityControl/Samples_quality_summary_hist.pdf",
    params:
        ReadBaseHist = SCRIPT_DIR + "/ReadBaseHist.R",
        width = 8,
        height = 4,
    log:
        IN_PATH + "/log/ReadBaseHist.log"
    run:
        # shell("source activate Rmeta && Rscript {params.ReadBaseHist} --input {input.stats} --pdf {output.hist}  --width {params.width} --height {params.height} >{log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s  --width %s --height %s > %s 2>&1" % (params.ReadBaseHist, input.stats, output.hist, params.width, params.height, log)
        os.system(cmd)




rule QualityHist:
    input:
        stats = IN_PATH + "/QualityControl/Samples_quality_summary.xls",
    output:
        hist1 = IN_PATH + "/QualityControl/Samples_read_length_N50_hist.pdf",
        hist2 = IN_PATH + "/QualityControl/Samples_read_quality_hist.pdf",
        hist3 = IN_PATH + "/QualityControl/Samples_read_depth_hist.pdf",
    params:
        QualityHist = SCRIPT_DIR + "/QualityHist.R",
        outPrefix = IN_PATH + "/QualityControl/Samples_read",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/QualityHist.log"
    run:
        shell("Rscript {params.QualityHist} --input {input.stats} --outPrefix {params.outPrefix} --width {params.width} --height {params.height} > {log} 2>&1")
#################################################################################################









################################# Error rate estimate #########################################

    
rule ErrorRateSummary:
    input:
        rate = expand(IN_PATH + "/mapping/minimap2/{sample}_error_rate.xls", sample=SAMPLES),
    output:
        rate = IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.xls",
    run:
        merge_multiple_sample_stats(input.rate, output.rate)


rule ErrorRatePlot:
    input:
        rate = IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.xls",
    output:
        pdf = IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.pdf",
        pdf2 = IN_PATH + "/mapping/minimap2/Samples_mapping_rate_stats.pdf",
    params:
        ErrorRateBar = SCRIPT_DIR + "/ErrorRateBar.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/ErrorRatePlot.log"
    run:
        shell("Rscript {params.ErrorRateBar} --input {input.rate} --pdf {output.pdf} --pdf2 {output.pdf2} --width {params.width} --height {params.height} > {log} 2>&1")

##############################################################################################



##################################### Assembly statistics  #####################################


# rule assemblyStats:
#     input:
#         stats = IN_PATH + "/Assembly/Assembly_stats.txt",
#     output:
#         pdf = IN_PATH + "/Assembly/Assembly_stats.pdf",
#         pdf2 = IN_PATH + "/Assembly/Assembly_N50_stats.pdf",
#     params:
#         AssemblyHist = SCRIPT_DIR + "/AssemblyHist.R",
#         AssemblyN50Hist = SCRIPT_DIR + "/AssemblyN50Hist.R",
#         width = 4,
#         height = 4,
#     log:
#         IN_PATH + "/log/assemblyStats.log"
#     run:
#         shell("Rscript {params.AssemblyHist} --input {input.stats} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
#         shell("Rscript {params.AssemblyN50Hist} --input {input.stats} --pdf {output.pdf2} --width {params.width} --height {params.height} >> {log} 2>&1")

################################################################################################





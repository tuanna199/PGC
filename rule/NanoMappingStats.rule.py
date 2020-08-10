
################################## mapping using minimap2 ######################
rule minimap2Align:
    input:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        sam = temp(IN_PATH + "/mapping/minimap2/{sample}.sam"),
    threads:
        THREADS * ThreadFold
    params:
        RefGenome = config["RefGenome"],
    log:
        IN_PATH + "/log/minimap2Align_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {params.RefGenome} {input.fastq} > {output.sam} 2>{log}")



rule SAM2BAM:
    input:
        sam = IN_PATH + "/mapping/minimap2/{sample}.sam",
    output:
        tempbam = temp(IN_PATH + "/mapping/minimap2/{sample}_temp.bam"),
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
        bai = IN_PATH + "/mapping/minimap2/{sample}.bam.bai",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/sam2bam_minimap2_{sample}.log"
    run:
        # shell('sambamba view --nthreads {threads} --sam-input --format bam  -o {output.tempbam} {input.sam} >{log} 2>&1' )
        # shell('sambamba sort --nthreads {threads} -o {output.bam} {output.tempbam} 2>>{log}')
        # shell("sambamba index --nthreads {threads} {output.bam} 2>>{log}")
        shell('samtools view -Sb --threads {threads}  -o {output.tempbam} {input.sam}   >{log} 2>&1')
        shell('samtools sort --threads {threads} -o {output.bam} {output.tempbam} 2>>{log}')
        shell("samtools index -@ {threads} {output.bam} 2>>{log}")
##############################################################################################




################################## BAM statistics ##########################################
rule BAMFlag:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        flag = IN_PATH + "/mapping/minimap2/{sample}_bam_flag.txt",
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMFlag_{sample}.log",
    run:
        shell(" samtools flagstat -@ {threads} {input.bam} > {output.flag} 2>{log}")
    


rule BAMFlagSum:
    input:
        flag = expand(IN_PATH + "/mapping/minimap2/{sample}_bam_flag.txt", sample=SAMPLES),
    output:
        summary = IN_PATH + "/mapping/minimap2/Samples_mapping_flag_summary.xls",
    params:
        bamFlagSum = SRC_DIR + "/bamFlagSum.py",
    log:
        IN_PATH + "/log/BAMFlagSum.log",
    run:
        FILES = ",".join(input.flag)
        shell("python {params.bamFlagSum} --file {FILES} --out {output.summary} > {log} 2>&1")






rule BAMStats:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        stat = IN_PATH + "/mapping/minimap2/{sample}_bam_stats.txt",
        summary = IN_PATH + "/mapping/minimap2/{sample}_bam_summary.xls",
    params:
        NanoBamStat = SRC_DIR + "/NanoBamStat.py",
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")
        shell("python {params.NanoBamStat} --stat {output.stat} --out {output.summary} --sample {wildcards.sample} 2>>{log}")


def read_file(file):
    in_h = open(file, "r")
    header = in_h.readline().strip()
    records = in_h.readline().strip().split("\t")
    record = "\t".join(records[1:])
    in_h.close()
    return header, record



def merge_bam_stats(bam_files, out_file):
    FileRecord = {}
    files = bam_files.split(",")
    for f in files:
        file_base = f.split("/")[-1].split("_")[0]
        header, record = read_file(f)
        FileRecord[file_base] = record
    out_h = open(out_file, "w")
    out_h.write("%s\n" % header)
    for ff in sorted(list(FileRecord.keys())):
        rec = FileRecord[ff]
        out_h.write("%s\t%s\n" % (ff, rec))
    out_h.close()
    





rule BAMStatsMerge:
    input:
        summary = expand(IN_PATH + "/mapping/minimap2/{sample}_bam_summary.xls", sample=SAMPLES),
    output:
        sumAll = IN_PATH + "/mapping/Samples_mapping_summary.xls",
    run:
        # bamSums = input.summary
        # merge_multiple_sample_stats(bamSums, output.sumAll)
        bamSums = ",".join(input.summary)
        merge_bam_stats(bamSums, output.sumAll)



rule mappingRateStats:
    input:
        summary = IN_PATH + "/mapping/Samples_mapping_summary.xls",
    output:
        pdf = IN_PATH + "/mapping/Samples_mapping_summary.pdf",
    params:
        mappingRateHist = SCRIPT_DIR + "/mappingRateHist.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/mappingRateStats.log"
    run:
        shell("Rscript {params.mappingRateHist} --input {input.summary} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")


##############################################################################################



################################# Error rate estimate #########################################
rule ErrorRate:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        rate = IN_PATH + "/mapping/minimap2/{sample}_error_rate.xls",
    threads:
        THREADS
    params:
        pysam_qc = SRC_DIR + "/pysam_qc.py",
    log:
        IN_PATH + "/log/ErrorRate_{sample}.log"
    run:
        shell("python {params.pysam_qc} --bam {input.bam} --out {output.rate} > {log} 2>&1")



# rule ErrorRateSummary:
#     input:
#         rate = expand(IN_PATH + "/mapping/minimap2/{sample}_error_rate.xls", sample=SAMPLES),
#     output:
#         rate = IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.xls",
#     run:
#         merge_multiple_sample_stats(input.rate, output.rate)


# rule ErrorRatePlot:
#     input:
#         rate = IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.xls",
#     output:
#         pdf = IN_PATH + "/mapping/minimap2/Samples_error_rate_stats.pdf",
#     params:
#         ErrorRateBar = SCRIPT_DIR + "/ErrorRateBar.R",
#         width = 4,
#         height = 4,
#     log:
#         IN_PATH + "/log/ErrorRatePlot.log"
#     run:
#         shell("Rscript {params.ErrorRateBar} --input {input.rate} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")

##############################################################################################


################################## summary of depth ##########################################

rule mosdepth:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        depth = IN_PATH + '/mapping/minimap2/{sample}/{sample}.regions.bed.gz',
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + '/mapping/minimap2/{sample}/{sample}',
    log:
        IN_PATH + "/log/mosdepth_{sample}.log"
    run:
        shell('mosdepth --threads {threads} --no-per-base --fast-mode  --by 5000  --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')




rule DepthChr:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        depth = temp(IN_PATH + '/mapping/minimap2/{sample}/{sample}_realigned.depth_{chr}.txt'),
    threads:
        THREADS
    run:
        shell('samtools depth -r {wildcards.chr} {input.bam} > {output.depth}')

rule DepthCoverageChr:
    input:
        depth = IN_PATH + '/mapping/minimap2/{sample}/{sample}_realigned.depth_{chr}.txt',
    output:
        temp = temp(IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage_{chr}_temp.txt'),
    threads:
        THREADS
    params:
        DepthCovStats = SRC_DIR + '/DepthCovStats.py',
        RefGenome = config['RefGenome'],
    run:
        shell('python {params.DepthCovStats} --input {input.depth} --fasta {params.RefGenome} --out {output.temp}')


rule DepthCoverageStats:
    input:
        depth = expand(IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage_{chr}_temp.txt', sample=SAMPLES, chr=CHRS),
    output:
        temp = IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage_temp.xls',
        coverage = IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage.xls',
    threads:
        THREADS
    params:
        DepthCovStats = SRC_DIR + '/DepthCovStats.py',
        RefGenome = config['RefGenome'],
        outdir = IN_PATH + '/mapping/minimap2/{sample}/',
    log:
        IN_PATH + "/log/{sample}_DepthCoverageStats.log",
    run:
        allFiles = input.depth
        files = []
        for f in allFiles:
            if params.outdir in f:
                files.append(f)
        chrLen = len(CHRS)
        merge_multiple_sample_stats(files, output.temp)
        shell("sed -n '1p' {output.temp} > {output.coverage}")
        shell("sort -k 3nr {output.temp}  | sed -n '1,{chrLen}p' | sort -k 1n >> {output.coverage}")
        # shell('python {params.DepthCovStats} --input {output.depth} --fasta {params.RefGenome} --out {output.temp} >{log} 2>&1')
        # shell("grep -v '^G' {output.temp} |  grep -v '^K' | grep -v '^M' > {output.out}")



rule DepthCoveragePlot:
    input:
        depth = IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage.xls',
    output:
        sortedChr = IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage_sortedChr.xls',
        ### the name is fixed and ends with ".Depth.Coverage.pdf"
        pdf = IN_PATH + '/mapping/minimap2/{sample}/{sample}.Depth.Coverage.pdf',
    params:
        SortChrDigitPart = SRC_DIR + '/SortChrDigitPart.py',
        DepthCoverPlot = SRC_DIR + '/DepthCoverPlot.py',
        odir = IN_PATH + '/mapping/minimap2/{sample}',
        maxDepth = config['maxDepth'],
    log:
        IN_PATH + "/log/{sample}.DepthCoveragePlot.log",    
    run:
        shell('python {params.SortChrDigitPart} --input {input.depth} --out {output.sortedChr} >{log} 2>&1')
        shell('python {params.DepthCoverPlot} --sample {wildcards.sample} --data {input.depth} --odir {params.odir} --depth {params.maxDepth} >{log} 2>&1')



rule DepthCoverageSum:
    input:
        depth = expand(IN_PATH + '/mapping/minimap2/{sample}/{sample}_depth_coverage_sortedChr.xls', sample=SAMPLES),
    output:
        depth = IN_PATH + '/mapping/minimap2/Samples_depth_coverage_average.xls',
    params:
        averageCoverage = SRC_DIR + "/averageCoverage.py",
    log:
        IN_PATH + "/log/DepthCoverageSum.log",    
    run:
        DEPTH = ",".join(input.depth)
        shell("python {params.averageCoverage} --files {DEPTH} --out {output.depth} > {log} 2>&1")



rule DepthCoverageAveragePlot:
    input:
        depth = IN_PATH + '/mapping/minimap2/Samples_depth_coverage_average.xls',
    output:
        pdf = IN_PATH + '/mapping/minimap2/Samples_depth_average.pdf',
        pdf2 = IN_PATH + '/mapping/minimap2/Samples_coverage_average.pdf',
    params:
        AverageDepth = SCRIPT_DIR + "/AverageDepth.R",
        width = 7,
        height = 4,
    log:
        IN_PATH + "/log/DepthCoveragePlot.log",
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --pdf2 %s --width %s --height %s > %s 2>&1" % (params.AverageDepth, input.depth, output.pdf, output.pdf2, params.width, params.height, log)
        os.system(cmd)


rule ReadDepthWinChr:
    input:
        depth = IN_PATH + '/mapping/minimap2/{sample}/{sample}_realigned.depth_{chr}.txt',
    output:
        depth = temp(IN_PATH + '/mapping/minimap2/{sample}/{sample}_region_depth_{chr}_win.txt'),
    params:
        SortReadDepthWin = SRC_DIR + '/SortReadDepthWin.py',
        windowSize = config['windowSize'],
    log:
        IN_PATH + "/log/{sample}_{chr}.ReadDepthWin.log",  
    run:
        #sambamba depth window  --window-size 100   --nthreads 10 -o M625-0_depth.txt  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M625-0.bam
        shell('python {params.SortReadDepthWin} --depth {input.depth} --out {output.depth} --window {params.windowSize} >{log} 2>&1')


rule ReadDepthWin:
    input:
        depth = expand(IN_PATH + '/mapping/minimap2/{sample}/{sample}_region_depth_{chr}_win.txt', chr=CHRS, sample=SAMPLES),
    output:
        depth = IN_PATH + '/mapping/minimap2/{sample}/{sample}_region_depth_window.txt',
    params:
        outdir = IN_PATH + '/mapping/minimap2/{sample}/',
    log:
        IN_PATH + "/log/{sample}.ReadDepthWin.log",  
    run:
        #sambamba depth window  --window-size 100   --nthreads 10 -o M625-0_depth.txt  /home/wuzhikun/Project/NanoTrio/mapping/minimap2/M625-0.bam
        allFiles = input.depth
        files = []
        for f in allFiles:
            if params.outdir in f:
                files.append(f)
        merge_multiple_sample_stats(files, output.depth)



rule ReadDepthWinPlot:
    input:
        depth = rules.ReadDepthWin.output.depth,
    output:
        pdf = IN_PATH + '/mapping/minimap2/{sample}/{sample}_region_depth.pdf',
    params:
        ReadDepthPlot = SCRIPT_DIR + '/ReadDepthPlot.R',
        windowSize = config['windowSize'],
        chromosomes = config['chromosome'],
        ylim = config['maxDepth'],
        width = config['width'],
        height = config['height'],
    log:
        IN_PATH + "/log/{sample}.ReadDepthWinPlot.log",
    run:
        # shell('source activate Rmeta && Rscript {params.ReadDepthPlot} --depth {input.depth} --window_size {params.windowSize} --chromosomes {params.chromosomes} --ylim {params.ylim} --width {params.width} --height {params.height} --pdf {output.pdf} > {log} 2>&1')
        cmd = "source activate Rmeta && Rscript %s --depth %s --window_size %s --chromosomes %s --ylim %s --width %s --height %s --pdf %s > %s 2>&1" % (params.ReadDepthPlot, input.depth, params.windowSize, params.chromosomes, params.ylim, params.width, params.height, output.pdf, log)
        os.system(cmd)

##############################################################################################


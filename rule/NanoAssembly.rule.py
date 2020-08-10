
################################ genome assembly ###########################
rule genomeAssembly:
    input:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        lay = IN_PATH + "/Assembly/{sample}/{sample}_assembly.ctg.lay.gz",
    threads:
        THREADS * ThreadFold
    params:
        outprefix = IN_PATH + "/Assembly/{sample}/{sample}_assembly",
    log:
        IN_PATH + "/log/genomeAssembly_{sample}.log"
    run:
        ### -S 3 can reduce the memory
        ### -g 3g
        cmd1 = "source activate Assembly && wtdbg2  -i %s  -fo %s  -t %s -p 19 -AS 2 -s 0.05 -L 1000 -S 4 > %s 2>&1" % (input.fastq, params.outprefix, threads, log)
        print(cmd1)
        os.system(cmd1)



rule genomeAssembly2:
    input:
        lay = IN_PATH + "/Assembly/{sample}/{sample}_assembly.ctg.lay.gz",
    output:
        fasta = IN_PATH + "/Assembly/{sample}/{sample}_assembly.fasta",
    threads:
        THREADS * ThreadFold
    params:
        outprefix = IN_PATH + "/Assembly/{sample}/{sample}_assembly",
    log:
        IN_PATH + "/log/genomeAssembly2_{sample}.log"
    run:
        cmd2 = "source activate Assembly && wtpoa-cns -t %s  -i  %s -fo %s > %s 2>&1" % (threads, input.lay, output.fasta, log)
        print(cmd2)
        os.system(cmd2)


rule assemblyMapping:
    input:
        fasta = IN_PATH + "/Assembly/{sample}/{sample}_assembly.fasta",
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        bam = IN_PATH + "/Assembly/mapping/{sample}_assembly.bam",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/assemblyMapping_{sample}.log"
    run:
        shell("minimap2 -t {threads} -ax map-ont -r2k {input.fasta} {input.fastq} | samtools sort -@ {threads} > {output.bam} 2>{log}")




rule assemblyPolish:
    input:
        bam = IN_PATH + "/Assembly/mapping/{sample}_assembly.bam",
        fasta = IN_PATH + "/Assembly/{sample}/{sample}_assembly.fasta",
    output:
        fasta = IN_PATH + "/Assembly/{sample}/{sample}_assembly_polish.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/assemblyPolish_{sample}.log"
    run:
        # shell("samtools view -F0x900 {input.bam} | wtpoa-cns -t {threads} -d {input.fasta} -i -  -fo {output.fasta} > {log} 2>&1")
        cmd = "source activate Assembly && samtools view --threads %s  -F0x900 %s | wtpoa-cns -t %s -d %s -i -  -fo %s > %s 2>&1" % (threads, input.bam, threads, input.fasta, output.fasta, log)
        os.system(cmd)



rule quastEvaluate:
    input:
        fasta = IN_PATH + "/Assembly/{sample}/{sample}_assembly_polish.fasta",
    output:
        # report = IN_PATH + '/Assembly/Quast/{sample}/polish/report.html',
        txt = IN_PATH + '/Assembly/Quast/{sample}/polish/report.txt',
    threads:
        THREADS
    params:
        REF = config["RefGenome"],
        GFF = config["GFF"],
        outdir = IN_PATH + '/Assembly/Quast/{sample}/polish',
        min_contig = 1000,
    log:
        IN_PATH + '/log/{sample}_quastpolish.log'         
    run:
        cmd = "source activate Assembly && quast.py --no-html --no-snps  -o %s -r %s -g %s -m %s -t %s %s > %s 2>&1" % (params.outdir, params.REF, params.GFF, params.min_contig, threads, input.fasta, log)
        os.system(cmd)

####################################################################################################



##################################### Assembly statistics  #####################################


rule assemblyStats:
    input:
        stats = IN_PATH + "/Assembly/Assembly_stats.txt",
    output:
        pdf = IN_PATH + "/Assembly/Assembly_stats.pdf",
        pdf2 = IN_PATH + "/Assembly/Assembly_N50_stats.pdf",
    params:
        AssemblyHist = SCRIPT_DIR + "/AssemblyHist.R",
        AssemblyN50Hist = SCRIPT_DIR + "/AssemblyN50Hist.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/assemblyStats.log"
    run:
        shell("Rscript {params.AssemblyHist} --input {input.stats} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
        shell("Rscript {params.AssemblyN50Hist} --input {input.stats} --pdf {output.pdf2} --width {params.width} --height {params.height} >> {log} 2>&1")

################################################################################################


##################################### call SV using Sniffles (minimap2) ###################################
rule Sniffles:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        vcf = IN_PATH + "/SVCall/Sniffles/minimap2/{sample}.vcf",
    threads:
        THREADS * ThreadFold
    params:
        min_support = config["min_support"],
        min_length  = config["min_length"],
        minmapping_qual = config["minmapping_qual"],
        num_reads_report = config["num_reads_report"],
        min_seq_size = config["min_seq_size"],
    log:
        IN_PATH + "/log/SnifflesNano_minimap2_{sample}.log"
    run:
        shell("sniffles --mapped_reads {input.bam} --vcf {output.vcf} --threads {threads}  --min_support {params.min_support} --min_length {params.min_length} --minmapping_qual {params.minmapping_qual} --num_reads_report {params.num_reads_report} --min_seq_size {params.min_seq_size}  --genotype --report_BND --report_seq  >{log} 2>&1")


##################################################################################################





################################### call SV using NanoSV (minimap2) ###############################
rule NanoSV:
    input:
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        vcf = IN_PATH + "/SVCall/NanoSV/minimap2/{sample}.vcf",
    threads:
        THREADS * ThreadFold
    params:
        nanosvConfig = config["nanosvConfig"],
    log:
        IN_PATH + "/log/NanoSV_minimap2_{sample}.log"
    run:
        ### depth_support = False
        shell("NanoSV --threads {threads} -c {params.nanosvConfig} -o {output.vcf} {input.bam} >{log} 2>&1")


rule NanoSVConvert:
    input:
        vcf = IN_PATH + "/SVCall/NanoSV/minimap2/{sample}.vcf",
    output:
        vcf = IN_PATH + "/SVCall/NanoSV/minimap2/{sample}_type.vcf",
    params:
        NanosvType = SRC_DIR + "/NanosvType.py",
    log:
        IN_PATH + "/log/NanoSVConvert_{sample}.log"
    run:
        shell("python {params.NanosvType} --vcf {input.vcf} --out {output.vcf} >{log} 2>&1")
############################################################################################################




################################### call SV using nanovar ######################################
rule nanovar:
    input:
        # fastq = IN_PATH + "/raw/{sample}.fastq.gz",
        bam = IN_PATH + "/mapping/minimap2/{sample}.bam",
    output:
        total = IN_PATH + "/SVCall/nanovar/{sample}/{sample}.nanovar.total.vcf",
        vcf = IN_PATH + "/SVCall/nanovar/{sample}/{sample}.nanovar.pass.vcf",
    threads:
        THREADS
    params:
        RefGenome = config["RefGenome"],
        outdir = IN_PATH + "/SVCall/nanovar/{sample}",
    log:
        IN_PATH + "/log/nanovar_{sample}.log"
    run:
        # shell("nanovar -r {params.RefGenome} -l {input.fastq} -t {threads} -o {params.outdir} >{log} 2>&1")
        ### nanovar was installed in another conda environment
        cmd = "source activate nanovar && nanovar -t %s --data_type ont --mincov 2 --minlen 50 %s %s %s > %s 2>&1" % (threads, input.bam, params.RefGenome, params.outdir, log)
        print(cmd)
        os.system(cmd)

################################################################################################



############################  HG002 INV ##################
rule HG002SV:
    input:
        sv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002_filt_centromere.vcf",
    output:
        Del = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_DEL.bed",
        Ins = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INS.bed",
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV.bed",
        length = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_SV_length.txt",
    params:
        outPrefix = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT",
        SVVcf2MultBed = SRC_DIR + "/SVVcf2MultBed.py",
    log:
        IN_PATH + "/log/HG002SV.log",
    run:
        shell("python {params.SVVcf2MultBed} --vcf {input.sv} --out {params.outPrefix} --lengthOut {output.length} > {log} 2>&1")



rule HG002HiFi:
    input:
        sv = "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002_filt_centromere.vcf",
    output:
        Del = "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002/HG002_ONT_DEL.bed",
        Ins = "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002/HG002_ONT_INS.bed",
        Inv = "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002/HG002_ONT_INV.bed",
        length = "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002/HG002_ONT_SV_length.txt",
    params:
        outPrefix = "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002/HG002_ONT",
        SVVcf2MultBed = SRC_DIR + "/SVVcf2MultBed.py",
    log:
        IN_PATH + "/log/HG002HiFi.log",
    run:
        shell("python {params.SVVcf2MultBed} --vcf {input.sv} --out {params.outPrefix} --lengthOut {output.length} > {log} 2>&1")






rule HG002Overlap:
    input:
        HG002_Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV.bed",
        HG002_HiFi =  "/home/wuzhikun/Project/HiFi/SVCall/FinalFilt/HG002/HG002_ONT_INV.bed",
    output:
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_INV_overlap.bed",
    run:
        shell("bedtools intersect -a {input.HG002_Inv} -b {input.HG002_HiFi}  -wa -wb > {output.Inv}")





rule HG002OverlapFilt:
    input:
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_INV_overlap.bed",
    output:
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_INV_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/LongOverlapFilt.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Inv} --ratioThreshold {params.ratioThreshold} --out {output.Inv} --method both 2>>{log}")








rule HG002HGSVC:
    input:
        HG002_Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV.bed",
        HGSVC = "/home/wuzhikun/database/HGSVC2/Science/HG002_freeze4inv_sv_inv-1.tsv",
    output:
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HGSVC_INV_overlap.bed",
    run:
        shell("bedtools intersect -a {input.HG002_Inv} -b {input.HGSVC}  -wa -wb > {output.Inv}")

    

rule HG002HGSVCFilt:
    input:
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HGSVC_INV_overlap.bed",
    output:
        Inv = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HGSVC_INV_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/LongOverlapFilt.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Inv} --ratioThreshold {params.ratioThreshold} --out {output.Inv} --method both 2>>{log}")





rule OverlapINV:
    input:
        HiFi = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_INV_overlap_filt.bed",
        HGSVC = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HGSVC_INV_overlap_filt.bed",
    output:
        merge = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_HGSVC_overlapped_INV.bed",
    run:
        shell("cat {input.HiFi} {input.HGSVC} | cut -f 1-4 | sort | uniq > {output.merge}")



rule unoverlap:
    input:
        HG002 = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV.bed",
        merge = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_HGSVC_overlapped_INV.bed",
    output:
        unoverlap = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV_unoverlaped.bed",
    run:
        shell("cat {input.HG002} {input.merge} | sort | uniq -c | sort -k 1nr | sed 's/^      //g' | grep '^1' | cut -f 2 -d ' ' | sort -k 1,1n -k 2,2n > {output.unoverlap}")






def extend_bed(bed_file, out_file, extend=0.1):
    extend = float(extend)
    bed_h = open(bed_file, "r")
    out_h = open(out_file, "w")
    for line in bed_h:
        lines = line.strip().split("\t")
        chrom, start, end, tag = lines[:4]
        start = int(start)
        end = int(end)
        length = end - start
        distance = int(extend * length)
        newStart = start - distance
        newEnd = end + distance
        out_h.write("%s\t%d\t%d\t%s\n" % (chrom, newStart, newEnd, tag))
    bed_h.close()
    out_h.close()
        

rule IGVsnapshot:
    input:
        bed = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV_unoverlaped.bed",
        bam1 = "/home/wuzhikun/Project/Revised/mapping/minimap2/HG002.bam",
        bam2 = "/home/wuzhikun/Project/HiFi/mapping/minimap2/HG002.bam",
    output:
        bed = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV_unoverlaped_extend.bed",
        batch = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_igv.batch",
    params:
        IGVsnapshot = SRC_DIR + "/IGVsnapshot.py",
        RefGenome = config["RefGenome"],
        outdir = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_HiFi_igv",
        height = 500,
    log:
        IN_PATH + "/log/IGVsnapshot.log",
    run:
        bamFile = ",".join([input.bam1, input.bam2])
        extend_bed(input.bed, output.bed, extend=0.1)
        shell("python {params.IGVsnapshot} --bed {output.bed} --out {output.batch} --outdir {params.outdir} --genome {params.RefGenome} --bam {bamFile} --sample HG002 --height {params.height} >{log} 2>&1")



rule extractBAM:
    input:
        bed = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV_unoverlaped.bed",
        bam1 = "/home/wuzhikun/Project/Revised/mapping/minimap2/HG002.bam",
        bam2 = "/home/wuzhikun/Project/HiFi/mapping/minimap2/HG002.bam",
    output:
        bed = "/home/wuzhikun/Project/Revised/SVCall/FinalFilt/HG002/HG002_ONT_INV_unoverlaped_extend1.bed",
        bam1 = "/home/wuzhikun/Project/Revised/mapping/minimap2/HG002_INV.bam",
        bam2 = "/home/wuzhikun/Project/HiFi/mapping/minimap2/HG002_INV.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/extractBAM.log",
    run:
        extend_bed(input.bed, output.bed, extend=1)
        shell("samtools view -@ {threads} -b -L {output.bed} -o {output.bam1} {input.bam1} >{log} 2>&1")
        shell("samtools index {output.bam1}")
        shell("samtools view -@ {threads} -b -L {output.bed} -o {output.bam2} {input.bam2} 2>>{log}")
        shell("samtools index {output.bam2}")
        

##########################################################





################## IBS ################################
rule plinkIBS:
    input:
        bed = IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
    output:
        genome = IN_PATH + "/population/subgroup/plink.genome",
    params:
        inPrefix = IN_PATH + "/population/subgroup/Sample_SV_geno",
        outDir = IN_PATH + "/population/subgroup",
    log:
        IN_PATH + "/log/plinkIBS.log",
    run:
        shell("plink --genome {params.inPrefix} 2> {log}")
        

rule IBSRegion:
    input:
        genome = IN_PATH + "/population/subgroup/plink.genome",
    output:
        ibs = IN_PATH + "/population/subgroup/Samples_region_ibs.txt",
    params:
        kinshipSummary = SRC_DIR + "/kinshipSummary.py",
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/IBSRegion.log",
    run:
        shell("python {params.kinshipSummary} --input {input.genome} --meta {params.metafile} --out {output.ibs} > {log} 2>&1")



rule IBSPlot:
    input:
        ibs = IN_PATH + "/population/subgroup/Samples_region_ibs.txt",
    output:
        ibs = IN_PATH + "/population/subgroup/Samples_region_ibs.pdf",
    params:
        IBSDensity = SCRIPT_DIR + "/IBSDensity.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/IBSPlot.log",
    run:
        shell("Rscript {params.IBSDensity} --input {input.ibs} --pdf {output.ibs} --width {params.width} --height {params.height} > {log} 2>&1")
#########################################################



############################## Merge #############################
rule MergeSV:
    input:
        vcf = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=African),
    output:
        vcf = IN_PATH + "/population/Merge/African_SV_common.vcf",
    params:
        clique_maxflow_SV = SRC_DIR + "/clique_maxflow_SV.py",
        allele_freq = 0.2,
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/MergeSV.log",
    run:
        Files = ",".join(sorted(input.vcf))
        shell("python {params.clique_maxflow_SV} --workers {threads} -v {Files} -o {output.vcf} --allele_freq {params.allele_freq} > {log} 2>&1")

rule MergeSV2:
    input:
        vcf = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=American),
    output:
        vcf = IN_PATH + "/population/Merge/American_SV_common.vcf",
    params:
        clique_maxflow_SV = SRC_DIR + "/clique_maxflow_SV.py",
        allele_freq = 0.2,
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/MergeSV2.log",
    run:
        Files = ",".join(sorted(input.vcf))
        shell("python {params.clique_maxflow_SV} --workers {threads} -v {Files} -o {output.vcf} --allele_freq {params.allele_freq} > {log} 2>&1")
    

rule MergeSV3:
    input:
        vcf = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=Chinese),
    output:
        vcf = IN_PATH + "/population/Merge/Chinese_SV_common.vcf",
    params:
        clique_maxflow_SV = SRC_DIR + "/clique_maxflow_SV.py",
        allele_freq = 0.2,
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/MergeSV3.log",
    run:
        Files = ",".join(sorted(input.vcf))
        shell("python {params.clique_maxflow_SV} --workers {threads} -v {Files} -o {output.vcf} --allele_freq {params.allele_freq} > {log} 2>&1")

rule MergeSV4:
    input:
        vcf = expand(IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf", sample=Tibetan),
    output:
        vcf = IN_PATH + "/population/Merge/Tibetan_SV_common.vcf",
    params:
        clique_maxflow_SV = SRC_DIR + "/clique_maxflow_SV.py",
        allele_freq = 0.2,
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/MergeSV4.log",
    run:
        Files = ",".join(sorted(input.vcf))
        shell("python {params.clique_maxflow_SV} --workers {threads} -v {Files} -o {output.vcf} --allele_freq {params.allele_freq} > {log} 2>&1")

###########################################################################################



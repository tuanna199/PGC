
########################################## get represented seq #####################################


rule SVSeq:
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf",
    output:
        seq = IN_PATH + "/SVCall/FinalFilt/sequence/{sample}_filt_centromere.fasta",
    threads:
        THREADS
    params:
        VCFSVSeq = SRC_DIR + "/VCFSVSeq.py",
    log:
        IN_PATH + "/log/SVSeq_{sample}.log"
    run:
        shell("python {params.VCFSVSeq} --vcf {input.vcf} --out {output.seq} > {log} 2>&1")



rule TagSequence:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
        seq = expand(IN_PATH + "/SVCall/FinalFilt/sequence/{sample}_filt_centromere.fasta", sample=SAMPLES),
    output:
        seq = IN_PATH + "/population/Merge/Sample_common_SV_tag_seq.fasta",
    params:
        CombineSVSeqID = SRC_DIR + "/CombineSVSeqID.py",
    log:
        IN_PATH + "/log/TagSequence.log"
    run:
        SEQ = ",".join(input.seq)
        shell("python {params.CombineSVSeqID} --vcf {input.vcf} --sequence {SEQ} --out {output.seq} > {log} 2>&1")


rule SplitFasta:
    input:
        seq = IN_PATH + "/population/Merge/Sample_common_SV_tag_seq.fasta",
    output:
        DEL = IN_PATH + "/population/Repeat/Sample_common_SV_tag_DEL.fasta",
        INS = IN_PATH + "/population/Repeat/Sample_common_SV_tag_INS.fasta",
    params:
        SplitTypeFasta = SRC_DIR + "/SplitTypeFasta.py",
        outPrefix = IN_PATH + "/population/Repeat/Sample_common_SV_tag",
    log:
        IN_PATH + "/log/SplitFasta.log"
    run:
        shell("python {params.SplitTypeFasta} --fasta {input.seq} --outPrefix {params.outPrefix} > {log} 2>&1")


rule DUPINVSeq:
    input:
        DUP = IN_PATH + "/population/bed/Sample_common_SV_DUP.bed",
        INV = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
    output:
        DUP = IN_PATH + "/population/Repeat/Sample_common_SV_tag_DUP.fasta",
        INV = IN_PATH + "/population/Repeat/Sample_common_SV_tag_INV.fasta",
    params:
        RefGenome = config["RefGenome"],
    log:
        IN_PATH + "/log/DUPINVSeq.log"
    run:
        shell("bedtools getfasta -fi {params.RefGenome} -bed {input.DUP} -fo {output.DUP} > {log} 2>&1")
        shell("bedtools getfasta -fi {params.RefGenome} -bed {input.INV} -fo {output.INV} >> {log} 2>&1")




rule repeatMask:
    input:
        fa = IN_PATH + "/population/Repeat/Sample_common_SV_tag_{SVtype}.fasta",
    output:
        out = IN_PATH + "/population/Repeat/{SVtype}/Sample_common_SV_tag_{SVtype}.fasta.out",
    params:
        outdir = IN_PATH + "/population/Repeat/{SVtype}",
    threads:
        THREADS * ThreadFold
    run:
        cmd = "source activate WGS && RepeatMasker -parallel %s -species human -html -gff -dir repeat %s" % (threads, input.fa)
        print(cmd)
        os.system(cmd)
        


##########################################################################################




################################### Dfam ########################################
rule SeqDfamHMM:
    input:
        seq = IN_PATH + "/population/Merge/Sample_common_SV_tag_seq.fasta",
    output:
        tblout = IN_PATH + "/population/Dfam/SV_seq_tblout.txt",
        domtblout = IN_PATH + "/population/Dfam/SV_seq_domtblout.txt",
        pfamtblout = IN_PATH + "/population/Dfam/SV_seq_pfamtblout.txt",
        hmm = IN_PATH + "/population/Dfam/SV_seq_hmm.txt",
    threads:
        THREADS * ThreadFold
    params:
        DfamHmm = config["DfamHmm"],
    log:
        IN_PATH + "/log/SeqDfamHMM.log"
    run:
        shell('echo "hmmsearch --cpu 60 -E 10 --tblout {output.tblout} --domtblout {output.domtblout} --pfamtblout {output.pfamtblout}  -o {output.hmm} {params.DfamHmm}  {input.seq} > {log}" ')


rule HmmDfamType:
    input:
        tblout = IN_PATH + "/population/Dfam/SV_seq_tblout.txt",
    output:
        dfam = IN_PATH + "/population/Dfam/SV_seq_dfam_type.xls",
    params:
        hmmTagClass = SRC_DIR + "/hmmTagClass.py",
        DfamType = config["DfamType"],
    log:
        IN_PATH + "/log/HmmDfamType.log"
    run:
        shell("python {params.hmmTagClass} --dfam {params.DfamType} --input {input.tblout} --out {output.dfam} > {log} 2>&1")


rule HmmDfamTypePlot:
    input:
        dfam = IN_PATH + "/population/Dfam/SV_seq_dfam_type.xls",
    output:
        INS = IN_PATH + "/population/Dfam/SV_seq_dfam_INS_pie.pdf",
        DEL = IN_PATH + "/population/Dfam/SV_seq_dfam_DEL_pie.pdf",
        length1 = IN_PATH + "/population/Dfam/SV_seq_dfam_length_10000.pdf",
        length2 = IN_PATH + "/population/Dfam/SV_seq_dfam_length_1000.pdf",
    params:
        SVLengthContentPie = SCRIPT_DIR + "/SVLengthContentPie.R",
        outPrefix = IN_PATH + "/population/Dfam/SV_seq_dfam",
    log:
        IN_PATH + "/log/HmmDfamTypePlot.log"
    run:
        shell("Rscript {params.SVLengthContentPie} --input {input.dfam} --outPrefix {params.outPrefix} > {log} 2>&1")

# #####################################################################################################
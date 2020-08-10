
###########################################################################
rule SVFiltStats:
    ### just output file types: DEL, INS, DUP, INV and TRA
    input:
        vcf = IN_PATH + "/SVCall/FinalFilt/{sample}_filt_centromere.vcf",
    output:
        summary = IN_PATH + "/SVCall/FiltStats/{sample}_summary.xls",
        record = IN_PATH + "/SVCall/FiltStats/{sample}_record.xls",
        position = IN_PATH + "/SVCall/FiltStats/{sample}_position.xls",
        tra = IN_PATH + "/SVCall/FiltStats/{sample}_tra.xls",
    threads:
        THREADS
    params:
        SnifflesSVStats = SRC_DIR + "/SnifflesSVStats.py",
    log:
        IN_PATH + "/log/SVFiltStats_{sample}.log", 
    run:
        shell("python {params.SnifflesSVStats} --vcf {input.vcf} --summary {output.summary} --record {output.record} --position {output.position} --sample {wildcards.sample} --tra {output.tra} >{log} 2>&1")



rule ReadStatsSummary:
    input:
        summary = expand(IN_PATH + "/SVCall/FiltStats/{sample}_summary.xls", sample=SAMPLES),
    output:
        summary = IN_PATH + "/SVCall/FiltStats/Samples_SV_type_number.xls",
    threads:
        THREADS
    run:
        merge_multiple_sample_stats(input.summary, output.summary)


################################################################################


##################################### SV Length  ############################################

def SV_tag_length(tag_file, outPrefix):
    """
    Tag     CN001   CN002   CN003   CN004   CN005   CN006   CN007   CN008   CN009   CN010   CN011   CN012   CN013   CN014   CN015   CN016   
    1_66288-1_66527-239.0-DEL       -       -       -       -       -       -       -       -       -       -       -       -       -       
    1_67910-1_68341-431.0-DEL       -       -       -       -       -       -       -       -       -       -       -       -       -  
    """
    outdir = "/".join(outPrefix.split("/")[:-1])
    outdir = outdir + "/"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    ins_h = open(outPrefix + "_INS_DEL.txt", "w")
    inv_h = open(outPrefix + "_INV_DUP.txt", "w")
    ins_h.write("Tag\tSVType\tSVLength\n")
    inv_h.write("Tag\tSVType\tSVLength\n")
    
    tag_h = open(tag_file, "r")
    header = tag_h.readline().strip()
    for line in tag_h:
        lines = line.strip().split("\t")
        tag = lines[0]
        tags = tag.split("-")
        length = tags[2]
        SVType = tags[3]
        if SVType == "INS" or SVType == "DEL":
            ins_h.write("%s\t%s\t%s\n" % (tag, SVType, length))
        elif SVType == "INV" or SVType == "DUP":
            inv_h.write("%s\t%s\t%s\n" % (tag, SVType, length))
        else:
            print("Please ckeck whether INS, DEL, INV or DUP is in description %s." % tag)
    tag_h.close()
    inv_h.close()
    ins_h.close()



rule SVTypeLength:
    input:
        tag = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    output:
        INS = IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL.txt",
        INV = IN_PATH + "/population/Merge/Length/Sample_common_SV_INV_DUP.txt",
    params:
        outPrefix = IN_PATH + "/population/Merge/Length/Sample_common_SV",
    run:
        SV_tag_length(input.tag, params.outPrefix)



rule SVTypeLengthPlot:
    input:
        INS = IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL.txt",
        INV = IN_PATH + "/population/Merge/Length/Sample_common_SV_INV_DUP.txt",
    output:
        pdf1 = IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL_length_1000.pdf",
        pdf2 = IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL_length_10000.pdf",
        pdf3 = IN_PATH + "/population/Merge/Length/Sample_common_SV_INV_DUP_length_1000.pdf",
        pdf4 = IN_PATH + "/population/Merge/Length/Sample_common_SV_INV_DUP_length_10000.pdf",
    params:
        LengthHistDist = SCRIPT_DIR + "/LengthHistDist.R",
        outPrefix = IN_PATH + "/population/Merge/Length/Sample_common_SV_INS_DEL",
        outPrefix2 = IN_PATH + "/population/Merge/Length/Sample_common_SV_INV_DUP",
    log:
        IN_PATH + "/log/SVTypeLengthPlot.log"
    run:
        shell("Rscript {params.LengthHistDist} --input {input.INS} --outPrefix {params.outPrefix} > {log} 2>&1")
        shell("Rscript {params.LengthHistDist} --input {input.INV} --outPrefix {params.outPrefix2} > {log} 2>&1")
#################################################################################################



#################################### individuals per SV ######################################
rule SVIndividuals:
    input:
        tag = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    output:
        stat = IN_PATH + "/population/Stats/Sample_SV_type_tag_individulas.xls",
        each = IN_PATH + "/population/Stats/Sample_SV_type_tag_indiv_each.xls",
    params:
        NumberPerSV = SRC_DIR + "/NumberPerSV.py",
    log:
        IN_PATH + "/log/SVIndividuals.log"
    run:
        shell("python {params.NumberPerSV} --tag {input.tag} --out {output.stat} --list {output.each} > {log} 2>&1")
    

rule SVIndividualsPlot:
    input:
        stat = IN_PATH + "/population/Stats/Sample_SV_type_tag_individulas.xls",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_SV_type_tag_individulas.pdf",
    params:
        IndividualsPerSV = SCRIPT_DIR + "/IndividualsPerSV.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/SVIndividualsPlot.log"
    run:
        shell("Rscript {params.IndividualsPerSV} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
##############################################################################################



######################################################################################
rule Pop2bed0:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        Del = IN_PATH + "/population/bed/Sample_common_SV_DEL.bed",
        Ins = IN_PATH + "/population/bed/Sample_common_SV_INS.bed",
        Dup = IN_PATH + "/population/bed/Sample_common_SV_DUP.bed",
        Inv = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
        Other = IN_PATH + "/population/bed/Sample_common_SV_OTHER.bed",
        Tra = IN_PATH + "/population/bed/Sample_common_SV_TRA.bed",
        bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
        length = IN_PATH + "/population/bed/Sample_common_SV_length.xls",
    params:
        SVvcf2bedType = SRC_DIR + "/SVvcf2bedType.py",
        distance = 0, #config["IGV_distance"],
        insDistance = 0,
        outPrefix = IN_PATH + "/population/bed/Sample_common_SV",
        outdir = IN_PATH + "/population/bed/",
    log:
        IN_PATH + "/log/Pop2bed0.log",
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        shell("python {params.SVvcf2bedType} --vcf {input.vcf} --out {params.outPrefix} --insDistance {params.insDistance}  --distance {params.distance} --lengthOut {output.length} >{log} 2>&1")
        shell("cat {output.Del} {output.Ins} {output.Dup} {output.Inv} {output.Other} | sort -k 1,1  -k 2,2n > {output.bed}")



rule INSPoint:
    input:
        Ins = IN_PATH + "/population/bed/Sample_common_SV_INS.bed",
        Del = IN_PATH + "/population/bed/Sample_common_SV_DEL.bed",
        Dup = IN_PATH + "/population/bed/Sample_common_SV_DUP.bed",
        Inv = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
        Other = IN_PATH + "/population/bed/Sample_common_SV_OTHER.bed",
    output:
        Ins = IN_PATH + "/population/bed/Sample_common_SV_INS_point.bed",
        bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP_point.bed",
    run:
        ### add 5 bp for start
        cmd = """awk '{print $1"\t"$2"\t"$2+3"\t"$4}' %s > %s""" % (input.Ins, output.Ins)
        os.system(cmd)
        shell("cat {input.Del} {output.Ins} {input.Dup} {input.Inv} {input.Other} | sort -k 1,1  -k 2,2n > {output.bed}")



############################################################################################



############################################  SV general Stats ###############################################
rule SVFiltStatsAll:
    ### just output file types: DEL, INS, DUP, INV and TRA
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        summary = IN_PATH + "/population/Stats/Sample_merge_summary.xls",
        record = IN_PATH + "/population/Stats/Sample_merge_record.xls",
        position = IN_PATH + "/population/Stats/Sample_merge_position.xls",
        tra = IN_PATH + "/population/Stats/Sample_merge_tra.xls",
    threads:
        THREADS
    params:
        SnifflesSVStats = SRC_DIR + "/SnifflesSVStats.py",
    log:
        IN_PATH + "/log/SVStatsAll.log", 
    run:
        shell("python {params.SnifflesSVStats} --vcf {input.vcf} --summary {output.summary} --record {output.record} --position {output.position} --sample Samples --tra {output.tra} >{log} 2>&1")


rule TypeDensityAllPlot:
    input:
        record = IN_PATH + "/population/Stats/Sample_merge_record.xls",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_merge_record_length_density.pdf",
    threads:
        THREADS
    params:
        SVDensity = SCRIPT_DIR + "/SVDensity.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/TypeDensityAll.log", 
    run:
        # shell("source activate Rmeta && Rscript {params.SVDensity} --input {input.record} --pdf {output.pdf} --width {params.width} --height {params.height} >{log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.SVDensity, input.record, output.pdf, params.width, params.height, log)
        print(cmd)
        os.system(cmd)



rule SVDisSummaryAll:
    input:
        position = IN_PATH + "/population/Stats/Sample_merge_position.xls",
    output:
        summary = IN_PATH + "/population/Stats/Sample_merge_dist_summary.xls",
    threads:
        THREADS
    params:
        SVDistribution = SRC_DIR + "/SVDistribution.py",
        chromosome = config["chromosome"],
        SVType = config["SVType"],
    log:
        IN_PATH + "/log/SVDisSummaryAll.log", 
    run:
        shell("python {params.SVDistribution} --input {input.position} --summary {output.summary} --chromosome {params.chromosome} --type {params.SVType} >{log} 2>&1")


rule SVSummaryPlotAll:
    input:
        summary = IN_PATH + "/population/Stats/Sample_merge_dist_summary.xls",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_merge_chroms_distribution.pdf",
    threads:
        THREADS
    params:
        SVTypeBarStack = SCRIPT_DIR + "/SVTypeBarStack.R",
        chromosome = config["chromosome"],
        width = 6,
        height = 4.5,
    log:
        IN_PATH + "/log/SVSummaryPlotAll.log", 
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s --chromosome %s > %s 2>&1" % (params.SVTypeBarStack, input.summary, output.pdf, params.width, params.height, params.chromosome, log)
        os.system(cmd)



rule SVChrLengthCor:
    input:
        summary = IN_PATH + "/population/Stats/Sample_merge_dist_summary.xls",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_merge_SV_number_ChrLength_cor.pdf",
    params:
        ChrLenSVCor = SCRIPT_DIR + "/ChrLenSVCor.R",
        GenomeLenChr = config["GenomeLenChr"],
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/SVChrLengthCor.log", 
    run:
        # shell("Rscript {params.ChrLenSVCor} --SV {input.summary} --length {params.GenomeLenChr} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
        ### something wrong with ggblur in NanoSV
        cmd = "source activate Rmeta && Rscript %s --SV %s --length %s --pdf %s --width %s --height %s > %s 2>&1" % (params.ChrLenSVCor, input.summary, params.GenomeLenChr, output.pdf, params.width, params.height, log)
        os.system(cmd)


rule SVDisPlotAll:
    input:
        position = IN_PATH + "/population/Stats/Sample_merge_position.xls",
    output:
        position = temp(IN_PATH + "/population/Stats/Sample_merge_position_temp.xls",),
        pdf = IN_PATH + "/population/Stats/Sample_merge_distribution.pdf",
    threads:
        THREADS
    params:
        SVDistribution = SCRIPT_DIR + "/SVDistribution.R",
        maxNum = 300,
        GenomeLen = config["GenomeLen"],
    log:
        IN_PATH + "/log/SVDisPlotAll.log", 
    run:
        # ### set the max length of the chromosome
        shell("cat {input.position} {params.GenomeLen} > {output.position}")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --maxNum %s > %s 2>&1" % (params.SVDistribution, output.position, output.pdf, params.maxNum, log)
        os.system(cmd)



rule TypeStackAll:
    input:
        # summary = IN_PATH + "/population/Stats/Sample_merge_dist_summary.xls",
        summary = IN_PATH + "/population/Stats/Sample_merge_summary.xls",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_merge_dist_summary.pdf",
    threads:
        THREADS
    params:
        SVTypeBar = SCRIPT_DIR + "/SVTypeBar.R",
        width = 3,
        height = 4,
    log:
        IN_PATH + "/log/TypeStackAll.log", 
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.SVTypeBar, input.summary, output.pdf, params.width, params.height, log)
        print(cmd)
        os.system(cmd)



rule SVLength:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        length = IN_PATH + "/population/Stats/Sample_common_SV_chr_length.xls",
    params:
        ChromSVLength = SRC_DIR + "/ChromSVLength.py",
        cutLength = 50000000,
    log:
        IN_PATH + "/log/SVLength.log",
    run:
         shell("python {params.ChromSVLength} --vcf {input.vcf} --out {output.length} --cutLength {params.cutLength} --LenTag 'AVGLEN' > {log} 2>&1")


rule SVLengthPlot:
    input:
        length = IN_PATH + "/population/Stats/Sample_common_SV_chr_length.xls",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_common_SV_chr_length.pdf",
    params:
        SVLengthBarStack = SCRIPT_DIR + "/SVLengthBarStack.R",
        width = 6,
        height = 4.5,
    log:
        IN_PATH + "/log/SVLengthPlot.log",
    run:
        shell("Rscript {params.SVLengthBarStack} --input {input.length} --pdf {output.pdf} --height {params.height} --width {params.width} > {log} 2>&1")


######################################################################################################


########################################## position circos #####################
rule positionWindow:
    input:
        position = IN_PATH + "/population/Stats/Sample_merge_position.xls",
    output:
        DEL = IN_PATH + "/population/Stats/Sample_DEL.xls",
        INS = IN_PATH + "/population/Stats/Sample_INS.xls",
        DUP = IN_PATH + "/population/Stats/Sample_DUP.xls",
        INV = IN_PATH + "/population/Stats/Sample_INV.xls",
    params:
        WindowSVNumber = SRC_DIR + "/WindowSVNumber.py",
        GenomeLen = config["GenomeLen"],
        window = 500000,
        sliding = 0,
        outprefix = IN_PATH + "/population/Stats/Sample",
    log:
        IN_PATH + "/log/positionWindow.log", 
    run:
        shell("python {params.WindowSVNumber} --genome {params.GenomeLen} --input {input.position} --window {params.window} --sliding {params.sliding} --out {params.outprefix} > {log} 2>&1")


rule Positions:
    input:
        DEL = IN_PATH + "/population/Stats/Sample_DEL.xls",
        INS = IN_PATH + "/population/Stats/Sample_INS.xls",
        DUP = IN_PATH + "/population/Stats/Sample_DUP.xls",
        INV = IN_PATH + "/population/Stats/Sample_INV.xls",
    output:
        DEL = IN_PATH + "/population/Stats/Sample_DEL_circos.xls",
        INS = IN_PATH + "/population/Stats/Sample_INS_circos.xls",
        DUP = IN_PATH + "/population/Stats/Sample_DUP_circos.xls",
        INV = IN_PATH + "/population/Stats/Sample_INV_circos.xls",
    threads:
        THREADS
    params:
        TRA2circos = SRC_DIR + "/TRA2circos.py",
        chromosome = config["chromosome"],
    log:
        IN_PATH + "/log/Positions.log", 
    run:
        shell("python {params.TRA2circos} --input  {input.DEL} --out {output.DEL} --chromosome {params.chromosome} --method one >{log} 2>&1")
        shell("python {params.TRA2circos} --input  {input.INS} --out {output.INS} --chromosome {params.chromosome} --method one 2>>{log}")
        shell("python {params.TRA2circos} --input  {input.DUP} --out {output.DUP} --chromosome {params.chromosome} --method one 2>>{log}")
        shell("python {params.TRA2circos} --input  {input.INV} --out {output.INV} --chromosome {params.chromosome} --method one 2>>{log}")




rule TRA2circos:
    input:
        tra = IN_PATH + "/population/Stats/Sample_merge_tra.xls",
    output:
        circos = IN_PATH + "/population/Stats/Sample_merge_tra_circos.xls",
    threads:
        THREADS
    params:
        TRA2circos = SRC_DIR + "/TRA2circos.py",
        chromosome = config["chromosome"],
    log:
        IN_PATH + "/log/TRA2circos.log", 
    run:
        shell("python {params.TRA2circos} --input  {input.tra} --out {output.circos} --chromosome {params.chromosome} --method two >{log} 2>&1")
################################################################################






######################################  SV density for different types ################################
rule TerminalDistanceAll:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        dist = IN_PATH + "/population/SVDensity/Sample_common_to_terminal_distance.xls",
    params:
        SVTerminalDistance = SRC_DIR + "/SVTerminalDistance.py",
        cytoband = config["cytoband"],
        window = 100000, ## 100000,
        sliding = 0,
    log:
        IN_PATH + "/log/TerminalDistanceAll.log", 
    run:
        shell("python {params.SVTerminalDistance} --vcf {input.vcf} --cytoband {params.cytoband} --out {output.dist} --window  {params.window} --sliding {params.sliding} >{log} 2>&1")


rule TerminalDistancePlotAll:
    input:
        dist = IN_PATH + "/population/SVDensity/Sample_common_to_terminal_distance.xls",
    output:
        pdf = IN_PATH + "/population/SVDensity/Sample_common_to_terminal_distance.pdf",
    params:
        SVTelomereDistance = SCRIPT_DIR + "/SVTelomereDistance.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/TerminalDistancePlotAll.log", 
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.SVTelomereDistance, input.dist, output.pdf, params.width, params.height, log)
        os.system(cmd)


rule SVTypeDensity:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_filt.vcf",
    output:
        DEL = IN_PATH + "/population/SVDensity/Sample_common_SV_{SVtype}.vcf",
    run:
        cmd = "grep '<%s>' %s > %s" % (wildcards.SVtype, input.vcf, output.DEL)
        os.system(cmd)




rule TerminalDistanceType:
    input:
        vcf = IN_PATH + "/population/SVDensity/Sample_common_SV_{SVtype}.vcf",
    output:
        dist = IN_PATH + "/population/SVDensity/Sample_{SVtype}_to_terminal_distance.xls",
    params:
        SVTerminalDistance = SRC_DIR + "/SVTerminalDistance.py",
        cytoband = config["cytoband"],
        window = 100000, ## 100000,
        sliding = 0,
    log:
        IN_PATH + "/log/CentroDistance_{SVtype}.log", 
    run:
        shell("python {params.SVTerminalDistance} --vcf {input.vcf} --cytoband {params.cytoband} --out {output.dist} --window  {params.window} --sliding {params.sliding} >{log} 2>&1")



rule TerminalDistancePlotType:
    input:
        dist = IN_PATH + "/population/SVDensity/Sample_{SVtype}_to_terminal_distance.xls",
    output:
        pdf = IN_PATH + "/population/SVDensity/Sample_{SVtype}_to_terminal_distance.pdf",
    params:
        SVTelomereDistance = SCRIPT_DIR + "/SVTelomereDistance.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/TerminalDistancePlot_{SVtype}.log", 
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.SVTelomereDistance, input.dist, output.pdf, params.width, params.height, log)
        os.system(cmd)
#####################################################################################################################




##############################  SV Distance ##################################


rule lofSVDistAll:
    input:
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
    output:
        temp1 = IN_PATH + "/population/Stats/Sample_SV_distance_temp1.txt",
        temp2 = IN_PATH + "/population/Stats/Sample_SV_distance_temp2.txt",
        dist = IN_PATH + "/population/Stats/Sample_SV_distance.txt",
    run:
        shell("sed '1d' {input.freq} | cut -f 1 > {output.temp1}")
        shell("sed '1d' {input.freq} | cut -f 1 | cut -f 3 -d '-' > {output.temp2}")
        shell("paste {output.temp1} {output.temp2} > {output.dist}")



rule CategoryLengthAll:
    input:
        dist = IN_PATH + "/population/Stats/Sample_SV_distance.txt",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        dist = IN_PATH + "/population/Stats/Sample_SV_distance_category.txt",
    run:
        catagory_SV_distance(input.category, input.dist, output.dist)



rule CategoryLengthAllPlot:
    input:
        dist = IN_PATH + "/population/Stats/Sample_SV_distance_category.txt",
    output:
        pdf = IN_PATH + "/population/Stats/Sample_SV_distance_category.pdf",
    params:
        lengthBox = SCRIPT_DIR + "/lengthBox.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/CategoryLengthAllPlot.log",
    run:
        shell("Rscript {params.lengthBox} --input {input.dist} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")




###############################################################################

################################ novel SV ##############################
rule NovelTag:
    input:
        tag1 = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap-1.xls",
        allTag = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    output:
        temp = temp(IN_PATH + "/population/bed/Sample_common_SV_Tag_temp.xls"),
        novel = IN_PATH + "/population/bed/SV_novel_tag.xls",
    run:
        shell("sed '1d' {input.allTag}  | cut -f 1 > {output.temp}")
        shell("cat {input.tag1} {output.temp} | sort | uniq -c |sed 's/^      //g' | grep '^1' | cut -f 2 -d ' ' > {output.novel}")




rule LongOverlap:
    input:
        Cell_Del = IN_PATH + "/preStudies/Cell2019_SV_dechr_DEL.bed",
        Cell_Ins = IN_PATH + "/preStudies/Cell2019_SV_dechr_INS_newEnd.bed",
        Cell_Inv = IN_PATH + "/preStudies/Cell2019_SV_dechr_INV.bed",
        HGSVC_Del =  IN_PATH + "/preStudies/HGSVC/HGSVC_combine_DEL.bed",
        HGSVC_Ins =  IN_PATH + "/preStudies/HGSVC/HGSVC_combine_INS.bed",
        HGSVC_Inv =  IN_PATH + "/preStudies/HGSVC/HGSVC_combine_INV.bed",
    output:
        Del = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_DEL.bed",
        Ins = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_INS.bed",
        Inv = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_INV.bed",
    run:
        shell("bedtools intersect -a {input.HGSVC_Del} -b {input.Cell_Del}  -wa -wb > {output.Del}")
        shell("bedtools intersect -a {input.HGSVC_Ins} -b {input.Cell_Ins}  -wa -wb > {output.Ins}")
        shell("bedtools intersect -a {input.HGSVC_Inv} -b {input.Cell_Inv}  -wa -wb > {output.Inv}")



rule LongOverlapFilt:
    input:
        Del = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_DEL.bed",
        Ins = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_INS.bed",
        Inv = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_INV.bed",
    output:
        Del = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_DEL_filt.bed",
        Ins = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_INS_filt.bed",
        Inv = IN_PATH + "/population/bed/DBOverlap/HGSVC_Cell_overlap_INV_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/LongOverlapFilt.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")
        shell("python {params.bedOverlapRecord} --bed {input.Ins} --ratioThreshold {params.ratioThreshold} --out {output.Ins} --method both 2>>{log}")
        shell("python {params.bedOverlapRecord} --bed {input.Inv} --ratioThreshold {params.ratioThreshold} --out {output.Inv} --method both 2>>{log}")

#########################################################################




################ INV overlap ###################
rule INVOverlap:
    input:
        Inv = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
        nstd152 = "/home/wuzhikun/database/HGSVC2/Chaisson_nstd152/bed/nstd152.GRCh38.variant_call_INV.bed",
        nstd169 = "/home/wuzhikun/database/HGSVC2/Giner-Delgado_nstd169/bed/nstd169.GRCh38_INV.bed",
        InvFEST = "/home/wuzhikun/database/HGSVC2/InvFEST/InvFEST_hg38_INV_noFalse_dechr.txt",
    output:
        nstd152 = IN_PATH + "/population/bed/INV/nstd152_INV_intersect.bed",
        nstd169 = IN_PATH + "/population/bed/INV/nstd169_INV_intersect.bed",
        InvFEST = IN_PATH + "/population/bed/INV/InvFEST_INV_intersect.bed",
    run:
        shell("bedtools intersect -a {input.Inv} -b {input.nstd152}  -wa -wb > {output.nstd152}")
        shell("bedtools intersect -a {input.Inv} -b {input.nstd169}  -wa -wb > {output.nstd169}")
        shell("bedtools intersect -a {input.Inv} -b {input.InvFEST}  -wa -wb > {output.InvFEST}")





rule LongOverlapFilt2:
    input:
        nstd152 = IN_PATH + "/population/bed/INV/nstd152_INV_intersect.bed",
        nstd169 = IN_PATH + "/population/bed/INV/nstd169_INV_intersect.bed",
        InvFEST = IN_PATH + "/population/bed/INV/InvFEST_INV_intersect.bed",
    output:
        nstd152 = IN_PATH + "/population/bed/INV/nstd152_INV_overlap.bed",
        nstd169 = IN_PATH + "/population/bed/INV/nstd169_INV_overlap.bed",
        InvFEST = IN_PATH + "/population/bed/INV/InvFEST_INV_overlap.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/LongOverlapFilt2.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.nstd152} --ratioThreshold {params.ratioThreshold} --out {output.nstd152} --method both > {log} 2>&1")
        shell("python {params.bedOverlapRecord} --bed {input.nstd169} --ratioThreshold {params.ratioThreshold} --out {output.nstd169} --method both 2>>{log}")
        shell("python {params.bedOverlapRecord} --bed {input.InvFEST} --ratioThreshold {params.ratioThreshold} --out {output.InvFEST} --method both 2>>{log}")


rule combineINV:
    input:
        # Inv = IN_PATH + "/population/bed/database/database_SV_INV.bed",
        LRS15 = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap_filt.bed",
        DGV = IN_PATH + "/population/bed/DGV/Sample_INV_overlap_filt.bed",
        gnomAD = IN_PATH + "/population/bed/gnomAD/Sample_INV_overlap_filt.bed",
        WGS911 = IN_PATH + "/population/bed/WGS911/Sample_INV_overlap_filt.bed",
        HGSVC = IN_PATH + "/population/bed/HGSVC/Sample_INV_overlap_filt.bed",
        nstd152 = IN_PATH + "/population/bed/INV/nstd152_INV_overlap.bed",
        nstd169 = IN_PATH + "/population/bed/INV/nstd169_INV_overlap.bed",
        InvFEST = IN_PATH + "/population/bed/INV/InvFEST_INV_overlap.bed",
    output:
        INV = IN_PATH + "/population/bed/INV/All_INV_overlap_unique.bed",
        length = IN_PATH + "/population/bed/INV/All_INV_overlap_unique_length.txt",
    run:
        shell("cat {input.LRS15} {input.DGV} {input.gnomAD} {input.WGS911} {input.HGSVC} {input.nstd152} {input.nstd169} {input.InvFEST} | cut -f 1-4 | sort | uniq | sort -k 1,1n -k 2,2n > {output.INV}")
        shell("cut -f 4  {output.INV} | cut -f 3 -d '-' > {output.length}")


################################################




################### Overlap Population #########
rule PopSV:
    input:
        pop = IN_PATH +  "/preStudies/population/combine/FivePop_{region}_SV.txt",
    output:
        bed = IN_PATH +  "/preStudies/population/bed/{region}_combine_SV_{svtype}.txt",
    run:
        shell(" grep {wildcards.svtype} {input.pop} | cut -f 1-4 | sort -k 1,1n -k 2,2n > {output.bed}")



rule PopSVOverlap:
    input:
        sv = IN_PATH + "/population/bed/Sample_common_SV_{svtype}.bed",
        bed = IN_PATH +  "/preStudies/population/bed/{region}_combine_SV_{svtype}.txt",
    output:
        bed = IN_PATH +  "/preStudies/population/bed/Overlap/{region}_combine_overlap_SV_{svtype}.txt",
        filt = IN_PATH +  "/preStudies/population/bed/Overlap/{region}_combine_overlap_SV_{svtype}_filt.txt",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/PopSVOverlap_{region}_{svtype}.log",
    run:
        shell("bedtools intersect -a {input.sv} -b {input.bed} -wa -wb > {output.bed}")
        shell("python {params.bedOverlapRecord} --bed {output.bed} --ratioThreshold {params.ratioThreshold} --out {output.filt} --method both > {log} 2>&1")




rule OverlapSum:
    input:
        DEL = IN_PATH +  "/preStudies/population/bed/Overlap/{region}_combine_overlap_SV_DEL_filt.txt",
        INS = IN_PATH +  "/preStudies/population/bed/Overlap/{region}_combine_overlap_SV_INS_filt.txt",
        DUP = IN_PATH +  "/preStudies/population/bed/Overlap/{region}_combine_overlap_SV_DUP_filt.txt",
        INV = IN_PATH +  "/preStudies/population/bed/Overlap/{region}_combine_overlap_SV_INV_filt.txt",
    output:
        tag = IN_PATH + "/preStudies/population/bed/Overlap/{region}_SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/preStudies/population/bed/Overlap/{region}_overlap_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/OverlapSum_{region}.log",
    run:
        BED = ",".join([input.DEL, input.INS, input.DUP, input.INV])
        shell("python {params.OverlapSVTypeLength} --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        BED2 = " ".join([input.DEL, input.INS, input.DUP, input.INV])
        shell("cat {BED2} | cut -f 4 | sort | uniq > {output.tag}")



rule OverlapTags:
    input:
        overlap = IN_PATH + "/preStudies/population/bed/Overlap/{region}_SV_overlap_filt_tags.xls",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        stats = IN_PATH + "/preStudies/population/bed/Overlap/{region}_kown_and_novel_tags_stats.txt",
        novel = IN_PATH + "/preStudies/population/bed/Overlap/{region}_novel_SV_category_tag.xls",
    params:
        CatNovelStats = SRC_DIR + "/CatNovelStats.py",
    log:
        IN_PATH + "/log/OverlapTags_{region}.log"
    run:
        shell("python {params.CatNovelStats} --category {input.category} --overlap {input.overlap} --stats {output.stats} --novel {output.novel} > {log} 2>&1")
################################################



############## Common category ###########
# rule commonBed:
#     input:
#         novel = IN_PATH + "/population/Category/Novel_SV_category_tag.xls",
#         Del = IN_PATH + "/population/bed/Sample_common_SV_DEL.bed",
#         Ins = IN_PATH + "/population/bed/Sample_common_SV_INS.bed",
#         Dup = IN_PATH + "/population/bed/Sample_common_SV_DUP.bed",
#         Inv = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
#     output:
#         common = IN_PATH + "/population/Category/Novel_common_tag.xls",
#         Del = IN_PATH + "/population/Category/Novel_common_DEL.bed",
#         Ins = IN_PATH + "/population/Category/Novel_common_INS.bed",
#         Dup = IN_PATH + "/population/Category/Novel_common_DUP.bed",
#         Inv = IN_PATH + "/population/Category/Novel_common_INV.bed",
#     run:
#         shell("grep Common  {input.novel} |cut -f 2 | sed 's/,/\n/g' > {output.common}")
#         shell("for i in `cat {input.novel}`; do grep $i {input.Del} > {output.Del}; done")


rule Chinesebed:
    input:
        vcf = "/home/wuzhikun/Project/Revised/population/Merge/Chinese_SV_common.vcf",
    output:
        Del = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_DEL.bed",
        Ins = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_INS.bed",
        Dup = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_DUP.bed",
        Inv = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_INV.bed",
        Other = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_OTHER.bed",
        Tra = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_TRA.bed",
        bed = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_DEL_INS_INV_DUP.bed",
        length = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_length.xls",
    params:
        SVvcf2bedType = SRC_DIR + "/SVvcf2bedType.py",
        distance = 0, #config["IGV_distance"],
        insDistance = 0,
        outPrefix = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV",
        outdir = "/home/wuzhikun/Project/Revised/population/Chinese",
    log:
        IN_PATH + "/log/Chinesebed.log",
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        shell("python {params.SVvcf2bedType} --vcf {input.vcf} --out {params.outPrefix} --insDistance {params.insDistance}  --distance {params.distance} --lengthOut {output.length} >{log} 2>&1")
        shell("cat {output.Del} {output.Ins} {output.Dup} {output.Inv} {output.Other} | sort -k 1,1  -k 2,2n > {output.bed}")



rule Africanbed:
    input:
        vcf = "/home/wuzhikun/Project/Revised/population/Merge/African_SV_common.vcf",
    output:
        Del = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_DEL.bed",
        Ins = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_INS.bed",
        Dup = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_DUP.bed",
        Inv = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_INV.bed",
        Other = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_OTHER.bed",
        Tra = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_TRA.bed",
        bed = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_DEL_INS_INV_DUP.bed",
        length = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_length.xls",
    params:
        SVvcf2bedType = SRC_DIR + "/SVvcf2bedType.py",
        distance = 0, #config["IGV_distance"],
        insDistance = 0,
        outPrefix = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV",
        outdir = "/home/wuzhikun/Project/Revised/population/African",
    log:
        IN_PATH + "/log/Africanbed.log",
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        shell("python {params.SVvcf2bedType} --vcf {input.vcf} --out {params.outPrefix} --insDistance {params.insDistance}  --distance {params.distance} --lengthOut {output.length} >{log} 2>&1")
        shell("cat {output.Del} {output.Ins} {output.Dup} {output.Inv} {output.Other} | sort -k 1,1  -k 2,2n > {output.bed}")




rule Americanbed:
    input:
        vcf = "/home/wuzhikun/Project/Revised/population/Merge/American_SV_common.vcf",
    output:
        Del = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_DEL.bed",
        Ins = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_INS.bed",
        Dup = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_DUP.bed",
        Inv = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_INV.bed",
        Other = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_OTHER.bed",
        Tra = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_TRA.bed",
        bed = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_DEL_INS_INV_DUP.bed",
        length = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_length.xls",
    params:
        SVvcf2bedType = SRC_DIR + "/SVvcf2bedType.py",
        distance = 0, #config["IGV_distance"],
        insDistance = 0,
        outPrefix = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV",
        outdir = "/home/wuzhikun/Project/Revised/population/American",
    log:
        IN_PATH + "/log/Americanbed.log",
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        shell("python {params.SVvcf2bedType} --vcf {input.vcf} --out {params.outPrefix} --insDistance {params.insDistance}  --distance {params.distance} --lengthOut {output.length} >{log} 2>&1")
        shell("cat {output.Del} {output.Ins} {output.Dup} {output.Inv} {output.Other} | sort -k 1,1  -k 2,2n > {output.bed}")



rule ChineseGenotype:
    input:
        vcf = "/home/wuzhikun/Project/Revised/population/Merge/Chinese_SV_common.vcf",
    output:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/Chinese_SV_genotype.txt",
        freq = "/home/wuzhikun/Project/Revised/population/Merge/genotype/Chinese_SV_genotype_numeric.txt",
    params:
        SVGenotype = SRC_DIR + "/SVGenotype.py",
    log:
        IN_PATH + "/log/SVGenotype.log", 
    run:
        shell("python {params.SVGenotype} --vcf {input.vcf} --out {output.geno} --frequency {output.freq} >{log} 2>&1")


rule ChineseFrequencyFilt:
    input:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/Chinese_SV_genotype.txt",
    output:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/Chinese_SV_genotype_filt.txt",
        freq = "/home/wuzhikun/Project/Revised/population/Merge/genotype/Chinese_SV_genotype_frequency.txt",
    params:
        SVGenoFreqFilt = SRC_DIR + "/SVGenoFreqFilt.py",
        freqThreshold = 0.05,
    log:
        IN_PATH + "/log/SVFrequencyFilt.log", 
    run:
        shell("python {params.SVGenoFreqFilt} --genotype {input.geno} --freqThreshold {params.freqThreshold} --out {output.geno} --frequency {output.freq} > {log} 2>&1")





rule AfricanGenotype:
    input:
        vcf = "/home/wuzhikun/Project/Revised/population/Merge/African_SV_common.vcf",
    output:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype.txt",
        freq = "/home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype_numeric.txt",
    params:
        SVGenotype = SRC_DIR + "/SVGenotype.py",
    log:
        IN_PATH + "/log/SVGenotype.log", 
    run:
        shell("python {params.SVGenotype} --vcf {input.vcf} --out {output.geno} --frequency {output.freq} >{log} 2>&1")




rule AfricanFrequencyFilt:
    input:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype.txt",
    output:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype_filt.txt",
        freq = "/home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype_frequency.txt",
    params:
        SVGenoFreqFilt = SRC_DIR + "/SVGenoFreqFilt.py",
        freqThreshold = 0.05,
    log:
        IN_PATH + "/log/SVFrequencyFilt.log", 
    run:
        shell("python {params.SVGenoFreqFilt} --genotype {input.geno} --freqThreshold {params.freqThreshold} --out {output.geno} --frequency {output.freq} > {log} 2>&1")





rule AmericanGenotype:
    input:
        vcf = "/home/wuzhikun/Project/Revised/population/Merge/American_SV_common.vcf",
    output:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype.txt",
        freq = "/home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype_numeric.txt",
    params:
        SVGenotype = SRC_DIR + "/SVGenotype.py",
    log:
        IN_PATH + "/log/SVGenotype.log", 
    run:
        shell("python {params.SVGenotype} --vcf {input.vcf} --out {output.geno} --frequency {output.freq} >{log} 2>&1")




rule AmericanFrequencyFilt:
    input:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype.txt",
    output:
        geno = "/home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype_filt.txt",
        freq = "/home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype_frequency.txt",
    params:
        SVGenoFreqFilt = SRC_DIR + "/SVGenoFreqFilt.py",
        freqThreshold = 0.05,
    log:
        IN_PATH + "/log/SVFrequencyFilt.log", 
    run:
        shell("python {params.SVGenoFreqFilt} --genotype {input.geno} --freqThreshold {params.freqThreshold} --out {output.geno} --frequency {output.freq} > {log} 2>&1")



rule PopulationOverlap:
    input:
        bed = IN_PATH + "/population/Category/Novel_common_{svtype}.bed",
        American = "/home/wuzhikun/Project/Revised/population/American/Sample_common_SV_{svtype}.bed",
        African = "/home/wuzhikun/Project/Revised/population/African/Sample_common_SV_{svtype}.bed",
        Chinese = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_common_SV_{svtype}.bed",
    output:
        American = "/home/wuzhikun/Project/Revised/population/American/Sample_overlap_SV_{svtype}.bed",
        African = "/home/wuzhikun/Project/Revised/population/African/Sample_overlap_SV_{svtype}.bed",
        Chinese = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_overlap_SV_{svtype}.bed",
        AmericanFilt = "/home/wuzhikun/Project/Revised/population/American/Sample_overlap_filt_SV_{svtype}.bed",
        AfricanFilt = "/home/wuzhikun/Project/Revised/population/African/Sample_overlap_filt_SV_{svtype}.bed",
        ChineseFilt = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_overlap_filt_SV_{svtype}.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/overlapCellFilt_{svtype}.log",
    run:
        shell("bedtools intersect -a {input.bed} -b {input.American} -wa -wb > {output.American}")
        shell("bedtools intersect -a {input.bed} -b {input.African} -wa -wb > {output.African}")
        shell("bedtools intersect -a {input.bed} -b {input.Chinese} -wa -wb > {output.Chinese}")
        shell("python {params.bedOverlapRecord} --bed {output.American} --ratioThreshold {params.ratioThreshold} --out {output.AmericanFilt} --method both > {log} 2>&1")
        shell("python {params.bedOverlapRecord} --bed {output.African} --ratioThreshold {params.ratioThreshold} --out {output.AfricanFilt} --method both > {log} 2>&1")
        shell("python {params.bedOverlapRecord} --bed {output.Chinese} --ratioThreshold {params.ratioThreshold} --out {output.ChineseFilt} --method both > {log} 2>&1")



rule PopulationTag:
    input:
        AmericanFilt = expand("/home/wuzhikun/Project/Revised/population/American/Sample_overlap_filt_SV_{svtype}.bed", svtype=SVTYPES),
        AfricanFilt = expand("/home/wuzhikun/Project/Revised/population/African/Sample_overlap_filt_SV_{svtype}.bed", svtype=SVTYPES),
        ChineseFilt = expand("/home/wuzhikun/Project/Revised/population/Chinese/Sample_overlap_filt_SV_{svtype}.bed", svtype=SVTYPES),
    output:
        American = "/home/wuzhikun/Project/Revised/population/American/Sample_overlap_tag.txt",
        African = "/home/wuzhikun/Project/Revised/population/African/Sample_overlap_tag.txt",
        Chinese = "/home/wuzhikun/Project/Revised/population/Chinese/Sample_overlap_tag.txt",
    run:
        American = " ".join(input.AmericanFilt)
        African = " ".join(input.AfricanFilt)
        Chinese = " ".join(input.ChineseFilt)
        shell("cat {American} | cut -f 8 | sort | uniq > {output.American}")
        shell("cat {African} | cut -f 8 | sort | uniq > {output.African}")
        shell("cat {Chinese} | cut -f 8 | sort | uniq > {output.Chinese}")
        ### for i in `cat /home/wuzhikun/Project/Revised/population/American/Sample_overlap_tag.txt`; do  grep $i /home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype_frequency.txt >> /home/wuzhikun/Project/Revised/population/Merge/genotype/American_SV_genotype_frequency_overlap.txt; done
        ### for i in `cat /home/wuzhikun/Project/Revised/population/African/Sample_overlap_tag.txt`; do  grep $i /home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype_frequency.txt >> /home/wuzhikun/Project/Revised/population/Merge/genotype/African_SV_genotype_frequency_overlap.txt; done


###############################################






################################## Overlap Asian ##############
rule OverlapAsian:
    input:
        DEL = IN_PATH + "/population/bed/Sample_common_SV_DEL.bed",
        INS = IN_PATH + "/population/bed/Sample_common_SV_INS.bed",
        INV = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
        DUP = IN_PATH + "/population/bed/Sample_common_SV_DUP.bed",
        CellDEL = IN_PATH + "/preStudies/population/Asian/Cell2019_Table_Asian_DEL.txt",
        CellINS = IN_PATH + "/preStudies/population/Asian/Cell2019_Table_Asian_INS.txt",
        CellINV = IN_PATH + "/preStudies/population/Asian/Cell2019_Table_Asian_INV.txt",
        gnomadDEL = IN_PATH + "/preStudies/population/Asian/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DEL.txt",
        gnomadINS = IN_PATH + "/preStudies/population/Asian/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INS.txt",
        gnomadINV = IN_PATH + "/preStudies/population/Asian/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INV.txt",
        gnomadDUP = IN_PATH + "/preStudies/population/Asian/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DUP.txt",
        HGSVCDEL = IN_PATH + "/preStudies/population/Asian/HGSVC2.East_Asian.reformat_DEL.txt",
        HGSVCINS = IN_PATH + "/preStudies/population/Asian/HGSVC2.East_Asian.reformat_INS.txt",
        HGDPDEL = IN_PATH + "/preStudies/population/Asian/WGS911_Asia_DEL.txt",
        HGDPINS = IN_PATH + "/preStudies/population/Asian/WGS911_Asia_INS.txt",
        HGDPINV = IN_PATH + "/preStudies/population/Asian/WGS911_Asia_INV.txt",
        HGDPDUP = IN_PATH + "/preStudies/population/Asian/WGS911_Asia_DUP.txt",
    output:
        CellDEL = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_DEL.txt",
        CellINS = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_INS.txt",
        CellINV = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_INV.txt",
        gnomadDEL = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DEL.txt",
        gnomadINS = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INS.txt",
        gnomadINV = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INV.txt",
        gnomadDUP = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DUP.txt",
        HGSVCDEL = IN_PATH + "/preStudies/population/Asian/overlap/HGSVC2.East_Asian.reformat_DEL.txt",
        HGSVCINS = IN_PATH + "/preStudies/population/Asian/overlap/HGSVC2.East_Asian.reformat_INS.txt",
        HGDPDEL = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_DEL.txt",
        HGDPINS = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_INS.txt",
        HGDPINV = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_INV.txt",
        HGDPDUP = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_DUP.txt",
    run:
        shell("bedtools intersect -a {input.DEL} -b {input.CellDEL}  -wa -wb > {output.CellDEL}")
        shell("bedtools intersect -a {input.DEL} -b {input.gnomadDEL}  -wa -wb > {output.gnomadDEL}")
        shell("bedtools intersect -a {input.DEL} -b {input.HGSVCDEL}  -wa -wb > {output.HGSVCDEL}")
        shell("bedtools intersect -a {input.DEL} -b {input.HGDPDEL}  -wa -wb > {output.HGDPDEL}")
        shell("bedtools intersect -a {input.INS} -b {input.CellINS}  -wa -wb > {output.CellINS}")
        shell("bedtools intersect -a {input.INS} -b {input.gnomadINS}  -wa -wb > {output.gnomadINS}")
        shell("bedtools intersect -a {input.INS} -b {input.HGSVCINS}  -wa -wb > {output.HGSVCINS}")
        shell("bedtools intersect -a {input.INS} -b {input.HGDPINS}  -wa -wb > {output.HGDPINS}")
        shell("bedtools intersect -a {input.INV} -b {input.CellINV}  -wa -wb > {output.CellINV}")
        shell("bedtools intersect -a {input.INV} -b {input.gnomadINV}  -wa -wb > {output.gnomadINV}")
        shell("bedtools intersect -a {input.INV} -b {input.HGDPINV}  -wa -wb > {output.HGDPINV}")
        shell("bedtools intersect -a {input.DUP} -b {input.gnomadDUP}  -wa -wb > {output.gnomadDUP}")
        shell("bedtools intersect -a {input.DUP} -b {input.HGDPDUP}  -wa -wb > {output.HGDPDUP}")


rule AsianFilt:
    input:
        CellDEL = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_DEL.txt",
        CellINS = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_INS.txt",
        CellINV = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_INV.txt",
        gnomadDEL = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DEL.txt",
        gnomadINS = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INS.txt",
        gnomadINV = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INV.txt",
        gnomadDUP = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DUP.txt",
        HGSVCDEL = IN_PATH + "/preStudies/population/Asian/overlap/HGSVC2.East_Asian.reformat_DEL.txt",
        HGSVCINS = IN_PATH + "/preStudies/population/Asian/overlap/HGSVC2.East_Asian.reformat_INS.txt",
        HGDPDEL = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_DEL.txt",
        HGDPINS = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_INS.txt",
        HGDPINV = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_INV.txt",
        HGDPDUP = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_DUP.txt",
    output:
        CellDEL = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_DEL_filt.txt",
        CellINS = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_INS_filt.txt",
        CellINV = IN_PATH + "/preStudies/population/Asian/overlap/Cell2019_Table_Asian_INV_filt.txt",
        gnomadDEL = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DEL_filt.txt",
        gnomadINS = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INS_filt.txt",
        gnomadINV = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_INV_filt.txt",
        gnomadDUP = IN_PATH + "/preStudies/population/Asian/overlap/gnomad_v2.1_sv.sites.EAS.reformat.hg38_DUP_filt.txt",
        HGSVCDEL = IN_PATH + "/preStudies/population/Asian/overlap/HGSVC2.East_Asian.reformat_DEL_filt.txt",
        HGSVCINS = IN_PATH + "/preStudies/population/Asian/overlap/HGSVC2.East_Asian.reformat_INS_filt.txt",
        HGDPDEL = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_DEL_filt.txt",
        HGDPINS = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_INS_filt.txt",
        HGDPINV = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_INV_filt.txt",
        HGDPDUP = IN_PATH + "/preStudies/population/Asian/overlap/WGS911_Asia_DUP_filt.txt",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/AsianFilt.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.CellDEL} --ratioThreshold {params.ratioThreshold} --out {output.CellDEL} --method both > {log} 2>&1")
        shell("python {params.bedOverlapRecord} --bed {input.CellINS} --ratioThreshold {params.ratioThreshold} --out {output.CellINS} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.CellINV} --ratioThreshold {params.ratioThreshold} --out {output.CellINV} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.gnomadDEL} --ratioThreshold {params.ratioThreshold} --out {output.gnomadDEL} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.gnomadINS} --ratioThreshold {params.ratioThreshold} --out {output.gnomadINS} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.gnomadINV} --ratioThreshold {params.ratioThreshold} --out {output.gnomadINV} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.gnomadDUP} --ratioThreshold {params.ratioThreshold} --out {output.gnomadDUP} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.HGSVCDEL} --ratioThreshold {params.ratioThreshold} --out {output.HGSVCDEL} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.HGSVCINS} --ratioThreshold {params.ratioThreshold} --out {output.HGSVCINS} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.HGDPDEL} --ratioThreshold {params.ratioThreshold} --out {output.HGDPDEL} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.HGDPINS} --ratioThreshold {params.ratioThreshold} --out {output.HGDPINS} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.HGDPINV} --ratioThreshold {params.ratioThreshold} --out {output.HGDPINV} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.HGDPDUP} --ratioThreshold {params.ratioThreshold} --out {output.HGDPDUP} --method both 2>> {log}")

###############################################################
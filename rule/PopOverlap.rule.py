
################################ overlap with and Cell2019 ###########################################
rule overlapCell:
    input:
        Del = IN_PATH + "/population/bed/Sample_common_SV_DEL.bed",
        Ins = IN_PATH + "/population/bed/Sample_common_SV_INS.bed",
        Inv = IN_PATH + "/population/bed/Sample_common_SV_INV.bed",
    output:
        Del = IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap.bed",
        Ins = IN_PATH + "/population/bed/LRS15/Sample_INS_overlap.bed",
        Inv = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap.bed",
    params:
        Del = IN_PATH + "/preStudies/Cell2019_SV_dechr_DEL.bed",
        Ins = IN_PATH + "/preStudies/Cell2019_SV_dechr_INS_newEnd.bed",
        Inv = IN_PATH + "/preStudies/Cell2019_SV_dechr_INV.bed",
    run:
        shell("bedtools intersect -a {input.Del} -b {params.Del} -wa -wb > {output.Del}")
        shell("bedtools intersect -a {input.Ins} -b {params.Ins} -wa -wb > {output.Ins}")
        shell("bedtools intersect -a {input.Inv} -b {params.Inv} -wa -wb > {output.Inv}")

rule overlapCellFilt:
    input:
        Del = IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap.bed",
        Ins = IN_PATH + "/population/bed/LRS15/Sample_INS_overlap.bed",
        Inv = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap.bed",
    output:
        Del = IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap_filt.bed",
        Ins = IN_PATH + "/population/bed/LRS15/Sample_INS_overlap_filt.bed",
        Inv = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/overlapCellFilt.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")
        shell("python {params.bedOverlapRecord} --bed {input.Inv} --ratioThreshold {params.ratioThreshold} --out {output.Inv} --method both 2>> {log}")
        shell("python {params.bedOverlapRecord} --bed {input.Ins} --ratioThreshold {params.ratioThreshold} --out {output.Ins} --method both 2>> {log}")

        


rule CellOverlapSum:
    input:
        Del = IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap_filt.bed",
        Ins = IN_PATH + "/population/bed/LRS15/Sample_INS_overlap_filt.bed",
        Inv = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap_filt.bed",
    output:
        tag = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/CellOverlapSum.log",
    run:
        BED = ",".join([input.Del, input.Inv, input.Ins])
        shell("python {params.OverlapSVTypeLength}  --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        shell("cat {input.Del} {input.Inv} {input.Ins} | cut -f 4 | sort | uniq > {output.tag}")


rule CellOverlapSumPlot:
    input:
        length = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_chr_length_summary.xls",
    output:
        number = IN_PATH + "/population/bed/LRS15/Sample_common_SV_number.pdf",
        length = IN_PATH + "/population/bed/LRS15/Sample_common_SV_length.pdf",
    params:
        OverlapSVLenChrBar = SCRIPT_DIR + "/OverlapSVLenChrBar.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/CellOverlapSumPlot.log",
    run:
        shell("Rscript {params.OverlapSVLenChrBar} --input {input.length} --pdf1 {output.number} --pdf2 {output.length} --width {params.width} --height {params.height} > {log} 2>&1")
#####################################################################################################










################################### gnomAD overlap ###############################################
rule overlapgnomAD:
    input:
        Del = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
        gnomAD = "/home/wuzhikun/database/gnomAD/nonredundant_200629/gnomad_v2.1_sv.sites.{SVtype}.from_vcf-4.bed",
    output:
        Del = IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap.bed",
    run:
        # annoBED = config["gnomAD"][wildcards.SVtype]
        shell("bedtools intersect -a {input.Del} -b {input.gnomAD} -wa -wb > {output.Del}")




rule gnomADFilt:
    input:
        Del = IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap.bed",
    output:
        Del = IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/gnomADFilt_{SVtype}.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")


rule gnomADOverlapSum:
    input:
        Del = expand(IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        tag = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/gnomAD/Sample_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/gnomADOverlapSum.log",
    run:
        BED = ",".join(input.Del)
        shell("python {params.OverlapSVTypeLength} --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        BED2 = " ".join(input.Del)
        shell("cat {BED2} | cut -f 4 | sort | uniq > {output.tag}")


rule gnomADOverlapSumPlot:
    input:
        length = IN_PATH + "/population/bed/gnomAD/Sample_overlap_filt_chr_length_summary.xls",
    output:
        number = IN_PATH + "/population/bed/gnomAD/Sample_common_overlap_SV_number.pdf",
        length = IN_PATH + "/population/bed/gnomAD/Sample_common_overlap_SV_length.pdf",
    params:
        OverlapSVLenChrBar = SCRIPT_DIR + "/OverlapSVLenChrBar.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/gnomADOverlapSumPlot.log",
    run:
        shell("Rscript {params.OverlapSVLenChrBar} --input {input.length} --pdf1 {output.number} --pdf2 {output.length} --width {params.width} --height {params.height} > {log} 2>&1")
#####################################################################################################


################################ overlap with DGV ##################################

rule overlapDGV:
    input:
        Del = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
        DGV = "/home/wuzhikun/database/DGV/type_200626/GRCh38_hg38_variants_{SVtype}.bed",
    output:
        Del = IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap.bed",
    log:
        IN_PATH + "/log/overlapDGV_{SVtype}.log",
    run:
        shell("bedtools intersect -a {input.Del} -b {input.DGV} -wa -wb > {output.Del} 2>{log}")




rule DGVFilt:
    input:
        Del = IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap.bed",
    output:
        Del = IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/DGVFilt_{SVtype}.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")


rule DGVOverlapSum:
    input:
        Del = expand(IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        tag = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/DGV/Sample_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/DGVOverlapSum.log",
    run:
        BED = ",".join(input.Del)
        shell("python {params.OverlapSVTypeLength} --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        BED2 = " ".join(input.Del)
        shell("cat {BED2} | cut -f 4 | sort | uniq > {output.tag}")

rule DGVOverlapSumPlot:
    input:
        length = IN_PATH + "/population/bed/DGV/Sample_overlap_filt_chr_length_summary.xls",
    output:
        number = IN_PATH + "/population/bed/DGV/Sample_common_overlap_SV_number.pdf",
        length = IN_PATH + "/population/bed/DGV/Sample_common_overlap_SV_length.pdf",
    params:
        OverlapSVLenChrBar = SCRIPT_DIR + "/OverlapSVLenChrBar.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/DGVOverlapSumPlot.log",
    run:
        shell("Rscript {params.OverlapSVLenChrBar} --input {input.length} --pdf1 {output.number} --pdf2 {output.length} --width {params.width} --height {params.height} > {log} 2>&1")

####################################################################################



################################ overlap with dbVar ##################################

rule overlapdbVar:
    input:
        Del = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
        dbVar = "/home/wuzhikun/database/dbVar/hg38_nr/GRCh38.nr_{SVtype}.bed",
    output:
        Del = IN_PATH + "/population/bed/dbVar/Sample_{SVtype}_overlap.bed",
    log:
        IN_PATH + "/log/overlapdbVar_{SVtype}.log",
    run:
        shell("bedtools intersect -a {input.Del} -b {input.dbVar} -wa -wb > {output.Del} 2>{log}")




rule dbVarFilt:
    input:
        Del = IN_PATH + "/population/bed/dbVar/Sample_{SVtype}_overlap.bed",
    output:
        Del = IN_PATH + "/population/bed/dbVar/Sample_{SVtype}_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/dbVarFilt_{SVtype}.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")


rule dbVarOverlapSum:
    input:
        Del = expand(IN_PATH + "/population/bed/dbVar/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        tag = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/dbVar/Sample_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/dbVarOverlapSum.log",
    run:
        BED = ",".join(input.Del)
        shell("python {params.OverlapSVTypeLength} --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        BED2 = " ".join(input.Del)
        shell("cat {BED2} | cut -f 4 | sort | uniq > {output.tag}")

rule dbVarOverlapSumPlot:
    input:
        length = IN_PATH + "/population/bed/dbVar/Sample_overlap_filt_chr_length_summary.xls",
    output:
        number = IN_PATH + "/population/bed/dbVar/Sample_common_overlap_SV_number.pdf",
        length = IN_PATH + "/population/bed/dbVar/Sample_common_overlap_SV_length.pdf",
    params:
        OverlapSVLenChrBar = SCRIPT_DIR + "/OverlapSVLenChrBar.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/dbVarOverlapSumPlot.log",
    run:
        shell("Rscript {params.OverlapSVLenChrBar} --input {input.length} --pdf1 {output.number} --pdf2 {output.length} --width {params.width} --height {params.height} > {log} 2>&1")

####################################################################################



################################ overlap with WGS17795 ##################################

rule overlapdWGS17795:
    input:
        Del = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
        WGS17795 = "/home/wuzhikun/database/WGS17795/nonredundant_200628/Build38_{SVtype}.bed",
    output:
        Del = IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap.bed",
    log:
        IN_PATH + "/log/overlapWGS17795_{SVtype}.log",
    run:
        shell("bedtools intersect -a {input.Del} -b {input.WGS17795} -wa -wb > {output.Del} 2>{log}")




rule WGS17795Filt:
    input:
        Del = IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap.bed",
    output:
        Del = IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/WGS17795Filt_{SVtype}.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")


rule WGS17795OverlapSum:
    input:
        Del = expand(IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        tag = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/WGS17795/Sample_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/WGS17795OverlapSum.log",
    run:
        BED = ",".join(input.Del)
        shell("python {params.OverlapSVTypeLength} --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        BED2 = " ".join(input.Del)
        shell("cat {BED2} | cut -f 4 | sort | uniq > {output.tag}")


rule WGS17795OverlapSumPlot:
    input:
        length = IN_PATH + "/population/bed/WGS17795/Sample_overlap_filt_chr_length_summary.xls",
    output:
        number = IN_PATH + "/population/bed/WGS17795/Sample_common_overlap_SV_number.pdf",
        length = IN_PATH + "/population/bed/WGS17795/Sample_common_overlap_SV_length.pdf",
    params:
        OverlapSVLenChrBar = SCRIPT_DIR + "/OverlapSVLenChrBar.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/WGS17795OverlapSumPlot.log",
    run:
        shell("Rscript {params.OverlapSVLenChrBar} --input {input.length} --pdf1 {output.number} --pdf2 {output.length} --width {params.width} --height {params.height} > {log} 2>&1")

####################################################################################




################################ overlap with WGS911 ##################################

rule overlapdWGS911:
    input:
        Del = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
        WGS911 = "/home/wuzhikun/database/WGS911/nonredundant/combine/{SVtype}_frequency.bed",
    output:
        Del = IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap.bed",
    log:
        IN_PATH + "/log/overlapWGS911_{SVtype}.log",
    run:
        shell("bedtools intersect -a {input.Del} -b {input.WGS911} -wa -wb > {output.Del} 2>{log}")




rule WGS911Filt:
    input:
        Del = IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap.bed",
    output:
        Del = IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/WGS911Filt_{SVtype}.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.Del} --ratioThreshold {params.ratioThreshold} --out {output.Del} --method both > {log} 2>&1")



rule WGS911OverlapSum:
    input:
        Del = expand(IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        tag = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/WGS911/Sample_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/WGS911OverlapSum.log",
    run:
        BED = ",".join(input.Del)
        shell("python {params.OverlapSVTypeLength} --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        BED2 = " ".join(input.Del)
        shell("cat {BED2} | cut -f 4 | sort | uniq > {output.tag}")


rule WGS911OverlapSumPlot:
    input:
        length = IN_PATH + "/population/bed/WGS911/Sample_overlap_filt_chr_length_summary.xls",
    output:
        number = IN_PATH + "/population/bed/WGS911/Sample_common_overlap_SV_number.pdf",
        length = IN_PATH + "/population/bed/WGS911/Sample_common_overlap_SV_length.pdf",
    params:
        OverlapSVLenChrBar = SCRIPT_DIR + "/OverlapSVLenChrBar.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/WGS911OverlapSumPlot.log",
    run:
        shell("Rscript {params.OverlapSVLenChrBar} --input {input.length} --pdf1 {output.number} --pdf2 {output.length} --width {params.width} --height {params.height} > {log} 2>&1")




def tag_overlap_record(tag_file, overlap_file, out_file):
    TagSets = set()
    tag_h = open(tag_file, "r")
    for line in tag_h:
        lines = line.strip().split("\t")
        if lines != []:
            tag = lines[0]
            TagSets.add(tag)
    tag_h.close()
    overlap_h = open(overlap_file, "r")
    out_h = open(out_file, "w")
    for line in overlap_h:
        line = line.strip()
        lines = line.split("\t")
        t = lines[3]
        if t in TagSets:
            out_h.write("%s\n" % line)
    overlap_h.close()
    out_h.close()



rule OverlapCDS:
    input:
        overlap = IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt.bed",
        cds = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS.bed",
    output:
        tag = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap_promoter_all_lof_CDS_tag_{SVtype}_temp.txt",
        record = IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt_cds.bed",
    run:
        shell("cut -f 4 {input.cds} | sort | uniq > {output.tag}")
        tag_overlap_record(output.tag, input.overlap, output.record)

####################################################################################






############################################  summary #########################################
rule overlapSummary:
    input:
        Del = IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap_filt.bed",
        Ins = IN_PATH + "/population/bed/LRS15/Sample_INS_overlap_filt.bed",
        Inv = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap_filt.bed",
        gnomAD = expand(IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        DGV = expand(IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        # dbVar = expand(IN_PATH + "/population/bed/dbVar/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        # WGS17795 = expand(IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
        WGS911 = expand(IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        summary = IN_PATH + "/population/bed/All/Sample_overlap_summary.txt",
    params:
        OverlapSVSummary = SRC_DIR + "/OverlapSVSummary.py",
    log:
        IN_PATH + "/log/overlapSummary.log"
    run:
        # + input.WGS17795
        Files = ",".join([input.Del, input.Ins, input.Inv] + input.gnomAD + input.DGV  + input.WGS911)
        shell("python {params.OverlapSVSummary} --file {Files} --out {output.summary} > {log} 2>&1")


rule overlapSummaryPlot:
    input:
        summary = IN_PATH + "/population/bed/All/Sample_overlap_summary.txt",
    output:
        pdf = IN_PATH + "/population/bed/All/Sample_overlap_summary.pdf",
    params:
        OverlapSummary = SCRIPT_DIR + "/OverlapSummary.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/overlapSummary.log"
    run:
        shell("Rscript {params.OverlapSummary} --input {input.summary} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")



rule knownAndNovelTags:
    input:
        overlap = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap.xls",
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
    output:
        stats = IN_PATH + "/population/Category/Kown_and_novel_tags_stats.txt",
        novel = IN_PATH + "/population/Category/Novel_SV_category_tag.xls",
    params:
        CatNovelStats = SRC_DIR + "/CatNovelStats.py",
    log:
        IN_PATH + "/log/knownAndNovelTags.log"
    run:
        shell("python {params.CatNovelStats} --category {input.category} --overlap {input.overlap} --stats {output.stats} --novel {output.novel} > {log} 2>&1")


rule knownAndNovelTagsPlot:
    input:
        stats = IN_PATH + "/population/Category/Kown_and_novel_tags_stats.txt",
    output:
        pdf = IN_PATH + "/population/Category/Kown_and_novel_tags_stats.pdf",
    params:
        KnownNovelStack = SCRIPT_DIR + "/KnownNovelStack.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/knownAndNovelTagsPlot.log"
    run:
        shell("Rscript {params.KnownNovelStack} --input {input.stats} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")






############################################################################################################





################################## venn plot  ########################################
rule multipleTags:
    input:
        LRS15 = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tags.xls",
        gnomAD = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tags.xls",
        # dbVar = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tags.xls",
        DGV = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tags.xls",
        WGS911 = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tags.xls",
        # WGS17795 = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tags.xls",
    output:
        tag = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap.xls",
        tag1 = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap-1.xls",
    params:
        multipleOverlap = SRC_DIR + "/multipleOverlap.py",
    log:
        IN_PATH + "/log/multipleTags.log"
    run:
        Files = ",".join([input.LRS15, input.gnomAD, input.DGV, input.WGS911])
        shell("python {params.multipleOverlap} --files {Files} --out {output.tag} > {log} 2>&1")
        shell("cat {input.LRS15} {input.gnomAD} {input.DGV} {input.WGS911} | sort | uniq > {output.tag1}")



def tag_type(tag_file):
    TypeCount = collections.Counter()
    tag_h = open(tag_file, "r")
    for line in tag_h:
        lines = line.strip().split("\t")
        if lines != []:
            tag = lines[0]
            SVType = tag.split("-")[-1]
            TypeCount[SVType] += 1
    tag_h.close()
    return TypeCount


def known_novel_SVtype(novel_tag, overlap_tag, out_file):
    novelCount = tag_type(novel_tag)
    knownCount = tag_type(overlap_tag)
    typeList = sorted(list(set(list(novelCount.keys()) + list(knownCount.keys()))))

    out_h = open(out_file, "w")
    out_h.write("Category\tKnownTags\tNovelTags\tKnownRatio\tNovelRatio\n")
    for t in typeList:
        if t in novelCount:
            novel = novelCount[t]
        else:
            novel = 0

        if t in knownCount:
            known = knownCount[t]
        else:
            known = 0

        both = known + novel

        if both != 0:
            knownRatio = known / both
            knownRatio = "%.3f" % knownRatio
            novelRatio = novel / both
            novelRatio = "%.3f" % novelRatio
            out_h.write("%s\t%d\t%d\t%s\t%s\n" % (t, known, novel, knownRatio, novelRatio))
    out_h.close()
    



rule SVtypeKnown:
    input:
        overlap = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap-1.xls",
        allTag = IN_PATH + "/population/Merge/Sample_common_SV_Tag.xls",
    output:
        allTag = IN_PATH + "/population/Merge/Sample_common_SV_Tag_name.xls",
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tags_overlap_novel.xls",
        stat = IN_PATH + "/population/bed/novel/SV_type_known_and_novwl_stats.txt",
    run:
        shell("sed '1d' {input.allTag} | cut -f 1 > {output.allTag}")
        shell("cat {output.allTag} {input.overlap} | sort | uniq -c | sort -k 1nr | sed 's/^      //g' | grep '^1' | cut -f 2 -d ' ' > {output.novel}")
        known_novel_SVtype(output.novel, input.overlap, output.stat)


rule SVtypeKnownPlot:
    input:
        stats = IN_PATH + "/population/bed/novel/SV_type_known_and_novwl_stats.txt",
    output:
        pdf = IN_PATH + "/population/bed/novel/SV_type_known_and_novwl_stats.pdf",
    params:
        KnownNovelStack = SCRIPT_DIR + "/KnownNovelStack.R",
        width = 4,
        height = 4,
    log:
        IN_PATH + "/log/SVtypeKnownPlot.log"
    run:
        shell("Rscript {params.KnownNovelStack} --input {input.stats} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")




rule multipleTagsPlot:
    ### edit the header of the file "Cell2019        DGV     dbVar   gnomAD"
    input:
        tag = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap.xls",
    output:
        png = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap.png",
    params:
        venn_4categories = SCRIPT_DIR + "/venn_5categories.R",
        width = 3000,
        height = 3000,
        cex = 0.8,
        outType = "png",
    log:
        IN_PATH + "/log/multipleTagsPlot.log"
    run:
        shell("Rscript {params.venn_4categories} {input.tag} {output.png}  {params.outType}  {params.height} {params.width} {params.cex} > {log} 2>&1")


rule venn2:
    input:
        tag = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap.xls",
    output:
        tag = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap_select2.xls",
    run:
        shell("cut -f 1-2 {input.tag} > {output.tag}")



rule multipleTagsPlot2:
    input:
        tag = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap_select2.xls",
    output:
        png = IN_PATH + "/population/bed/SV_overlap_filt_tags_overlap_select2.png",
    params:
        venn_4categories = SCRIPT_DIR + "/venn_2categories.R",
        width = 3000,
        height = 3000,
        cex = 0.8,
        outType = "png",
    log:
        IN_PATH + "/log/multipleTagsPlot2.log"
    run:
        shell("Rscript {params.venn_4categories} {input.tag} {output.png}  {params.outType}  {params.height} {params.width} {params.cex} > {log} 2>&1")


###################################################################################



####################################################################################

rule overlapFreq:
    input:
        LRS15 = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tags.xls",
        gnomAD = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tags.xls",
        # dbVar = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tags.xls",
        DGV = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tags.xls",
        # WGS17795 = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tags.xls",
        WGS911 = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tags.xls",
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
    output:
        LRS15 = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tag_freq.xls",
        gnomAD = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tag_freq.xls",
        # dbVar = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tag_freq.xls",
        DGV = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tag_freq.xls",
        # WGS17795 = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tag_freq.xls",
        WGS911 = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tag_freq.xls",
    run:
        target_tag_freq(input.LRS15, input.freq, output.LRS15, "include")
        target_tag_freq(input.gnomAD, input.freq, output.gnomAD, "include")
        target_tag_freq(input.DGV, input.freq, output.DGV, "include")
        # target_tag_freq(input.dbVar, input.freq, output.dbVar, "include")
        # target_tag_freq(input.WGS17795, input.freq, output.WGS17795, "include")
        target_tag_freq(input.WGS911, input.freq, output.WGS911, "include")


rule FreqStat:
    input:
        LRS15 = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tag_freq.xls",
        gnomAD = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tag_freq.xls",
        # dbVar = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tag_freq.xls",
        DGV = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tag_freq.xls",
        # WGS17795 = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tag_freq.xls",
        WGS911 = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tag_freq.xls",
    output:
        LRS15 = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tag_freq_stat.xls",
        gnomAD = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tag_freq_stat.xls",
        # dbVar = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tag_freq_stat.xls",
        DGV = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tag_freq_stat.xls",
        # WGS17795 = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tag_freq_stat.xls",
        WGS911 = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tag_freq_stat.xls",
    params:
        TypeAlleleFrequency = SRC_DIR + "/TypeAlleleFrequency.py",
    log:
        IN_PATH + "/log/FreqStat.log", 
    run:
        shell("python {params.TypeAlleleFrequency} --input {input.LRS15} --out {output.LRS15} > {log} 2>&1")
        shell("python {params.TypeAlleleFrequency} --input {input.gnomAD} --out {output.gnomAD} >> {log} 2>&1")
        # shell("python {params.TypeAlleleFrequency} --input {input.dbVar} --out {output.dbVar} >> {log} 2>&1")
        shell("python {params.TypeAlleleFrequency} --input {input.DGV} --out {output.DGV} >> {log} 2>&1")
        # shell("python {params.TypeAlleleFrequency} --input {input.WGS17795} --out {output.WGS17795} >> {log} 2>&1")
        shell("python {params.TypeAlleleFrequency} --input {input.WGS911} --out {output.WGS911} >> {log} 2>&1")




rule FreqSummary:
    input:
        LRS15 = IN_PATH + "/population/bed/LRS15/SV_overlap_filt_tag_freq_stat.xls",
        gnomAD = IN_PATH + "/population/bed/gnomAD/SV_overlap_filt_tag_freq_stat.xls",
        # dbVar = IN_PATH + "/population/bed/dbVar/SV_overlap_filt_tag_freq_stat.xls",
        DGV = IN_PATH + "/population/bed/DGV/SV_overlap_filt_tag_freq_stat.xls",
        # WGS17795 = IN_PATH + "/population/bed/WGS17795/SV_overlap_filt_tag_freq_stat.xls",
        WGS911 = IN_PATH + "/population/bed/WGS911/SV_overlap_filt_tag_freq_stat.xls",
    output:
        summary = IN_PATH + "/population/bed/SV_overlap_filt_tag_freq_stat_summary.xls",
    run:
        Files = ",".join([input.LRS15, input.gnomAD,  input.DGV,  input.WGS911])
        multiple_file_SV_summary(Files, output.summary)



rule FreqSummaryPlot:
    input:
        summary = IN_PATH + "/population/bed/SV_overlap_filt_tag_freq_stat_summary.xls",
    output:
        pdf = IN_PATH + "/population/bed/SV_overlap_filt_tag_freq_stat_summary.pdf",
    params:
        TypeFreqBarSingle = SCRIPT_DIR + "/TypeFreqBarSingle.R",
        width = 6,
        height = 4,
    log:
        IN_PATH + "/log/FreqSummaryPlot.log", 
    run:
        shell("Rscript {params.TypeFreqBarSingle} --input {input.summary} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
######################################################################################









########################################### all overlap and novel sV ##########################################
rule database:
    input:
        gnomAD = "/home/wuzhikun/database/gnomAD/nonredundant_200629/gnomad_v2.1_sv.sites.{SVtype}.from_vcf-4.bed",
        DGV = "/home/wuzhikun/database/DGV/type_200626/GRCh38_hg38_variants_{SVtype}.bed",
        # WGS17795 = "/home/wuzhikun/database/WGS17795/nonredundant_200628/Build38_{SVtype}.bed",
        WGS911 = "/home/wuzhikun/database/WGS911/nonredundant/combine/{SVtype}_frequency.bed",
    output:
        SV = temp(IN_PATH + "/population/bed/database/database_SV_{SVtype}_temp.bed"),
    run:
        shell("cat {input.DGV} {input.gnomAD} {input.WGS911} > {output.SV}")


rule database2:
    input:
        Del = IN_PATH + "/population/bed/database/database_SV_DEL_temp.bed",
        Ins = IN_PATH + "/population/bed/database/database_SV_INS_temp.bed",
        Inv = IN_PATH + "/population/bed/database/database_SV_INV_temp.bed",
        Dup = IN_PATH + "/population/bed/database/database_SV_DUP_temp.bed",
        CellDel = IN_PATH + "/preStudies/Cell2019_SV_dechr_DEL.bed",
        CellIns = IN_PATH + "/preStudies/Cell2019_SV_dechr_INS_newEnd.bed",
        CellInv = IN_PATH + "/preStudies/Cell2019_SV_dechr_INV.bed",
    output:
        Del = IN_PATH + "/population/bed/database/database_SV_DEL.bed",
        Ins = IN_PATH + "/population/bed/database/database_SV_INS.bed",
        Inv = IN_PATH + "/population/bed/database/database_SV_INV.bed",
        Dup = IN_PATH + "/population/bed/database/database_SV_DUP.bed",
    run:
        shell("cat {input.CellDel} {input.Del} | cut -f 1-4 | sort -k 1,1n -k 2,2n > {output.Del}")
        shell("cat {input.CellIns} {input.Ins} | cut -f 1-4 | sort -k 1,1n -k 2,2n > {output.Ins}")
        shell("cat {input.CellInv} {input.Inv} | cut -f 1-4 | sort -k 1,1n -k 2,2n  > {output.Inv}")
        shell("cat {input.Dup} | cut -f 1-4 | sort -k 1,1n -k 2,2n >  {output.Dup}")


rule dbOverlap:
    input:
        SV = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
        db = IN_PATH + "/population/bed/database/database_SV_{SVtype}.bed",
    output:
        db = IN_PATH + "/population/bed/database/Sample_{SVtype}_overlap.bed",
    log:
        IN_PATH + "/log/dbOverlap_{SVtype}.log",
    run:
        shell("bedtools intersect -a {input.SV} -b {input.db} -wa -wb > {output.db} 2>{log}")

rule dbOverlapFilt:
    input:
        db = IN_PATH + "/population/bed/database/Sample_{SVtype}_overlap.bed",
    output:
        db = IN_PATH + "/population/bed/database/Sample_{SVtype}_overlap_filt.bed",
    params:
        bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
        ratioThreshold = 0.5,
    log:
        IN_PATH + "/log/dbOverlapFilt_{SVtype}.log",
    run:
        shell("python {params.bedOverlapRecord} --bed {input.db} --ratioThreshold {params.ratioThreshold} --out {output.db} --method both > {log} 2>&1")




rule allOverlapSum:
    input:
        overlap = expand(IN_PATH + "/population/bed/database/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        tag = IN_PATH + "/population/bed/database/SV_overlap_filt_tags.xls",
        lengthSummary = IN_PATH + "/population/bed/database/SV_overlap_filt_chr_length_summary.xls",
    params:
        OverlapSVTypeLength = SRC_DIR + "/OverlapSVTypeLength.py",
    log:
        IN_PATH + "/log/allOverlapSum.log",
    run:
        BED = ",".join(input.overlap)
        shell("python {params.OverlapSVTypeLength}  --bed {BED} --out {output.lengthSummary} > {log} 2>&1")
        Files = " ".join(input.overlap)
        shell("cat {Files} | cut -f 4 | sort | uniq > {output.tag}")


rule dataFreq:
    input:
        overlap = IN_PATH + "/population/bed/database/SV_overlap_filt_tags.xls",
        freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
    output:
        overlap = IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq.xls",
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq.xls",
        novelTag = IN_PATH + "/population/bed/novel/SV_overlap_filt_tags.xls",
    run:
        target_tag_freq(input.overlap, input.freq, output.overlap, "include")
        target_tag_freq(input.overlap, input.freq, output.novel, "exclude")
        shell("cut -f 1 {output.novel} > {output.novelTag}")





rule dataFreqStat:
    input:
        overlap = IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq.xls",
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq.xls",
    output:
        overlap = IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq_stat.xls",
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq_stat.xls",
    params:
        TypeAlleleFrequency = SRC_DIR + "/TypeAlleleFrequency.py",
    log:
        IN_PATH + "/log/dataFreqStat.log", 
    run:
        shell("python {params.TypeAlleleFrequency} --input {input.overlap} --out {output.overlap} > {log} 2>&1")
        shell("python {params.TypeAlleleFrequency} --input {input.novel} --out {output.novel} >> {log} 2>&1")



rule NovelFrequencyPlot:
    input:
        overlap = IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq_stat.xls",
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq_stat.xls",
    output:
        overlap = IN_PATH + "/population/bed/database/SV_overlap_filt_tag_freq_stat.pdf",
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq_stat.pdf",
    params:
        TypeFreqBarSingle = SCRIPT_DIR + "/TypeFreqBarSingle-1.R",
        width = 6, 
        height = 4,
    log:
        IN_PATH + "/log/NovelFrequencyPlot.log", 
    run:
        shell("Rscript {params.TypeFreqBarSingle} --input {input.overlap} --pdf {output.overlap} --width {params.width} --height {params.height} > {log} 2>&1")
        shell("Rscript {params.TypeFreqBarSingle} --input {input.novel} --pdf {output.novel} --width {params.width} --height {params.height} >> {log} 2>&1")


###############################################################################################






####################################################### novel gene annotation #############################
rule Targettag:
    input:
        novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq.xls",
    output:
        common = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common.xls",
        rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare.xls",
    run:
        cmd1 = """awk '{if($4>0.05){print $0}}' %s | cut -f 1 | grep -v Tag > %s""" % (input.novel, output.common)
        print(cmd1)
        os.system(cmd1)
        cmd2 = """awk '{if($4<0.01){print $0}}' %s | cut -f 1 | grep -v Tag > %s""" % (input.novel, output.rare)
        print(cmd2)
        os.system(cmd2)







rule tagGene:
    input:
        common = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common.xls",
        rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare.xls",
        annotation = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
    output:
        common = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_gene.xls",
        rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_gene.xls",
        allGene = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_gene.xls",
        commonRecord = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_record.xls",
        rareRecord = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_record.xls",
        allRecord = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_record.xls",
    params:
        TagGeneFeature = SRC_DIR + "/TagGeneFeature.py",
        CodingGene = config["CodingGene"],
    log:
        IN_PATH + "/log/tagGene.log", 
    run:
        ### select coding genes
        shell("python {params.TagGeneFeature} --coding {params.CodingGene} --annotation {input.annotation} --tag {input.common} --out {output.common} --record {output.commonRecord} --iscoding true > {log} 2>&1")
        shell("python {params.TagGeneFeature} --coding {params.CodingGene} --annotation {input.annotation} --tag {input.rare} --out {output.rare} --record {output.rareRecord} --iscoding true >> {log} 2>&1")
        shell("cat {output.common} {output.rare} | sort | uniq > {output.allGene}")
        shell("cat {output.commonRecord} {output.rareRecord} > {output.allRecord}")




# rule DELCDS:
#     input:
#         commonRecord = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_record.xls",
#     output:
#         DEL = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_DEL_CDS_gene.xls",
#         INS = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_INS_CDS_gene.xls",
#     run:
#         shell("grep  DEL  {input.commonRecord} | grep CDS | cut -f 8 | sort | uniq > {output.DEL}")
#         shell("grep  INS  {input.commonRecord} | grep CDS | cut -f 8 | sort | uniq > {output.INS}")


# rule DELCDSRich:
#     input:
#         DEL = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_DEL_CDS_gene.xls",
#         INS = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_INS_CDS_gene.xls",
#     output:
#         DEL = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_DEL_INS_CDS_gene.xls",
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/DELCDS/KEGG_2019_Human..enrichr.reports.txt", 
#         gwas = IN_PATH + "/population/bed/novel/Enrichment/DELCDS/GWAS_Catalog_2019..enrichr.reports.txt",
#     params:
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/DELCDS",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/RareEnrich.log",
#     run:
#         shell("cat {input.DEL} {input.INS} > {output.DEL}")
#         shell("python {params.geneGSEA} --annotation {output.DEL} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")




# rule NovelEnrichGeneGeno:
#     input:
#         anno = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_all_record.xls",
#         geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
#         omim_expand = IN_PATH + "/population/bed/novel/Enrichment/DELCDS/OMIM_Expanded..enrichr.reports.txt",
#         gwas = IN_PATH + "/population/bed/novel/Enrichment/DELCDS/GWAS_Catalog_2019..enrichr.reports.txt",
#     output:
#         gene_expand = IN_PATH + "/population/bed/novel/Enrichment/DELCDS/OMIM_Expand_gene_tag_genotypes.txt", 
#         gwas = IN_PATH + "/population/bed/novel/Enrichment/DELCDS/GWAS_gene_tag_genotypes.txt", 
#     log:
#         IN_PATH + "/log/NovelEnrichGeneGeno.log",
#     params:
#         DiseaseGeneTagGeno = SRC_DIR + "/DiseaseGeneTagGeno.py",
#     run:
#         shell("python {params.DiseaseGeneTagGeno} --enrich {input.omim_expand} --annotation {input.anno} --genotype {input.geno} --out {output.gene_expand} 2>>{log}")
#         shell("python {params.DiseaseGeneTagGeno} --enrich {input.gwas} --annotation {input.anno} --genotype {input.geno} --out {output.gwas} 2>>{log}")



#################################################################################################################







######################################################## novel rare enrichment ###############################################


# rule RareEnrich:
#     input:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_gene.xls",
#     output:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_gene_cds_temp.xls",
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/KEGG_2019_Human..enrichr.reports.txt", 
#         gwas = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/GWAS_Catalog_2019..enrichr.reports.txt",
#     params:
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/rareCDS",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/RareEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("grep CDS {input.rare} | cut -f 1 > {output.rare}")
#         shell("python {params.geneGSEA} --annotation {output.rare} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")



# rule EnriPlot:
#     input:
#         gwas = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/GWAS_Catalog_2019..enrichr.reports.txt",
#     output:
#         pdf = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/GWAS_Catalog_2019..enrichr.reports.pdf",
#     params:
#         enrichBarPlot = SCRIPT_DIR + "/enrichBarPlot.R",
#         width = 8,
#         height = 4,
#         termNum = 12,
#     log:
#         IN_PATH + "/log/EnriPlot.log",
#     run:
#         shell("Rscript {params.enrichBarPlot} --input {input.gwas} --pdf {output.pdf} --width {params.width} --height {params.height} --termNum {params.termNum} > {log} 2>&1")



# rule RareEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/",
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/rareCDS/Plots/",
#         width = 8,
#         height = 4,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)









# rule RareEnrichCDS:
#     input:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_gene.xls",
#     output:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_gene_temp.xls",
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/rare/KEGG_2019_Human..enrichr.reports.txt", 
#         gwas = IN_PATH + "/population/bed/novel/Enrichment/rare/GWAS_Catalog_2019..enrichr.reports.txt",
#     params:
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/rare",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/RareEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("cut -f 1 {input.rare} > {output.rare}")
#         shell("python {params.geneGSEA} --annotation {output.rare} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")




# rule RareEnrichPlotCDS:
#     input:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/rare/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/rare/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/bed/novel/Enrichment/rare/",
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/rare/Plots/",
#         width = 8,
#         height = 4,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)








# rule NovelDELCDS:
#     input:
#         novel = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_freq.xls",
#     output:
#         common = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_DEL.xls",
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_rare_DEL.xls",
#     run:
#         cmd1 = """awk '{if($4>0.05){print $0}}' %s | cut -f 1 | grep -v INV | grep -v DUP | grep -v Tag > %s""" % (input.novel, output.common)
#         print(cmd1)
#         os.system(cmd1)
#         cmd2 = """awk '{if($4<0.01){print $0}}' %s | cut -f 1 | grep -v INV | grep -v DUP | grep -v Tag > %s""" % (input.novel, output.rare)
#         print(cmd2)
#         os.system(cmd2)
##############################################################################################################



######################################### enrich gene record #########################################
# rule EnrichRecord:
#     input:
#         anno = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_record.xls",
#         enrich = IN_PATH + "/population/bed/novel/Enrichment/rare/GWAS_Catalog_2019..enrichr.reports.txt", 
#     output:
#         anno = IN_PATH + "/population/bed/novel/Enrichment/rare/GWAS_Catalog_2019..enrichr.reports_gene_record.txt",
#     params:
#         EnrichGeneTag = SRC_DIR + "/EnrichGeneTag.py",
#         term = " 'Post bronchodilator FEV1/FVC ratio' ",
#     log:
#         IN_PATH + "/log/EnrichRecord.log",
#     run:
#         shell("python {params.EnrichGeneTag} --annotation {input.anno} --enrichment {input.enrich} --out {output.anno} --term {params.term} > {log} 2>&1")
#########################################################################################################




####################################### novel common gene enrichment #####################################
# rule CommonEnrich:
#     input:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_gene.xls",
#     output:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_gene_temp.xls",
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/common/KEGG_2019_Human..enrichr.reports.txt", 
#     params:
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/common",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/CommonEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("cut -f 1 {input.rare} > {output.rare}")
#         shell("python {params.geneGSEA} --annotation {output.rare} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


# rule CommonEnrichPlot:
#     input:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/common/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/common/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/bed/novel/Enrichment/common/",
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/common/Plots/",
#         width = 8,
#         height = 4,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)




# rule CommonEnrichCDS:
#     input:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_gene.xls",
#     output:
#         rare = IN_PATH + "/population/bed/novel/SV_overlap_filt_tag_common_gene_cds_temp.xls",
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/commonCDS/KEGG_2019_Human..enrichr.reports.txt", 
#     params:
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/commonCDS",
#         geneGSEA = SRC_DIR + "/geneGSEA.py",
#         EnrichLibrary = config["EnrichLibrary"],
#         number = 1,
#     log:
#         IN_PATH + "/log/CommonEnrich.log",
#     run:
#         ### mu01
#         ### --number {params.number} 
#         shell("grep CDS {input.rare} | cut -f 1 > {output.rare}")
#         shell("python {params.geneGSEA} --annotation {output.rare} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


# rule CommonEnrichPlotCDS:
#     input:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/commonCDS/KEGG_2019_Human..enrichr.reports.txt",
#     output:
#         kegg = IN_PATH + "/population/bed/novel/Enrichment/commonCDS/Plots/KEGG_2019_Human..enrichr.reports.pdf",
#     params:
#         geneEnrichment = SCRIPT_DIR + "/geneEnrichment.R",
#         indir = IN_PATH + "/population/bed/novel/Enrichment/commonCDS/",
#         outdir = IN_PATH + "/population/bed/novel/Enrichment/commonCDS/Plots/",
#         width = 8,
#         height = 4,
#         selectNum = 12,
#     log:
#         IN_PATH + "/log/CommonEnrichPlot.log",
#     run:
#         # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, input.kegg, output.kegg, params.selectNum, params.width, params.height, log)
#         # print(cmd)
#         # os.system(cmd)
#         files = os.listdir(params.indir)
#         for f in files:
#             if f.endswith("enrichr.reports.txt"):
#                 fname = params.indir + f
#                 pdf = params.outdir + f.rstrip("txt") + "pdf"
#                 cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --selectNum %s --width %s --height %s > %s 2>&1" % (params.geneEnrichment, fname, pdf, params.selectNum, params.width, params.height, log)
#                 os.system(cmd)

######################################################################################################

































# ############################  Category for overlap to Cell2019 ######################
rule CatStats:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        Del = IN_PATH + "/population/bed/LRS15/Sample_DEL_overlap_filt.bed",
        Ins = IN_PATH + "/population/bed/LRS15/Sample_INS_overlap_filt.bed",
        Inv = IN_PATH + "/population/bed/LRS15/Sample_INV_overlap_filt.bed",
    output:
        stat = IN_PATH + "/population/bed/LRS15/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/categoryTag2.log"
    run:
        BEDS = ",".join([input.Del, input.Ins, input.Inv])
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")


rule CatStatsPlot:
    input:
        stat = IN_PATH + "/population/bed/LRS15/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/LRS15/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/CatStatsPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")





rule gnomADCat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        bed = expand(IN_PATH + "/population/bed/gnomAD/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        stat = IN_PATH + "/population/bed/gnomAD/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/gnomADCat.log"
    run:
        BEDS = ",".join(input.bed)
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")



rule gnomADCatPlot:
    input:
        stat = IN_PATH + "/population/bed/gnomAD/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/gnomAD/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/gnomADCatPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")




rule DGVCat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        bed = expand(IN_PATH + "/population/bed/DGV/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        stat = IN_PATH + "/population/bed/DGV/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/DGVCat.log"
    run:
        BEDS = ",".join(input.bed)
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")



rule DGVCatPlot:
    input:
        stat = IN_PATH + "/population/bed/DGV/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/DGV/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/DGVCatPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")







rule dbVarCat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        bed = expand(IN_PATH + "/population/bed/dbVar/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        stat = IN_PATH + "/population/bed/dbVar/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/dbVarCat.log"
    run:
        BEDS = ",".join(input.bed)
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")



rule dbVarCatPlot:
    input:
        stat = IN_PATH + "/population/bed/dbVar/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/dbVar/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/dbVarCatPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")


rule WGS911Cat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        bed = expand(IN_PATH + "/population/bed/WGS911/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        stat = IN_PATH + "/population/bed/WGS911/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/WGS911Cat.log"
    run:
        BEDS = ",".join(input.bed)
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")



rule WGS911CatPlot:
    input:
        stat = IN_PATH + "/population/bed/WGS911/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/WGS911/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/WGS911CatPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")





rule WGS17795Cat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        bed = expand(IN_PATH + "/population/bed/WGS17795/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        stat = IN_PATH + "/population/bed/WGS17795/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/WGS17795Cat.log"
    run:
        BEDS = ",".join(input.bed)
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")



rule WGS17795CatPlot:
    input:
        stat = IN_PATH + "/population/bed/WGS17795/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/WGS17795/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/WGS17795CatPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
# #############################################################################################







# ##################################### all overlap summary ########################################
rule AllCat:
    input:
        category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
        bed = expand(IN_PATH + "/population/bed/database/Sample_{SVtype}_overlap_filt.bed", SVtype=SVTYPES),
    output:
        stat = IN_PATH + "/population/bed/database/Sample_SV_overlap_category_stats.xls",
    params:
        OverlapCategory = SRC_DIR + "/OverlapCategory.py",
    log:
        IN_PATH + "/log/AllCat.log"
    run:
        BEDS = ",".join(input.bed)
        shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")



rule AllCatPlot:
    input:
        stat = IN_PATH + "/population/bed/database/Sample_SV_overlap_category_stats.xls",
    output:
        pdf = IN_PATH + "/population/bed/database/Sample_SV_overlap_category_stats.pdf",
    params:
        OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/AllCatPlot.log"
    run:
        shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")










# rule NovelSV:
#     input:
#         SV = IN_PATH + "/population/bed/Sample_common_SV_{SVtype}.bed",
#         overlap = IN_PATH + "/population/bed/database/Sample_{SVtype}_overlap_filt.bed",
#     output:
#         temp = temp(IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag_temp.txt"),
#         novel =  IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag.bed",
#     run:
#         shell("cut -f 1-4 {input.overlap} | sort | uniq > {output.temp}")
#         shell("cat {input.SV} {output.temp} | cut -f 1-4 | sort | uniq -c | sort -k 1nr |  cut -f 7-8 -d ' ' | grep '^1' | cut -f 2 -d ' ' > {output.novel}")


# rule NovelSV2:
#     input:
#         category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
#         novel =  expand(IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag.bed", SVtype=SVTYPES),
#     output:
#         stat = IN_PATH + "/population/bed/database/Sample_novel_SV_category_stats.xls",
#     params:
#         OverlapCategory = SRC_DIR + "/OverlapCategory.py",
#     log:
#         IN_PATH + "/log/NovelSV2.log"
#     run:
#         BEDS = ",".join(input.novel)
#         shell("python {params.OverlapCategory} --category {input.category} --bed {BEDS} --out {output.stat} > {log} 2>&1")

# rule NovelCatPlot:
#     input:
#         stat = IN_PATH + "/population/bed/database/Sample_novel_SV_category_stats.xls",
#     output:
#         pdf = IN_PATH + "/population/bed/database/Sample_novel_SV_category_stats.pdf",
#     params:
#         OverlapCategory = SCRIPT_DIR + "/OverlapCategory.R",
#         width = 5,
#         height = 4,
#     log:
#         IN_PATH + "/log/NovelCatPlot.log"
#     run:
#         shell("Rscript {params.OverlapCategory} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")




# rule NovelAnno:
#     input:
#         novel =  expand(IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag.bed", SVtype=SVTYPES),
#         category = IN_PATH + "/population/Category/Sample_SV_category_tag.xls",
#     output:
#         novel = IN_PATH + "/population/bed/database/Sample_all_novel_tag.bed",
#         overlap = IN_PATH + "/population/bed/database/Sample_all_novel_tag_overlap.bed",
#         gene = IN_PATH + "/population/bed/database/Sample_all_novel_tag_anno_gene.txt",
#     params:
#         genePred = config["genePred"],
#         annoCategoryRegion = SRC_DIR + "/annoCategoryRegion.py",
#         targetCategory = "major,common",
#         targetRegion = "CDS",
#     log:
#         IN_PATH + "/log/NovelAnno.log",
#     run:
#         NOVEL = " ".join(input.novel)
#         shell("cat {NOVEL} > {output.novel}")
#         shell("bedtools intersect -a {output.novel} -b {params.genePred} -wa -wb > {output.overlap} 2>{log}")
#         shell("python {params.annoCategoryRegion} --category  {input.category} --annotation {output.overlap} --out {output.gene} --targetCategory {params.targetCategory}  --targetRegion {params.targetRegion} > {log} 2>&1")



# rule GeneDisease:
#     input:
#         gene = IN_PATH + "/population/bed/database/Sample_all_novel_tag_anno_gene.txt",
#     output:
#         disease = IN_PATH + "/population/bed/database/Sample_all_novel_tag_anno_gene_disease.txt",
#     params:
#         geneCardDisease = config["geneCardDisease"],
#     run:
#         cmd = """cat %s  | while read marker ; do sed -n '/'"$marker"'\t.*/p' %s ; done > %s """ % (input.gene, params.geneCardDisease, output.disease)
#         os.system(cmd)





# rule NovelTagFreq:
#     input:
#         bed = IN_PATH + "/population/bed/database/Sample_all_novel_tag.bed",
#         freq = IN_PATH + "/population/genotype/Sample_SV_genotype_frequency.txt",
#     output:
#         freq = IN_PATH + "/population/bed/database/Sample_all_novel_tag_freq.txt",
#     params:
#         extractTag_freq = SRC_DIR + "/extractTag_freq.py",
#     log:
#         IN_PATH + "/log/NovelTagFreq.log",
#     run:
#         shell("python {params.extractTag_freq} --freq {input.freq} --bed {input.bed} --out {output.freq} > {log} 2>&1")



###############################################################################################



#################################### Novel Seq ##################################
# rule NovelTagSeq:
#     input:
#         tag = IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag.bed",
#         seq = IN_PATH + "/population/Repeat/Sample_common_SV_tag_{SVtype}.fasta",
#     output:
#         seq = IN_PATH + "/population/bed/database/Sample_{SVtype}_novel_tag_seq.fasta",
#     run:
#         Tag_fasta(input.tag, input.seq, output.seq)




# def Tag_fasta(bed_file, seq_fasta, out_file):
#     from tinyfasta import FastaParser
#     TagSeq = {}
#     for record in FastaParser(seq_fasta):
#         desc = str(record.description)
#         desc = desc.lstrip(">")
#         seq = str(record.sequence)
#         TagSeq[desc] = seq
#     bed_h = open(bed_file, "r")
#     out_h = open(out_file, "w")
#     for line in bed_h:
#         lines = line.strip()
#         tag = lines[3]
#         if tag in TagSeq:
#             seq = TagSeq[tag]
#             out_h.write(">%s\n%s\n" % (tag, seq))
#     bed_h.close()
#     out_h.close()

################################################################################
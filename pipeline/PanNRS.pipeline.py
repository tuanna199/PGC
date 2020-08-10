import os
import yaml


OUT_PATH = config["OUT_PATH"]
SAMPLES = config["SAMPLES"]
QUAST_REF = config["QUAST_REF"]
QUAST_GFF = config["QUAST_GFF"]
PROTEIN_CODING_GFF = config["PROTEIN_CODING_GFF"]
MINIMAP2_MMI = config["MINIMAP2_MMI"]
DNA_BRNN_MODEL = config["DNA_BRNN_MODEL"]
KAIJU_NODES = config["KAIJU_NODES"]
KAIJU_DATABASE = config["KAIJU_DATABASE"]
BLAST_DB = config["BLAST_DB"]
DIAMOND_DB = config["DIAMOND_DB"]
DIAMOND_NODES = config["DIAMOND_NODES"]
DIAMOND_MAPS = config["DIAMOND_MAPS"]
CHORDATE_TAXID = config["CHORDATE_TAXID"]
MINCOV = config["MINCOV"]
MINDIV = config["MINDIV"]
HX1_REF = config["HX1_REF"]
HUPAN = config["HUPAN"]
SRC_PATH = config["SRC_PATH"]
THREADS = config["THREADS"]



rule all:
    input:
        expand(OUT_PATH + "/quast/{sample}_noalt/report.txt", sample=SAMPLES),
        expand(OUT_PATH + "/one_NRS/{sample}_NRS_one.fa", sample=SAMPLES),
        expand(OUT_PATH + "/quast/{sample}_protein_fraction.txt", sample=SAMPLES),
        expand(OUT_PATH + "/one_NRS/{sample}_NRS_one.paf", sample=SAMPLES),
        expand(OUT_PATH + "/two_NRS/{sample}_NRS_two.fa", sample=SAMPLES),
        OUT_PATH + "/dna-brnn/merge_NRS_two.fa.dna-brnn.bed",
        OUT_PATH + "/dna-brnn/sorted_merge_bed.bed",
        expand(OUT_PATH + "/dna-brnn/{sample}_NRS_nonrep.fa", sample=SAMPLES),
        OUT_PATH + "/kaiju/merge_NRS_nonrep.fa",
        OUT_PATH + "/blast/blast_taxid.tab",
        OUT_PATH + "/diamond/diamond_taxid.out",
        expand(OUT_PATH + "/non_contam/{sample}_NRS_nonrep_nc.fa", sample=SAMPLES),
        expand(OUT_PATH + "/anchor/{sample}_flank_left.fa", sample=SAMPLES),
        expand(OUT_PATH + "/anchor/{sample}_flank_left.paf", sample=SAMPLES),
        expand(OUT_PATH + "/anchor/{sample}_flank_anchor.bed", sample=SAMPLES),
        expand(OUT_PATH + "/non_contam/{sample}_NRS_nonrep_nc.bed", sample=SAMPLES),
        OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.fa",
        OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.can_anchor.fa",
        OUT_PATH + "/cdhit/two_end_anchor.anhit.fa",
        OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.paf",
        OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.minihit.fa",
        OUT_PATH + "/cdhit/combine_2type.minihit.fa",
        OUT_PATH + "/cdhit/final_3type.NRS.fa",
        OUT_PATH + "/cdhit/final_3type.NRS.unplaced.fa",
        expand(OUT_PATH + "/pav/{sample}_pav.pav", sample=SAMPLES),
        OUT_PATH + "/other_human_genome/HX1.paf",
        OUT_PATH + "/other_human_genome/HX1_pav.pav",
        OUT_PATH + "/other_pan_genome/cpg_275Chinese.paf",
        OUT_PATH + "/other_pan_genome/275Chinese_pav.pav"





# Using QUAST to extract non-reference sequences
rule quast_first_extract:
    input:
        polished_genome = OUT_PATH + "/Assembly/Fasta/{sample}_assembly_polish.fasta.gz"
    output:
        quast_report = OUT_PATH + "/quast/{sample}_noalt/report.txt",
        one_NRS = OUT_PATH + "/one_NRS/{sample}_NRS_one.fa",
        one_gff = OUT_PATH + "/one_NRS/{sample}_NRS_one.gff"
    params:
        ref = QUAST_REF,
        gff = QUAST_GFF,
        quast_dir = OUT_PATH + "/quast/{sample}_noalt",
        script = SRC_PATH + "/quast_unmap_stat.py",
        unalign_info = OUT_PATH + "/quast/{sample}_noalt/contigs_reports/contigs_report_{sample}_assembly_polish.unaligned.info",
    threads:
        THREADS
    log:
        OUT_PATH + "/quast/{sample}_quast.log"
    run:
        cmd1 = "source activate pan_genome && quast --no-html --no-snps -o %s -r %s -g %s -t %s %s > %s 2>&1" % (params.quast_dir, params.ref, params.gff, threads, input.polished_genome, log)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "source activate base && python %s -quast_unalign_info %s -assembly_fasta %s -out_NRS %s -out_gff_NRS %s -min_contigs_len 10000" % (params.script, params.unalign_info, input.polished_genome, output.one_NRS, output.one_gff)
        print(cmd2)
        os.system(cmd2)


# Assessing the completeness of protein coding genes
rule assess_gene:
    input:
        feature_from_quast = OUT_PATH + "/quast/{sample}_noalt/genome_stats/{sample}_assembly_polish_genomic_features_any.txt"
    output:
        OUT_PATH + "/quast/{sample}_protein_fraction.txt"
    params:
        gff = PROTEIN_CODING_GFF,
        script = SRC_PATH + "/protein_coding_evaluate_quast.py"
    run:
        cmd1 = "source activate base && python %s -q %s -gff %s > %s" % (params.script, input.feature_from_quast, params.gff, output)
        print(cmd1)
        os.system(cmd1)


# Align NRS against the reference genome again and filter
rule minimap2_align:
    input:
        one_NRS = OUT_PATH + "/one_NRS/{sample}_NRS_one.fa"
    output:
        one_paf = OUT_PATH + "/one_NRS/{sample}_NRS_one.paf",
        two_NRS = OUT_PATH + "/two_NRS/{sample}_NRS_two.fa"
    params:
        ref_mmi = MINIMAP2_MMI,
        script = SRC_PATH + "/stat_realign_to_ref.py"
    threads:
        THREADS
    run:
        cmd1 = "minimap2 -DP -t %s %s %s > %s" % (threads, params.ref_mmi, input.one_NRS, output.one_paf)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "python %s -NRS_one %s -in_paf %s -out_fasta %s" % (params.script, input.one_NRS, output.one_paf, output.two_NRS)
        print(cmd2)
        os.system(cmd2)


def cat_func(OUT_PATH, SAMPLES):
    temp_list = []
    for i in SAMPLES:
        temp_list.append(OUT_PATH + "/two_NRS/" + i + "_NRS_two.fa")
    return temp_list


# Remove Alpha and HSat2,3 repeat sequences
rule repeatmasker_and_dna_brnn:
    input:
        all_fa = unpack(cat_func(OUT_PATH, SAMPLES))
    output:
        merge_fa = OUT_PATH + "/repeatmasker/merge_NRS_two.fa",
        alpha_hsat23_bed = OUT_PATH + "/repeatmasker/masker_all/merge_NRS_two.fa.bed",
        dna_brnn_out = OUT_PATH + "/dna-brnn/merge_NRS_two.fa.dna-brnn.bed",
        combine_bed = OUT_PATH + "/dna-brnn/combine_bed.bed",
        sorted_merge_bed = OUT_PATH + "/dna-brnn/sorted_merge_bed.bed"
    params:
        repeatmasker_dir = OUT_PATH + "/repeatmasker/masker_all",
        script = SRC_PATH + "/repeatmasker_extract_alpha_hsat23.py",
        repeatmasker_out = OUT_PATH + "/repeatmasker/masker_all/merge_NRS_two.fa.out",
        dna_brnn_model = DNA_BRNN_MODEL
    threads:
        THREADS
    log:
        OUT_PATH + "/repeatmasker/repma_all.log"
    run:
        cmd1 = "cat %s > %s" % (input.all_fa, output.merge_fa)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "source activate repeatmasker && RepeatMasker -species human -pa %s -gff -dir %s %s > %s 2>&1" % (threads, params.repeatmasker_dir, output.merge_fa, log)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "source activate base && python %s -re_in %s -re_bed %s" % (params.script, params.repeatmasker_out, output.alpha_hsat23_bed)
        print(cmd3)
        os.system(cmd3)
        cmd4 = "dna-brnn -Ai %s -t%s %s > %s" % (params.dna_brnn_model, threads, output.merge_fa, output.dna_brnn_out)
        print(cmd4)
        os.system(cmd4)
        cmd5 = "cat %s %s > %s" % (output.dna_brnn_out, output.alpha_hsat23_bed, output.combine_bed)
        print(cmd5)
        os.system(cmd5)
        cmd6 = "bedtools sort -i %s | bedtools merge -i - > %s" % (output.combine_bed, output.sorted_merge_bed)
        print(cmd6)
        os.system(cmd6)



rule remove_and_cut:
    input:
        sorted_merge_bed = OUT_PATH + "/dna-brnn/sorted_merge_bed.bed",
        two_NRS = OUT_PATH + "/two_NRS/{sample}_NRS_two.fa"
    output:
        non_rep_NRS = OUT_PATH + "/dna-brnn/{sample}_NRS_nonrep.fa"
    params:
        script = SRC_PATH + "/cut_seq_dna-brnn.py"
    run:
        cmd1 = "python %s -bed %s -NRS_fasta %s -out_NRS %s" % (params.script, input.sorted_merge_bed, input.two_NRS, output.non_rep_NRS)
        print(cmd1)
        os.system(cmd1)



def cat_func2(OUT_PATH, SAMPLES):
    temp_list = []
    for i in SAMPLES:
        temp_list.append(OUT_PATH + "/dna-brnn/" + i + "_NRS_nonrep.fa")
    return temp_list



# Remove contaminants base on Kaiju, BLAST and DIAMOND
rule taxid_gain:
    input:
        all_fa = unpack(cat_func2(OUT_PATH, SAMPLES))
    output:
        merge_fa = OUT_PATH + "/kaiju/merge_NRS_nonrep.fa",
        kaiju_out = OUT_PATH + "/kaiju/kaiju_taxid_greedy.out",
        blast_out = OUT_PATH + "/blast/blast_taxid.tab",
        diamond_out = OUT_PATH + "/diamond/diamond_taxid.out"
    params:
        kaiju_nodes = KAIJU_NODES,
        kaiju_database = KAIJU_DATABASE,
        outfmt = "'6 qseqid sseqid bitscore evalue length pident mismatch gapopen qstart qend sstart send qlen slen qcovs qcovhsp stitle staxids sscinames sskingdom'",
        blast_db = BLAST_DB,
        diamond_db = DIAMOND_DB,
        diamond_nodes = DIAMOND_NODES,
        diamond_maps = DIAMOND_MAPS
    threads:
        THREADS
    run:
        cmd1 = "cat %s > %s" % (input.all_fa, output.merge_fa)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "kaiju -t %s -f %s -i %s -z %s -o %s -v" % (params.kaiju_nodes, params.kaiju_database, output.merge_fa, threads, output.kaiju_out)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "blastn -query %s -db %s -outfmt %s -max_hsps 1 -max_target_seqs 5 -out %s -num_threads %s" % (output.merge_fa, params.blast_db, params.outfmt, output.blast_out, threads)
        print(cmd3)
        os.system(cmd3)
        cmd4 = "diamond blastx --db %s --query %s -t /dev/shm/ --threads %s --max-target-seqs 5 --out %s --taxonnodes %s --taxonmap %s --outfmt 6 qseqid sseqid evalue qlen length pident mismatch gapopen qstart qend sstart send stitle staxids" % (params.diamond_db, output.merge_fa, threads, output.diamond_out, params.diamond_nodes, params.diamond_maps)
        print(cmd4)
        os.system(cmd4)


rule remove_non_chordate:
    input:
        kaiju_out = OUT_PATH + "/kaiju/kaiju_taxid_greedy.out",
        blast_out = OUT_PATH + "/blast/blast_taxid.tab",
        diamond_out = OUT_PATH + "/diamond/diamond_taxid.out",
        non_rep_NRS = OUT_PATH + "/dna-brnn/{sample}_NRS_nonrep.fa"
    output:
        nc_nr_NRS = OUT_PATH + "/non_contam/{sample}_NRS_nonrep_nc.fa"
    params:
        script = SRC_PATH + "/kaiju_blastn_diamond_recover.py",
        chordate_taxid = CHORDATE_TAXID
    run:
        cmd1 = "python %s -dia %s -bo %s -ko %s -taxid %s -i %s -o %s" % (params.script, input.diamond_out, input.blast_out, input.kaiju_out, params.chordate_taxid, input.non_rep_NRS, output.nc_nr_NRS)
        print(cmd1)
        os.system(cmd1)


# Extract flanked sequences and anchor
rule anchor_flanked:
    input:
        one_gff = OUT_PATH + "/one_NRS/{sample}_NRS_one.gff",
        nc_nr_NRS = OUT_PATH + "/non_contam/{sample}_NRS_nonrep_nc.fa",
        polished_genome = "/home/wuzhikun/Project/Population/Assembly/Fasta/{sample}_assembly_polish.fasta.gz"
    output:
        left_flank = OUT_PATH + "/anchor/{sample}_flank_left.fa",
        right_flank = OUT_PATH + "/anchor/{sample}_flank_right.fa",
        left_paf = OUT_PATH + "/anchor/{sample}_flank_left.paf",
        right_paf = OUT_PATH + "/anchor/{sample}_flank_right.paf",
        anchor_bed = OUT_PATH + "/anchor/{sample}_flank_anchor.bed",
        filter_bed = OUT_PATH + "/non_contam/{sample}_NRS_nonrep_nc.bed"
    params:
        script = SRC_PATH + "/flank_extract_ctg_per.py",
        ref_mmi = MINIMAP2_MMI,
        script2 = SRC_PATH + "/anchor_from_minimap2.py",
        script3 = SRC_PATH + "/anchor_nc_frow_raw_bed.py"
    threads:
        THREADS
    run:
        cmd1 = "python %s -gff %s -i %s -ctg %s -ol %s -or %s" % (params.script, input.one_gff, input.nc_nr_NRS, input.polished_genome, output.left_flank, output.right_flank)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "minimap2 -DP -t %s %s %s > %s" % (threads, params.ref_mmi, output.left_flank, output.left_paf)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "minimap2 -DP -t %s %s %s > %s" % (threads, params.ref_mmi, output.right_flank, output.right_paf)
        print(cmd3)
        os.system(cmd3)
        cmd4 = "python %s -left %s -right %s -bed %s" % (params.script2, output.left_paf, output.right_paf, output.anchor_bed)
        print(cmd4)
        os.system(cmd4)
        cmd5 = "python %s -bed %s -out %s -fasta %s" % (params.script3, output.anchor_bed, output.filter_bed, input.nc_nr_NRS)
        print(cmd5)
        os.system(cmd5)



def cat_func3(OUT_PATH, SAMPLES):
    temp_list = []
    for i in SAMPLES:
        temp_list.append(OUT_PATH + "/non_contam/" + i + "_NRS_nonrep_nc.fa")
    return temp_list


def cat_func4(OUT_PATH, SAMPLES):
    temp_list = []
    for i in SAMPLES:
        temp_list.append(OUT_PATH + "/non_contam/" + i + "_NRS_nonrep_nc.bed")
    return temp_list


# Clustering of placed sequences
rule de_redundant_placed:
    input:
        all_fa = unpack(cat_func3(OUT_PATH, SAMPLES)),
        all_bed = unpack(cat_func4(OUT_PATH, SAMPLES))
    output:
        merge_fa = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.fa",
        merge_bed = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.bed",
        can_anchor = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.can_anchor.fa",
        cant_anchor = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.fa",
        two_end_clstr = OUT_PATH + "/cdhit/two_end_anchor.clstr",
        one_end_clstr = OUT_PATH + "/cdhit/one_end_anchor.clstr",
        two_end_fa = OUT_PATH + "/cdhit/two_end_anchor.anhit.fa",
        one_end_fa = OUT_PATH + "/cdhit/one_end_anchor.anhit.fa"
    params:
        script = SRC_PATH + "/can_anchor.py",
        script2 = SRC_PATH + "/anchor_cluster.py"
    run:
        cmd1 = "cat %s > %s" % (input.all_fa, output.merge_fa)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "cat %s > %s" % (input.all_bed, output.merge_bed)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "python %s -merge_fa %s -merge_bed %s -can_anchor %s -cant_anchor %s" % (params.script, output.merge_fa, output.merge_bed, output.can_anchor, output.cant_anchor)
        print(cmd3)
        os.system(cmd3)
        cmd4 = "python %s -b %s -f %s -two_end_clstr %s -one_end_clstr %s -two_fa %s -one_fa %s" % (params.script2, params.merge_bed, output.can_anchor, output.two_end_clstr, output.one_end_clstr, output.two_end_fa, output.one_end_fa)
        print(cmd4)
        os.system(cmd4)



# All-versus-all alignments and remove similar
rule de_redundant_unplaced:
    input:
        cant_anchor = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.fa",
        merge_bed = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.bed",
        one_end_fa = OUT_PATH + "/cdhit/one_end_anchor.anhit.fa",
        two_end_fa = OUT_PATH + "/cdhit/two_end_anchor.anhit.fa"
    output:
        cant_anchor_paf = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.paf",
        cant_anchor_paf_part = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.paf.part",
        minihit_fa = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.minihit.fa",
        minihit_clstr = OUT_PATH + "/cdhit/merge_NRS_nonrep_nc.cant_anchor.minihit.fa.clstr",
        combine_two_type = OUT_PATH + "/cdhit/combine_2type.fa",
        combine_two_type_paf = OUT_PATH + "/cdhit/combine_2type.paf",
        combine_two_type_paf_part = OUT_PATH + "/cdhit/combine_2type.paf.part",
        combine_two_type_mihihit_fa = OUT_PATH + "/cdhit/combine_2type.minihit.fa",
        combine_two_type_mihihit_clstr = OUT_PATH + "/cdhit/combine_2type.minihit.clstr",
        two_type_2_two_end = OUT_PATH + "/cdhit/2type_2_twoend.paf",
        maphit_fa = OUT_PATH + "/cdhit/2type_2_twoend.maphit.fa",
        final_NRS = OUT_PATH + "/cdhit/final_3type.NRS.fa",
        final_two = OUT_PATH + "/cdhit/final_3type.NRS.two_end.fa",
        final_one = OUT_PATH + "/cdhit/final_3type.NRS.one_end.fa",
        final_unplaced = OUT_PATH + "/cdhit/final_3type.NRS.unplaced.fa"
    params:
        script = SRC_PATH + "/redun_remove_minimap2_like_cdhit_bedtools.py",
        mincov = MINCOV,
        mindiv = MINDIV,
        script2 = SRC_PATH + "/stat_realign_to_ref.py",
        script3 = SRC_PATH + "/classify_non_redun.py"
    threads:
        THREADS
    run:
        cmd1 = "minimap2 -DP -t %s %s %s > %s" % (threads, input.cant_anchor, input.cant_anchor, output.cant_anchor_paf)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "python %s -paf %s -paf_part %s -i %s -o %s -cl %s -mincov %s -mindiv %s" % (params.script, output.cant_anchor_paf, output.cant_anchor_paf_part, input.cant_anchor, output.minihit_fa, output.minihit_clstr, params.mincov, params.mindiv)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "cat %s %s > %s" % (output.minihit_fa, input.one_end_fa, output.combine_two_type)
        print(cmd3)
        os.system(cmd3)
        cmd4 = "minimap2 -DP -t %s %s %s > %s" % (threads, output.combine_two_type, output.combine_two_type, output.combine_two_type_paf)
        print(cmd4)
        os.system(cmd4)
        cmd5 = "python %s -paf %s -paf_part %s -i %s -o %s -cl %s -mincov %s -mindiv %s" % (params.script, params.combine_two_type_paf, output.combine_two_type_paf_part, output.combine_two_type, output.combine_two_type_mihihit_fa, output.combine_two_type_mihihit_clstr, params.mincov, params.mindiv)
        print(cmd5)
        os.system(cmd5)
        cmd6 = "minimap2 -DP -t %s %s %s > %s" % (threads, input.two_end_fa, output.combine_two_type_mihihit_fa, output.two_type_2_two_end)
        print(cmd6)
        os.system(cmd6)
        cmd7 = "python %s -NRS_one %s -in_paf %s -out_fasta %s" % (params.script2, output.combine_two_type_mihihit_fa, output.two_type_2_two_end, output.maphit_fa)
        print(cmd7)
        os.system(cmd7)
        cmd8 = "cat %s %s > %s" % (output.maphit_fa, input.two_end_fa, output.final_NRS)
        print(cmd8)
        os.system(cmd8)
        cmd9 = "python %s -all_nrs %s -all_bed %s -two %s -one %s -unplaced %s" % (params.script3, output.final_NRS, input.merge_bed, output.final_two, output.final_one, output.final_unplaced)
        print(cmd9)
        os.system(cmd9)



# Calling presence/absence per sample
rule pav_analyze:
    input:
        final_NRS = OUT_PATH + "/cdhit/final_3type.NRS.fa",
        nc_nr_NRS = OUT_PATH + "/non_contam/{sample}_NRS_nonrep_nc.fa",
        final_two = OUT_PATH + "/cdhit/final_3type.NRS.two_end.fa",
        two_end_clstr = OUT_PATH + "/cdhit/two_end_anchor.clstr",
        one_end_clstr = OUT_PATH + "/cdhit/one_end_anchor.clstr"
    output:
        pav_paf = OUT_PATH + "/pav/{sample}_pav.paf",
        pav_file = OUT_PATH + "/pav/{sample}_pav.pav"
    params:
        script = SRC_PATH + "/pav_minimap2_stat.py",
        mincov = MINCOV,
        mindiv = MINDIV
    threads:
        THREADS
    run:
        cmd1 = "minimap2 -DP -t %s %s %s > %s" % (threads, input.final_NRS, input.nc_nr_NRS, output.pav_paf)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "python %s -paf %s -pan %s -two_end %s -mincov %s -mindiv %s -pav %s -one_clstr %s -two_clstr %s" % (params.script, output.pav_paf, input.final_NRS, input.final_two, params.mincov, params.mindiv, output.pav_file, input.two_end_clstr, input.one_end_clstr)
        print(cmd2)
        os.system(cmd2)



# Comparisons to other genomes
rule map_genome:
    input:
        final_NRS = OUT_PATH + "/cdhit/final_3type.NRS.fa"
    output:
        HX1_paf = OUT_PATH + "/other_human_genome/HX1.paf",
        HX1_paf_part = OUT_PATH + "/other_human_genome/HX1.paf.part",
        HX1_pav = OUT_PATH + "/other_human_genome/HX1_pav.pav"
    params:
        HX1_ref = HX1_REF,
        script = SRC_PATH + "/intersection_other_genome_paf.py",
    log:
        HX1_log = OUT_PATH + "/other_human_genome/HX1_pav.log"
    threads:
        THREADS
    run:
        cmd1 = "minimap2 -DP -t %s %s %s > %s" % (threads, params.HX1_ref, input.final_NRS, output.HX1_paf)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "python %s -cpg %s -in_paf %s -pav_out %s -paf_part %s > %s" % (params.script, input.final_NRS, output.HX1_paf, output.HX1_pav, output.HX1_paf_part, log.HX1_log)
        print(cmd2)
        os.system(cmd2)




# Comparisons to the other pan genomes
rule map_pan:
    input:
        final_NRS = OUT_PATH + "/cdhit/final_3type.NRS.fa"
    output:
        hupan_cpg = OUT_PATH + "/other_pan_genome/275Chinese_cpg.paf",
        cpg_hupan = OUT_PATH + "/other_pan_genome/cpg_275Chinese.paf",
        hupan_cpg_part = OUT_PATH + "/other_pan_genome/275Chinese_cpg.paf.part",
        cpg_hupan_part = OUT_PATH + "/other_pan_genome/cpg_275Chinese.paf.part",
        hupan_pav = OUT_PATH + "/other_pan_genome/275Chinese_pav.pav"
    params:
        HUPAN = HUPAN,
        script = SRC_PATH + "/intersection_pan_reciprocal_paf.py",
    log:
        hupan_log = OUT_PATH + "/other_pan_genome/275Chinese_pav.log"
    threads:
        THREADS
    run:
        cmd1 = "minimap2 -DP -t %s %s %s > %s" % (threads, input.final_NRS, params.HUPAN, output.hupan_cpg)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "minimap2 -DP -t %s %s %s > %s" % (threads, params.HUPAN, input.final_NRS, output.cpg_hupan)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "python %s -cpg %s -cpg_ref %s -other_ref %s -pav_out %s -paf_part_cpg %s -paf_part_other %s > %s" % (params.script, input.final_NRS, output.hupan_cpg, output.cpg_hupan, output.hupan_pav, output.hupan_cpg_part, output.cpg_hupan_part, log.hupan_log)
        print(cmd3)
        os.system(cmd3)



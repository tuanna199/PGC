###################################### sub-group  genotype ################################
rule SubGenotype:
    input:
        vcf = IN_PATH + "/population/Merge/Sample_SV_common_CN329.vcf",
    output:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        freq = IN_PATH + "/population/subgroup/Sample_SV_genotype_numeric.txt",
    params:
        SVGenotype = SRC_DIR + "/SVGenotype.py",
    log:
        IN_PATH + "/log/SubGenotype.log", 
    run:
        shell("python {params.SVGenotype} --vcf {input.vcf} --out {output.geno} --frequency {output.freq} >{log} 2>&1")


rule SubFrequencyFilt:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
    output:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype_filt.txt",
        freq = IN_PATH + "/population/subgroup/Sample_SV_genotype_frequency.txt",
    params:
        SVGenoFreqFilt = SRC_DIR + "/SVGenoFreqFilt.py",
        freqThreshold = 0, #0.05
    log:
        IN_PATH + "/log/SubFrequencyFilt.log", 
    run:
        shell("python {params.SVGenoFreqFilt} --genotype {input.geno} --freqThreshold {params.freqThreshold} --out {output.geno} --frequency {output.freq} > {log} 2>&1")


rule SubMeta2PED:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype_filt.txt",
    output:
        Ped = IN_PATH + "/population/subgroup/Sample_SV_geno.ped",
        Map = IN_PATH + "/population/subgroup/Sample_SV_geno.map",
    params:
        ### need chromosome X
        # Triometa2ped = SRC_DIR + "/Triometa2ped.py",
        meta2ped = SRC_DIR + "/meta2ped.py",
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/SubMeta2PED.log", 
    run:
        ### just get autosomal SV
        shell("python {params.meta2ped} --meta {params.metafile} --genotype {input.geno} --map {output.Map} --ped {output.Ped} >{log} 2>&1")




rule Subped2bed:
    input:
        Ped = IN_PATH + "/population/subgroup/Sample_SV_geno.ped",
        Map = IN_PATH + "/population/subgroup/Sample_SV_geno.map",
    output:
        bed = IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
        bim = IN_PATH + "/population/subgroup/Sample_SV_geno.bim",
        fam = IN_PATH + "/population/subgroup/Sample_SV_geno.fam",
    log:
        IN_PATH + "/log/Subped2bed.log", 
    run:
        ### plink run in mu01
        pedPrefix = input.Ped.rstrip(".ped")
        shell("plink --file {pedPrefix}  --make-bed --out  {pedPrefix} >{log} 2>&1")



############################################################################################




######################################### Admixture #########################################
rule admixture:
    input:
        bed = IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
        bim = IN_PATH + "/population/subgroup/Sample_SV_geno.bim",
        fam = IN_PATH + "/population/subgroup/Sample_SV_geno.fam",
    output:
        log = IN_PATH + "/population/admixture/log1.out",
        cv = IN_PATH + "/population/admixture/admixture_cv.txt",
        pop = IN_PATH + "/population/admixture/Sample_SV_geno.2.Q",
        pop3 = IN_PATH + "/population/admixture/Sample_SV_geno.3.Q",
    threads:
        THREADS * ThreadFold
    params:
        outdir = IN_PATH + "/population/admixture",
        admixture = config["admixture"],
    log:
        IN_PATH + "/log/admixture.log", 
    run:
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        ### 1 2 3 4 5 6 7 8 9 10
        cmd = "cd %s && for K in 1 2 3 4 5 6 7 8 9 10 ; do %s  --cv  %s  $K  -j%d | tee log${K}.out >> %s; done " % (params.outdir, params.admixture, input.bed, threads, log)
        os.system(cmd)
        cmd1 = "grep -h CV %s/log*.out > %s" % (params.outdir, output.cv)
        os.system(cmd1)


rule AssignPopulation:
    input:
        cv = IN_PATH + "/population/admixture/admixture_cv.txt",
        pop = IN_PATH + "/population/admixture/Sample_SV_geno.2.Q",
        ped = IN_PATH + "/population/subgroup/Sample_SV_geno.ped",
    output:
        pop = IN_PATH + "/population/admixture/Sample_assign_population.xls",
    params:
        admixturePop = SRC_DIR + "/admixturePop.py",
        number = 2, ### the given population number
    log:
        IN_PATH + "/log/AssignPopulation.log",
    run:
        shell("python {params.admixturePop} --cv {input.cv} --ped {input.ped} --population {input.pop} --out {output.pop} --number {params.number} >{log} 2>&1")



rule StructurePlot:
    input:
        pop = IN_PATH + "/population/admixture/Sample_assign_population.xls",
    output:
        pdf = IN_PATH + "/population/admixture/Sample_assign_population.pdf",
    params:
        structurePlot = SCRIPT_DIR + "/structurePlot.R",
        width = 8,
        height = 4,
    log:
        IN_PATH + "/log/StructurePlot.log",
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.structurePlot, input.pop, output.pdf, params.width, params.height, log)
        os.system(cmd)

##############################################################################################



########################################### eigenstrat ##############################################
rule geno2eigenstrat:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype_filt.txt",
        # geno = IN_PATH + "/population/genotype/Sample_SV_genotype.txt",
    output:
        ind = IN_PATH + "/population/eigenstrat/Sample_SV.ind",
        snp = IN_PATH + "/population/eigenstrat/Sample_SV.snp",
        geno = IN_PATH + "/population/eigenstrat/Sample_SV.geno",
        snp_temp = IN_PATH + "/population/eigenstrat/Sample_SV_temp.snp",
        snp_number = IN_PATH + "/population/eigenstrat/Sample_SV_number.snp",
        snp_convert = IN_PATH + "/population/eigenstrat/Sample_SV_convert.snp",
    params:
        meta2eigenstrat = SRC_DIR + "/meta2eigenstrat.py",
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/geno2eigenstrat.log", 
    run:
        shell("python {params.meta2eigenstrat} --meta {params.metafile} --genotype {input.geno} --ind {output.ind} --geno {output.geno} --snp {output.snp} >{log} 2>&1")
        cmd = """ awk '{print NR"\t"$s}' %s > %s """ % (output.snp, output.snp_temp)
        os.system(cmd)
        shell("cut -f 1-2 {output.snp_temp} > {output.snp_number}")
        shell("cut -f 1,3- {output.snp_temp} > {output.snp_convert}")


rule smartpca:
    input:
        ind = IN_PATH + "/population/eigenstrat/Sample_SV.ind",
        snp = IN_PATH + "/population/eigenstrat/Sample_SV_convert.snp",
        geno = IN_PATH + "/population/eigenstrat/Sample_SV.geno",
    output:
        pca = IN_PATH + "/population/eigenstrat/Sample_SV.pca.txt",
        eva = IN_PATH + "/population/eigenstrat/Sample_SV.eval.txt",
    params:
        smartpca = config["smartpca"],
        plot = IN_PATH + "/population/eigenstrat/Sample_SV.plot.txt",
    log:
        IN_PATH + "/log/smartpca.log", 
    run:
        shell('echo "perl {params.smartpca} -i {input.geno} -a {input.snp} -b {input.ind} -k 10 -o {output.pca} -p {params.plot} -e {output.eva} -l {log}" ')
        ### export PATH=/home/wuzhikun/software/Eigensoft/EIG-7.2.1/bin:$PATH
        ### mu01
        shell("export PATH=/home/wuzhikun/software/Eigensoft/EIG-7.2.1/bin:$PATH && perl {params.smartpca} -i {input.geno} -a {input.snp} -b {input.ind} -k 10 -o {output.pca} -p {params.plot} -e {output.eva} -l {log}")


rule IndPca:
    input:
        ind = IN_PATH + "/population/eigenstrat/Sample_SV.ind",
        pca = IN_PATH + "/population/eigenstrat/Sample_SV.pca.txt",
    output:
        pca = IN_PATH + "/population/eigenstrat/Sample_PCA.txt",
    params:
        eigenstrat2PCA = SRC_DIR + "/eigenstrat2PCA.py",
    log:
        IN_PATH + "/log/IndPca.log", 
    run:
        shell("python {params.eigenstrat2PCA} --ind {input.ind} --pca {input.pca} --out {output.pca} >{log} 2>&1")


rule PCA1:
    input:
        pca = IN_PATH + "/population/eigenstrat/Sample_PCA.txt",
        # structure = IN_PATH + "/population/admixture/Sample_assign_population.xls",
    output:
        pca0 = IN_PATH + "/population/eigenstrat/Sample_PCA-0.txt",
        # pca1 = IN_PATH + "/population/eigenstrat/Sample_PCA-1.txt",
        pca_province = IN_PATH + "/population/eigenstrat/Sample_PCA-province.txt",
        pca_province23 = IN_PATH + "/population/eigenstrat/Sample_PCA23-province.txt",
    params:
        PCAaddGroup = SRC_DIR + "/PCAaddGroup.py",
        metafile = config["metafile"],
    log:
        IN_PATH + "/log/PCA1.log", 
    run:
        cmd = """awk '{print $1"\t"$1"\t"$2"\t"$3}' %s | sed '1d' | sed '1i FID IID PCA1 PCA2' | sed 's/ /\t/g' > %s""" % (input.pca, output.pca0)
        os.system(cmd)
        # shell("python {params.PCAaddGroup} --structure {input.structure} --pca {output.pca0} --out {output.pca1} --column group > {log} 2>&1")
        shell("python {params.PCAaddGroup} --structure {params.metafile} --pca {output.pca0} --out {output.pca_province} --column province >> {log} 2>&1")
        cmd = """awk '{print $1"\t"$1"\t"$3"\t"$4}' %s | sed '1d' | sed '1i FID IID PCA2 PCA3' | sed 's/ /\t/g' > %s""" % (input.pca, output.pca0)
        os.system(cmd)
        shell("python {params.PCAaddGroup} --structure {params.metafile} --pca {output.pca0} --out {output.pca_province23} --column province >> {log} 2>&1")


rule PCARegion:
    input:
        pca = IN_PATH + "/population/eigenstrat/Sample_PCA-province.txt",
        pca23 = IN_PATH + "/population/eigenstrat/Sample_PCA23-province.txt",
    output:
        pca = IN_PATH + "/population/eigenstrat/Sample_PCA-province_region.txt",
        pca23 = IN_PATH + "/population/eigenstrat/Sample_PCA23-province_region.txt",
    params:
        metafile = config["metafile"],
    run:
        PCA_sample_region(params.metafile, input.pca, output.pca)
        PCA_sample_region(params.metafile, input.pca23, output.pca23)



rule PCAPlot:
    input:
        # pca = IN_PATH + "/population/eigenstrat/Sample_PCA-1.txt",
        pca_province = IN_PATH + "/population/eigenstrat/Sample_PCA-province_region.txt",
        pca23 = IN_PATH + "/population/eigenstrat/Sample_PCA23-province_region.txt",
    output:
        # pdf = IN_PATH + "/population/eigenstrat/Sample_PCA12.pdf",
        pdf2 = IN_PATH + "/population/eigenstrat/Sample_PCA12-province.pdf",
        pca23 = IN_PATH + "/population/eigenstrat/Sample_PCA23-province.pdf",
    params:
        genoPCA = SCRIPT_DIR + "/genoPCA.R",
        width = 5,
        height = 4,
        width2 = 5,
    log:
        IN_PATH + "/log/PCAPlot.log", 
    run:
        # cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.genoPCA, input.pca, output.pdf, params.width, params.height, log)
        # print(cmd)
        # os.system(cmd)
        cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.genoPCA, input.pca_province, output.pdf2, params.width2, params.height, log)
        print(cmd2)
        os.system(cmd2)
        cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.genoPCA, input.pca23, output.pca23, params.width2, params.height, log)
        print(cmd2)
        os.system(cmd2)





################################################################################################



#################################### shuffle Fst #############################
rule shuffleGroup:
    input:
        group = "/home/wuzhikun/Project/Population/geno_group.txt",
    output:
        group = IN_PATH + "/population/ShuffleFst/group/Sample_SV_genotype_fst_{number}.txt",
    params:
        shuffleGroup = SRC_DIR + "/shuffleGroup.py",
    threads:
        THREADS
    run:
        shell("python {params.shuffleGroup} --input {input.group} --out {output.group}")
    

rule shuffleFst:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        structure = IN_PATH + "/population/ShuffleFst/group/Sample_SV_genotype_fst_{number}.txt",
    output:
        fst = IN_PATH + "/population/ShuffleFst/Sample_SV_genotype_fst_{number}.txt",
    params:
        PopGenoFst = SRC_DIR + "/PopGenoFst.py",
        target = "South,North",
    threads:
        THREADS
    log:
        IN_PATH + "/log/shuffleFst/{number}.log", 
    run:
        shell("python {params.PopGenoFst} --structure {input.structure} --genotype {input.geno} --target {params.target} --out {output.fst} > {log} 2>&1")



rule FstSummary:
    input:
        fst = expand(IN_PATH + "/population/ShuffleFst/Sample_SV_genotype_fst_{number}.txt", number=SHUFFLENUM),
    output:
        fst = IN_PATH + "/population/Fst/Fst_permutation_values.txt",
        fstSort = IN_PATH + "/population/Fst/Fst_permutation_values_sorted.txt",
    run:
        FST = input.fst
        for i in range(len(FST)):
            f = FST[i]
            if i == 0:
                cmd = "sort -k 7nr %s | cut -f 7 | sed -n '1p' > %s" % (f, output.fst)
            else:
                cmd = "sort -k 7nr %s | cut -f 7 | sed -n '1p' >> %s" % (f, output.fst)
            os.system(cmd)
        shell("sort -k 1nr {output.fst} > {output.fstSort}")


##############################################################################





###################################### sub-population Fst ########################
rule PopFst:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        # structure =  IN_PATH + "/population/admixture/Sample_assign_population.xls",
        structure = "/home/wuzhikun/Project/Population/geno_group.txt",
    output:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst.txt",
    params:
        PopGenoFst = SRC_DIR + "/PopGenoFst.py",
        # target = "P1,P2",
        target = "South,North",
    log:
        IN_PATH + "/log/PopFst.log", 
    run:
        shell("python {params.PopGenoFst} --structure {input.structure} --genotype {input.geno} --target {params.target} --out {output.fst} > {log} 2>&1")


rule FstHist:
    input:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst.txt",
    output:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_temp.txt",
        pdf = IN_PATH + "/population/Fst/Sample_SV_genotype_fst.pdf",
    params:
        FstHist = SCRIPT_DIR + "/FstHist.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/FstHist.log", 
    run:
        cmd = """awk '{if($7!="0.000"){print $0}}' %s | grep -v 'nan'  > %s """ % (input.fst, output.fst)
        os.system(cmd)
        shell("Rscript {params.FstHist} --input {output.fst} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")


rule FstPlot:
    input:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst.txt",
    output:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_table.txt",
        manhattan = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_manhattan.pdf",
    params:
        FstPlot = SCRIPT_DIR + "/FstPlot.R",
    log:
        IN_PATH + "/log/FstPlot.log", 
    run:
        shell("""python -c 'print("CHR\tBP\tP\tSNP")' > {output.fst}""")
        cmd = """awk '{print $1"\t"$2"\t"$7"\t"$1"_"$2"-"$3"_"$4"-"$5"-"$6}' %s | grep -v 'nan' | sed 's/^X/23/g' >> %s """ % (input.fst, output.fst)
        os.system(cmd)
        ### manhattan plot
        shell("Rscript {params.FstPlot} --input {output.fst} --manhattan {output.manhattan} > {log} 2>&1")



rule FstSig:
    input:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst.txt",
    output:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig.txt",
    log:
        IN_PATH + "/log/FstSig.log", 
    run:
        # cmd = """awk '{if($7>0.04){print $1"\t"$2"\t"$4"\t"$1"_"$2"-"$3"_"$4"-"$5"-"$6}}' %s > %s 2>%s""" % (input.fst, output.fst, log)
        ### $7>0.04
        cmd = """awk '{if($7!="nan" && $7>0.052){print $1"\t"$2"\t"$4"\t"$7}}' %s > %s 2>%s""" % (input.fst, output.fst, log)
        print(cmd)
        os.system(cmd)



rule FstPredOverlap:
    input:
        fst = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig.txt",
    output:
        genepred = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig_overlap.bed",
    params:
        genePred = config["genePred"],
    log:
        IN_PATH + "/log/FstPredOverlap.log",
    run:
        shell("bedtools intersect -a {input.fst} -b {params.genePred} -wa -wb > {output.genepred} 2>{log}")



rule FstGene:
    input:
        overlap = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig_overlap.bed",
    output:
        gene = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig_overlap_gene.txt",
    run:
        shell("cut -f 8 {input.overlap} | sort | uniq > {output.gene}")



rule FstEnrich:
    input:
        gene = IN_PATH + "/population/Fst/Sample_SV_genotype_fst_sig_overlap_gene.txt",
    output:
        kegg = IN_PATH + "/population/Fst/Enrichment/KEGG_2019_Human..enrichr.reports.txt", 
    params:
        outdir = IN_PATH + "/population/Fst/Enrichment",
        geneGSEA = SRC_DIR + "/geneGSEA.py",
        EnrichLibrary = config["EnrichLibrary"],
        number = 1,
    log:
        IN_PATH + "/log/SpecialEnrich.log",
    run:
        ### mu01
        ### --number {params.number} 
        shell("python {params.geneGSEA} --annotation {input.gene} --outdir {params.outdir} --library {params.EnrichLibrary} --method website > {log} 2>&1")


##################################################################################################








#################################### age different #####################################

rule AgeFst:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        structure = "/home/wuzhikun/Project/Population/geno_group_age.txt",
    output:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst.txt",
    params:
        PopGenoFst = SRC_DIR + "/PopGenoFst.py",
        target = "Young,Old",
    log:
        IN_PATH + "/log/AgeFst.log", 
    run:
        shell("python {params.PopGenoFst} --structure {input.structure} --genotype {input.geno} --target {params.target} --out {output.fst} > {log} 2>&1")


rule AgeFstPlot:
    input:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst.txt",
    output:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_table.txt",
        manhattan = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst.pdf",
    params:
        FstPlot = SCRIPT_DIR + "/FstPlot.R",
    log:
        IN_PATH + "/log/AgeFstPlot.log", 
    run:
        shell("""python -c 'print("CHR\tBP\tP\tSNP")' > {output.fst}""")
        cmd = """awk '{print $1"\t"$2"\t"$7"\t"$1"_"$2"-"$3"_"$4"-"$5"-"$6}' %s | grep -v 'nan' | sed 's/^X/23/g' >> %s """ % (input.fst, output.fst)
        os.system(cmd)
        ### manhattan plot
        shell("Rscript {params.FstPlot} --input {output.fst} --manhattan {output.manhattan} > {log} 2>&1")




rule AgeFstSig:
    input:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst.txt",
    output:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_sig.txt",
    log:
        IN_PATH + "/log/FstSig.log", 
    run:
        # cmd = """awk '{if($7>0.04){print $1"\t"$2"\t"$4"\t"$1"_"$2"-"$3"_"$4"-"$5"-"$6}}' %s > %s 2>%s""" % (input.fst, output.fst, log)
        ### $7>0.04
        cmd = """awk '{if($7!="nan" && $7>0.052){print $1"\t"$2"\t"$4"\t"$7}}' %s > %s 2>%s""" % (input.fst, output.fst, log)
        print(cmd)
        os.system(cmd)



rule AgeFstPredOverlap:
    input:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_sig.txt",
    output:
        genepred = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_sig_overlap.bed",
    params:
        genePred = config["genePred"],
    log:
        IN_PATH + "/log/FstPredOverlap.log",
    run:
        shell("bedtools intersect -a {input.fst} -b {params.genePred} -wa -wb > {output.genepred} 2>{log}")





rule AgeFst2:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        structure = "/home/wuzhikun/Project/Population/geno_group_age_30_60.txt",
    output:
        fst = IN_PATH + "/population/Fst/Age/Sample_SV_genotype_age_fst_30_60.txt",
    params:
        PopGenoFst = SRC_DIR + "/PopGenoFst.py",
        target = "Young,Old",
    log:
        IN_PATH + "/log/AgeFst.log", 
    run:
        shell("python {params.PopGenoFst} --structure {input.structure} --genotype {input.geno} --target {params.target} --out {output.fst} > {log} 2>&1")

#########################################################################################



#################################### shuffle Age Fst #############################
rule AgeshuffleGroup:
    input:
        group = "/home/wuzhikun/Project/Population/geno_group_age.txt",
    output:
        group = IN_PATH + "/population/ShuffleFstAge/group/Sample_SV_genotype_fst_{number}.txt",
    params:
        shuffleGroup = SRC_DIR + "/shuffleGroup.py",
    threads:
        THREADS
    run:
        shell("python {params.shuffleGroup} --input {input.group} --out {output.group}")
    

rule AgeshuffleFst:
    input:
        geno = IN_PATH + "/population/subgroup/Sample_SV_genotype.txt",
        structure = IN_PATH + "/population/ShuffleFstAge/group/Sample_SV_genotype_fst_{number}.txt",
    output:
        fst = IN_PATH + "/population/ShuffleFstAge/Sample_SV_genotype_fst_{number}.txt",
    params:
        PopGenoFst = SRC_DIR + "/PopGenoFst.py",
        target = "Young,Old",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ShuffleFstAge/{number}.log", 
    run:
        shell("python {params.PopGenoFst} --structure {input.structure} --genotype {input.geno} --target {params.target} --out {output.fst} > {log} 2>&1")



rule AgeFstSummary:
    input:
        fst = expand(IN_PATH + "/population/ShuffleFstAge/Sample_SV_genotype_fst_{number}.txt", number=SHUFFLENUM),
    output:
        fst = IN_PATH + "/population/ShuffleFstAge/Fst_permutation_values.txt",
        fstSort = IN_PATH + "/population/ShuffleFstAge/Fst_permutation_values_sorted.txt",
    run:
        FST = input.fst
        for i in range(len(FST)):
            f = FST[i]
            if i == 0:
                cmd = "sort -k 7nr %s | cut -f 7 | sed -n '1p' > %s" % (f, output.fst)
            else:
                cmd = "sort -k 7nr %s | cut -f 7 | sed -n '1p' >> %s" % (f, output.fst)
            os.system(cmd)
        shell("sort -k 1nr {output.fst} > {output.fstSort}")


##############################################################################



####################################### Gemma GWAS ################################




# rule kinship:
#     input:
#         bed = IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
#         bim = IN_PATH + "/population/subgroup/Sample_SV_geno.bim",
#         fam = IN_PATH + "/population/subgroup/Sample_SV_geno.fam",
#     output:
#         ks = IN_PATH + "/population/GWAS/output/SV_geno_kS.cXX.txt",
#     params:
#         inPrefix = IN_PATH + "/population/subgroup/Sample_SV_geno",
#         outDir = IN_PATH + "/population/GWAS",
#         outPrefix = "SV_geno_kS",
#     log:
#         IN_PATH + "/log/kinship.log", 
#     run:
#         if not os.path.exists(params.outPrefix):
#             os.mkdir(params.outPrefix)
#         # cmd = "source activate WGS && cd %s && gemma -bfile %s -gk 1 -o %s > %s 2>&1" % (params.outDir, params.inPrefix, params.outPrefix, log)
#         # print(cmd)
#         # os.system(cmd)
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -gk 1 -o {params.outPrefix} > {log} 2>&1")



# rule association:
#     input:
#         bed = IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
#         bim = IN_PATH + "/population/subgroup/Sample_SV_geno.bim",
#         fam = IN_PATH + "/population/subgroup/Sample_SV_geno.fam",
#         ks = IN_PATH + "/population/GWAS/output/SV_geno_kS.cXX.txt",
#         pheno = IN_PATH + "/phenotype/CN329_phenotype_value.txt",
#     output:
#         pheno = IN_PATH + "/phenotype/pheno/phenotype_{trait}.txt",
#         asso = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.txt",
#     threads:
#         THREADS
#     params:
#         inPrefix = IN_PATH + "/population/subgroup/Sample_SV_geno",
#         outDir = IN_PATH + "/population/GWAS",
#         # outPrefix = "pheno1",
#     log:
#         IN_PATH + "/log/association_{trait}.log", 
#     run:
#         ### gemma -bfile /home/wuzhikun/Project/Population/population/plink/Sample_SV_geno -k /home/wuzhikun/Project/Population/population/GWAS/output/SV_geno_kS.cXX.txt -lmm 1 -o pheno1 -p /home/wuzhikun/Project/Population/phenotype/CN329_phenotype_value-1.txt
#         shell("cut -f {wildcards.trait} {input.pheno} > {output.pheno}")
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o pheno_{wildcards.trait} -p {output.pheno} > {log} 2>&1")


# rule ManhattanPlot:
#     input:
#         pheno = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.txt",
#     output:
#         pheno = IN_PATH + "/population/GWAS/output/plot/pheno_{trait}.assoc.txt",
#         manhattan = IN_PATH + "/population/GWAS/output/plot/pheno_{trait}.assoc_mahattan.pdf",
#         qq = IN_PATH + "/population/GWAS/output/plot/pheno_{trait}.assoc_qq.pdf",
#     params:
#         GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
#     log:
#         IN_PATH + "/log/gwasPlot_{trait}.log", 
#     run:
#         shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.pheno}")
#         cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.pheno, output.pheno)
#         print(cmd)
#         os.system(cmd)
#         shell("Rscript {params.GWASPlot} --input {output.pheno} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")


# rule SigLoci:
#     input:
#         asso = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.txt",
#     output:
#         sig = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.sig.txt",
#     run:
#         # cmd = "awk '{if ($12< 5e-8){print $0}}'  %s > %s" % (input.asso, output.sig)
#         cmd = "awk '{if ($12< 1.59e-6){print $0}}'  %s > %s" % (input.asso, output.sig)
#         print(cmd)
#         os.system(cmd)



# rule sigRegion:
#     input:
#         sig = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.sig.txt",
#         anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
#         bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
#     output:
#         bed = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.sig_region.bed",
#         anno = IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.sig_anno.txt",
#     params:
#         TagBedAnno = SRC_DIR + "/TagBedAnno.py",
#     log:
#         IN_PATH + "/log/sigRegion_{trait}.log", 
#     run:
#         shell("python {params.TagBedAnno} --bed {input.bed} --annotation {input.anno} --significant {input.sig} --bedOut {output.bed} --annoOut {output.anno} --tagColumn 2 > {log} 2>&1")




# rule MultipleAsso:
#     ### sig: 2,3,4,5,11,12,19,21,24,31,38,39,40,42,45
#     input:
#         asso = expand(IN_PATH + "/population/GWAS/output/pheno_{trait}.assoc.txt", trait=TRAITS),
#     output:
#         asso = IN_PATH + "/population/GWAS/multiple_traits.assoc.txt",
#     params:
#         MergeManhattan = SRC_DIR + "/MergeManhattan.py",
#         target = "2,3,4,5,11,12,19,21,24,31,38,39,40,42,45",
#     log:
#         IN_PATH + "/log/MultipleAsso.log", 
#     run:
#         Files = ",".join(input.asso)
#         shell("python {params.MergeManhattan} --file {Files} --out {output.asso} --target {params.target} > {log} 2>&1")



    


######################################################################################




# ###################################  GWAS test ###################################

# rule associationTest:
#     input:
#         bed = IN_PATH + "/population/subgroup/Sample_SV_geno.bed",
#         bim = IN_PATH + "/population/subgroup/Sample_SV_geno.bim",
#         fam = IN_PATH + "/population/subgroup/Sample_SV_geno.fam",
#         ks = IN_PATH + "/population/GWAS/output/SV_geno_kS.cXX.txt",
#         pheno = IN_PATH + "/phenotype_test/CN329_phenotype_value.txt",
#     output:
#         pheno = IN_PATH + "/phenotype_test/pheno/phenotype_{trait}.txt",
#         asso = IN_PATH + "/population/GWASTest/output/pheno_{trait}.assoc.txt",
#     threads:
#         THREADS
#     params:
#         inPrefix = IN_PATH + "/population/subgroup/Sample_SV_geno",
#         outDir = IN_PATH + "/population/GWASTest",
#         # outPrefix = "pheno1",
#     log:
#         IN_PATH + "/log/associationTest_{trait}.log", 
#     run:
#         ### gemma -bfile /home/wuzhikun/Project/Population/population/plink/Sample_SV_geno -k /home/wuzhikun/Project/Population/population/GWAS/output/SV_geno_kS.cXX.txt -lmm 1 -o pheno1 -p /home/wuzhikun/Project/Population/phenotype/CN329_phenotype_value-1.txt
#         shell("cut -f {wildcards.trait} {input.pheno} > {output.pheno}")
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o pheno_{wildcards.trait} -p {output.pheno} > {log} 2>&1")


# rule ManhattanPlotTest:
#     input:
#         pheno = IN_PATH + "/population/GWASTest/output/pheno_{trait}.assoc.txt",
#     output:
#         pheno = IN_PATH + "/population/GWASTest/output/plot/pheno_{trait}.assoc.txt",
#         manhattan = IN_PATH + "/population/GWASTest/output/plot/pheno_{trait}.assoc_mahattan.pdf",
#         qq = IN_PATH + "/population/GWASTest/output/plot/pheno_{trait}.assoc_qq.pdf",
#     params:
#         GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
#     log:
#         IN_PATH + "/log/ManhattanPlotTest_{trait}.log", 
#     run:
#         shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.pheno}")
#         cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.pheno, output.pheno)
#         print(cmd)
#         os.system(cmd)
#         shell("Rscript {params.GWASPlot} --input {output.pheno} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")


# ##################################################################################




######################################### Plink GWAS ############################################
# rule gwas:
#     input:
#         Ped = IN_PATH + "/population/subgroup/Sample_SV_geno.ped",
#         Map = IN_PATH + "/population/subgroup/Sample_SV_geno.map",
#         pca = IN_PATH + "/population/eigenstrat/Sample_PCA-0.txt",
#     output:
#         gwas = IN_PATH + "/population/GWAS/Sample_gwas.assoc", ### --assoc 
#         # gwas = IN_PATH + "/population/GWAS/Sample_gwas.assoc.logistic", ### ‐‐logistic
#         gwas1 = IN_PATH + "/population/GWAS/Sample_gwas.assoc.logistic.tab",
#     params:
#         prefix = IN_PATH + "/population/subgroup/Sample_SV_geno",
#         outPreifx = IN_PATH + "/population/GWAS/Sample_gwas",
#     log:
#         IN_PATH + "/log/gwas.log", 
#     run:
#         ### ‐‐logistic  --assoc  --adjust
#         ### assoc
#         shell("plink --file {params.prefix}  --out {params.outPreifx} --assoc > {log} 2>&1")
#         ### logistic
#         shell("plink --file {params.prefix}  --out {params.outPreifx} --logistic > {log} 2>&1")
#         ### logistic PCA
#         # shell("plink --file {params.prefix}  --covar {input.pca} --covar-name PCA1,PCA2  --out {params.outPreifx} --logistic > {log} 2>&1")
#         shell("sed 's/ \+/\t/g' {output.gwas} | cut -f 2- | grep -v 'NA' > {output.gwas1}")


# rule gwasPlot:
#     input:
#         gwas = IN_PATH + "/population/GWASPlink/output/ALT.assoc.linear.tab",
#     output:
#         manhattan = IN_PATH + "/population/GWASPlink/output/ALT.assoc.linear.tab_manhattan.pdf",
#         qq = IN_PATH + "/population/GWASPlink/output/ALT.assoc.linear.tab_qq.pdf",
#     params:
#         GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
#     log:
#         IN_PATH + "/log/gwasPlot.log", 
#     run:
#         shell("Rscript {params.GWASPlot} --input {input.gwas} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")


rule gwas:
    input:
        Ped = IN_PATH + "/population/subgroup/Sample_SV_geno.ped",
        Map = IN_PATH + "/population/subgroup/Sample_SV_geno.map",
        covariant = IN_PATH + "/phenotype/pheno_covariantPCA.txt",
        pheno = IN_PATH + "/phenotype/quant_phenotype.txt",
    output:
        gwas = IN_PATH + "/population/GWASPlink/output/{pheno}.assoc.linear", 
        gwas1 = IN_PATH + "/population/GWASPlink/output/{pheno}.assoc.linear.tab",
    params:
        prefix = IN_PATH + "/population/subgroup/Sample_SV_geno",
        # covarName = "AGE,BMI,SEX,PCA1,PCA2",
        covarName = "AGE,SEX,PCA1,PCA2", ### for BMI
        outPreifx = IN_PATH + "/population/GWASPlink/output/{pheno}",
    log:
        IN_PATH + "/log/gwas_{pheno}.log", 
    run:
        shell("plink --file {params.prefix} --pheno {input.pheno} --pheno-name {wildcards.pheno} --out {params.outPreifx}  --linear --covar {input.covariant} --covar-name  {params.covarName}  > {log} 2>&1")
        ### logistic PCA
        # shell("plink --file {params.prefix}  --covar {input.pca} --covar-name PCA1,PCA2  --out {params.outPreifx} --logistic > {log} 2>&1")
        shell("sed 's/ \+/\t/g' {output.gwas} | cut -f 2- | grep -v 'NA' | grep ADD > {output.gwas1}")



rule ManPlot1:
    input:
        pheno = IN_PATH + "/population/GWASPlink/output/{pheno}.assoc.linear.tab",
    output:
        pheno = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.linear.tab.txt",
    run:
        shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.pheno}")
        cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$9}' >> %s"""  % (input.pheno, output.pheno)
        print(cmd)
        os.system(cmd)



rule ManPlot2:
    input:
        pheno = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.linear.tab.txt",
    output:
        manhattan = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.linear_mahattan.pdf",
        qq = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.linear_qq.pdf",
    params:
        GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
    log:
        IN_PATH + "/log/ManPlot_{pheno}.log", 
    run:
        shell("Rscript {params.GWASPlot} --input {input.pheno} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")




rule SigLoci2:
    input:
        asso = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.linear.tab.txt",
    output:
        sig = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig.txt",
    run:
        # cmd = "awk '{if ($12< 5e-8){print $0}}'  %s > %s" % (input.asso, output.sig)
        cmd = "awk '{if ($4< 1.59e-6){print $0}}'  %s > %s" % (input.asso, output.sig)
        print(cmd)
        os.system(cmd)






rule sigRegion2:
    input:
        sig = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig.txt",
        anno = IN_PATH + "/population/Annotation/Sample_common_SV_genepred_overlap.bed",
        bed = IN_PATH + "/population/bed/Sample_common_SV_DEL_INS_INV_DUP.bed",
    output:
        bed = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig_region.bed",
        anno = IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig_anno.txt",
    params:
        TagBedAnno = SRC_DIR + "/TagBedAnno.py",
    log:
        IN_PATH + "/log/sigRegion2_{pheno}.log", 
    run:
        shell("python {params.TagBedAnno} --bed {input.bed} --annotation {input.anno} --significant {input.sig} --bedOut {output.bed} --annoOut {output.anno} --tagColumn 1 > {log} 2>&1")


def Merge_significant_loci(fileStr, out_file):
    out_h = open(out_file, "w")
    files = fileStr.split(",")
    for f in files:
        f = f.strip()
        fbase = f.split("/")[-1].split(".")[0]
        in_h = open(f, "r")
        for line in in_h:
            line = line.strip()
            if line != "":
                out_h.write("%s\t%s\n" % (fbase, line))
        in_h.close()
    out_h.close()


rule MergeSigLoci:
    input:
        sig = expand(IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig.txt", pheno=PHENOS),
        bed = expand(IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig_region.bed", pheno=PHENOS),
        anno = expand(IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.sig_anno.txt", pheno=PHENOS),
    output:
        sig = IN_PATH + "/population/GWASPlink/output/plot/All_phenotypes.assoc.sig.txt",
        bed = IN_PATH + "/population/GWASPlink/output/plot/All_phenotypes.assoc.sig_region.bed",
        anno = IN_PATH + "/population/GWASPlink/output/plot/All_phenotypes.assoc.sig_anno.txt",
    run:
        SIG = ",".join(sorted(input.sig))
        Merge_significant_loci(SIG, output.sig)
        BED = ",".join(sorted(input.bed))
        Merge_significant_loci(BED, output.bed)
        ANNO = ",".join(sorted(input.anno))
        Merge_significant_loci(ANNO, output.anno)


rule MultipleAsso:
    ### sig: 2,3,4,5,11,12,19,21,24,31,38,39,40,42,45
    input:
        asso = expand(IN_PATH + "/population/GWASPlink/output/plot/{pheno}.assoc.linear.tab.txt", pheno=PHENOS),
    output:
        asso = IN_PATH + "/population/GWASPlink/output/plot/multiple_significant_traits.assoc.txt",
    params:
        MergeManhattan = SRC_DIR + "/MergeManhattan.py",
        target = "BACT,AST,XTAL,RDW-CV,MUCS,ASAL,ALT,CAST,MCH,MCHC,SRC,GLU,LDL-C,MOABS,TG,GCT,NEABS,WBC,HGB,AFP",
    log:
        IN_PATH + "/log/MultipleAsso.log", 
    run:
        Files = ",".join(input.asso)
        shell("python {params.MergeManhattan} --file {Files} --out {output.asso} --target {params.target} > {log} 2>&1")



#############################################################################################


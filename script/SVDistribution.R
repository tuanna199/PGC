#!/usr/bin/env Rscript
# library(CMplot)
### modified CMplot.r
source("/home/wuzhikun/github/NanoHub/script/CMplot.r") 
library(argparser)


### package
### https://cloud.r-project.org/src/contrib/CMplot_3.3.3.tar.gz

#usage: Rscript ~/github/TrioWGS/script/SVDistribution.R --input /home/wuzhikun/Project/NanoTrio/SVStats/Sniffles/minimap2/M628-2_position.xls --pdf /home/wuzhikun/Project/NanoTrio/SVStats/Sniffles/minimap2/M628-2_position.pdf


arg <- arg_parser('Density plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--maxNum', help='The max number for heatmap.')
argv <- parse_args(arg)



SV_distribution <- function(in_file, out_file, maxNum){
    ### in_file
    # 1       76683   76789   DEL
    # 1       83940   84004   DEL
    # 1       136275  136340  INS
    # 1       136905  137040  INS


    SV_pos <- read.table(in_file, sep="\t", header=FALSE)
    colnames(SV_pos) <- c("Chromosome", "Position", "Position2", "Type")
    SV_pos$SNP <- paste(SV_pos$Chromosome, SV_pos$Position, SV_pos$Position2, sep="_")

    SV_position <- data.frame(SNP=SV_pos$SNP, Chromosome=SV_pos$Chromosome, Position=SV_pos$Position, Position2=SV_pos$Position2)

    # chrom <- "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
    # chromAll <- gsub(",", "|", chrom)
    # chroms <- strsplit(chrom, split=",")
    chroms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "X")

    New_position <- SV_position[SV_position$Chromosome %in% chroms,]

    maxNum <- as.numeric(maxNum)
    CMplot(New_position, outName=out_file, maxNum=maxNum, plot.type="d",bin.size=1e6,col=c("darkgreen", "gold2", "red"),file="pdf",memo="",dpi=300, file.output=TRUE, verbose=TRUE, chr.labels=NULL)

}


SV_distribution(argv$input, argv$pdf, argv$maxNum)

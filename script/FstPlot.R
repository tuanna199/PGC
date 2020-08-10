#!/usr/bin/env Rscript
source("/home/wuzhikun/github/NanoHub/script/manhattanModify.R")
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/NanoHub/script/FstPlot.R --input /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype_fst_table.txt --manhattan /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype_fst_table.pdf

arg <- arg_parser('Manhattan for Fst result.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--manhattan', help='Output manhattan plot with pdf format.')
argv <- parse_args(arg)


Manhattan_QQ_plot <- function(infile, manhattan_pdf){
    ### infile:
    #   SNP CHR BP      P
    # rs1   1  1 0.9148
    # rs2   1  2 0.9371
    # rs3   1  3 0.2861
    
    data <- read_tsv(infile)

    pdf(manhattan_pdf, width=10, height=5) #width=12
    ### manhattan(data, cex = 0.5, cex.axis = 0.8, highlight = snpsOfInterest, col = c("blue4", "orange3"), ymax = 12, genomewideline = T)
    ### manhattan(subset(gwasResults, CHR == 1))
    manhattan(data, main = "Fst")
    dev.off()


}


Manhattan_QQ_plot(argv$input, argv$manhattan)

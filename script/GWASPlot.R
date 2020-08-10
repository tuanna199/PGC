#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)
library(qqman)


arg <- arg_parser('Manhattan and QQ plot for GWAS result.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--manhattan', help='Output manhattan plot with pdf format.')
arg <- add_argument(arg, '--qq', help='Output QQ plot with pdf format.')
argv <- parse_args(arg)


Manhattan_QQ_plot <- function(infile, manhattan_pdf, qq_pdf){
    ### infile:
    #   SNP CHR BP      P
    # rs1   1  1 0.9148
    # rs2   1  2 0.9371
    # rs3   1  3 0.2861
    
    data <- read_tsv(infile)

    pdf(manhattan_pdf, width=12, height=5)
    ### manhattan(data, cex = 0.5, cex.axis = 0.8, highlight = snpsOfInterest, col = c("blue4", "orange3"), ymax = 12, genomewideline = T)
    ### manhattan(subset(gwasResults, CHR == 1))
    manhattan(data, genomewideline=-log10(1.59e-6),  suggestiveline=FALSE, main = "Manhattan Plot")
    dev.off()

    pdf(qq_pdf, width=3.5, height=4)
    qq(data$P, main = "Q-Q plot of GWAS p-values")
    dev.off()

}


Manhattan_QQ_plot(argv$input, argv$manhattan, argv$qq)


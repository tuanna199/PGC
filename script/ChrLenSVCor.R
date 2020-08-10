#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(readr)
library(ggpubr)


#usage: Rscript ~/github/NanoHub/script/ChrLenSVCor.R  --SV /home/wuzhikun/Project/NanoTrio/population/Stats/Sample_merge_dist_summary.xls  --length  /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38_length_chr_chr.txt --pdf temp.pdf --width 5 --height 4



arg <- arg_parser('Correlation for SV number and chromosome length.')
arg <- add_argument(arg, '--length', help='The file with chromosome length.')
arg <- add_argument(arg, '--SV', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



Chr_length_SV_correlation <- function(chr_length, SV_number, pdf_file, width, height){
    ### SV_number
    # Chr DEL INS DUP INV TRA Total
    # 1   3814    4325    133 130 537 8939
    # 2   4062    3776    104 295 266 8503
    # 3   2810    2605    82  332 241 6070

    ### chr_length
    # Chr Length  LengthWithoutN
    # 1   248956422   230481012
    # 10  133797422   133262962

    width <- as.numeric(width)
    height <- as.numeric(height)

    SVNumber <- read_tsv(SV_number) 
    rownames(SVNumber) <- SVNumber$Chr

    ChrLength <- read_tsv(chr_length)
    rownames(ChrLength) <- ChrLength$Chr

    mergedData <- merge(SVNumber, ChrLength, by="row.names", all.x=FALSE)
    mergedData$Length <- mergedData$Length / 1000000

    # mergedData$Length <- as.integer(mergedData$Length / 1000000)
    # mergedData$Total <- as.numeric(mergedData$Total)

    print(mergedData)
    correlationPlot <- ggscatter(mergedData, x = "Length", y = "Total", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "Chromosome length (Mb)", ylab = "SV number", cex.axis=0.9, cex.lab=1.0)  +  theme(plot.margin = margin(0.5, 1, 1, 1, "cm")) # cor.coef = TRUE,
    
    correlationPlot
    ggsave(pdf_file, width=width, height=height)
}


Chr_length_SV_correlation(argv$length, argv$SV, argv$pdf, argv$width, argv$height)


#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/NanoHub/script/structurePlot.R --input Sample_assign_population.xls --pdf pop_stru.pdf --height 4 --width 8

arg <- arg_parser('Stack bar plot for population structure.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



structurePlot <- function(infile, pdf_file, width, height){
    colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "darkolivegreen3", "mediumorchid3",  "gold3", "deeppink3", "green4", "lightblue3", "royalblue1",  "thistle1", "skyblue", "purple2", "orangered")
    height <- as.numeric(height)
    width <- as.numeric(width)

    
    data <- read_tsv(infile)

    sortedData <- data[with(data, order(-P1, -P2)), ]

    melt_dt <- reshape2::melt(data)

    colnames(melt_dt) <- c("Sample", "Assign", "Group", "Value")


    melt_dt$Sample <- factor(melt_dt$Sample, order=TRUE, levels=as.vector(sortedData$geno))


    barPlot <- ggplot(melt_dt, aes(x=Sample, y=Value, fill=Group)) +
        theme_bw() +
        geom_bar(stat="identity",width=0.75) +
        xlab("Samples") + ylab("K=2") +
        scale_y_continuous(expand = c(0,0)) + 
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_fill_manual(values=colors)

    barPlot
    ggsave(pdf_file, width=width, height=height)
}

structurePlot(argv$input, argv$pdf, argv$width, argv$height)

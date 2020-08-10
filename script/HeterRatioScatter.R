#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/NanoHub/script/HeterRatioScatter.R --input Sample_SV_heter2homo.txt --pdf Sample_SV_heter2homo.pdf --width 6 --height 4


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


RatioBoxPlot <- function(in_file, pdf_file, height, width){
    height <- as.numeric(height)
    width <- as.numeric(width)

    CHRS <- c("1","2","3","4","5","6","7","8","9","10","11", "12","13","14","15","16","17","18","19","20","21","22")
    
    data <- read_tsv(in_file)
    melt_dt <- reshape2::melt(data)
    colnames(melt_dt) <- c("Chr", "Sample", "Value")

    melt_dt <- melt_dt[which(melt_dt$Chr %in% CHRS), ]

    melt_dt$Chr=factor(melt_dt$Chr,levels=c("1","2","3","4","5","6","7","8","9","10","11", "12","13","14","15","16","17","18","19","20","21","22"))

    
    BoxPlot <- ggplot(melt_dt, aes(x=Chr, y=Value)) +
        geom_boxplot() + 
        ylim(0, 15) +
        xlab("Chromosome") + ylab("Heter/Homo ratio") +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +  ##panel.border=element_blank(),
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + #angle = 90, hjust = 1
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme_bw()
        
    BoxPlot
    ggsave(pdf_file, width=width, height=height)
}

RatioBoxPlot(argv$input, argv$pdf, argv$height, argv$width)


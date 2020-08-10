#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/NanoHub/script/OverlapSVLenChrBar.R --input /home/wuzhikun/Project/NanoTrio/population/bed/overlap/SV_overlap_filt_chr_length_summary.xls --pdf1 temp1.pdf  --pdf2 temp2.pdf --width 6 --height 4


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf1', help='output file with pdf format.')
arg <- add_argument(arg, '--pdf2', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)




SV_length_chr_bar_Plot <- function(infile, pdf1, pdf2, width, height){

    height <- as.numeric(height)
    width <- as.numeric(width)



    levels=c("DEL", "INS", "DUP", "INV", "TRA", "INVDUP")
    colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")



    data <- read_tsv(infile)
    
    colors <- c("gold2", "limegreen", "royalblue1", "tomato")
    data$SVType <- factor(data$SVType, order=TRUE, levels=c("INV", "DEL", "INS", "DUP"))

    # colors <- c("tomato", "gold2", "limegreen", "royalblue1")
    # data$SVType <- factor(data$SVType, order=TRUE, levels=c("DUP", "INV", "DEL", "INS"))

    data$Chr=factor(data$Chr,levels=c("1","2","3","4","5","6","7","8","9","10","11", "12"
    ,"13","14","15","16","17","18","19","20","21","22","X"))


    barPlot <- ggplot(data, aes(x=Chr, y=log2(Number), fill=SVType)) +
        theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("Chromosome") + ylab("SV Number (log10)") +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank()) +
        theme(axis.text.x=element_blank()) + # axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())
        # theme(legend.position="none")

    barPlot
    ggsave(pdf1, width=width, height=height)


    barPlot2 <- ggplot(data, aes(x=Chr, y=log10(TotalLength), fill=SVType)) +
        theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("Chromosome") + ylab("Total length (log10(bp))") +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme(panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),plot.background=element_blank()) +
        theme(axis.text.x=element_blank()) + # axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())
        # theme(legend.position="none")

    barPlot2
    ggsave(pdf2, width=width, height=height)
}
    



SV_length_chr_bar_Plot(argv$input, argv$pdf1, argv$pdf2, argv$width, argv$height)
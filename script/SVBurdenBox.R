#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


#usage: Rscript ~/github/NanoHub/script/LengthHeterBox.R --input temp_del.txt --pdf temp_del.pdf --wdith 4 --height 4

arg <- arg_parser('Box plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output box file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)

SVBurden_boxplot <- function(in_file, pdf_file, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)
    melt_dt <- reshape2::melt(data[1:nrow(data)-1, ])
    colnames(melt_dt) <- c("Category", "Sample", "Type", "Value")

    BoxPlot <- ggplot(melt_dt, aes(x=Category, y=Value)) +
        geom_violin(aes(fill=Type), trim=FALSE) +
        geom_boxplot(width=0.2) +
        xlab("Category") + ylab("SV burden") +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +  ##panel.border=element_blank(),
        theme(axis.text.x=element_blank()) + #axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + #angle = 90, hjust = 1
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        #scale_fill_manual(values=colors)  + 
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme_bw()
    BoxPlot <- BoxPlot + facet_wrap( ~ Type, scales="free")

    ggsave(pdf_file, width=width, height=height)

}

SVBurden_boxplot(argv$input, argv$pdf, argv$width, argv$height)

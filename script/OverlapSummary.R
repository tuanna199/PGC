#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(scales)
library(readr)
library(plyr)

#usage: Rscript ~/github/NanoHub/script/OverlapSummary.R --input /home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_overlap_summary.txt --pdf /home/wuzhikun/Project/Population/population/bed/Cell2019/Sample_overlap_summary.pdf --width 5 --height 4


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


SV_type_freq <- function(in_file, pdf_file, height, width){



    # colors <- c( "slategray3", "mediumseagreen",  "dodgerblue3", "gold2", "orchid2",  "tomato2")
    colors <- c(  "mediumseagreen",  "dodgerblue3", "gold2", "tomato2")

    type_data <- read_tsv(in_file)
    type_data$Number <- log10(type_data$Number)
    type_data$Number <- ifelse(type_data$Number == -Inf, 0, type_data$Number)
    print(type_data)

    # type_data$Database <- factor(type_data$Database, levels=c("dbVar", "DGV", "LRS15", "WGS911",  "WGS17795", "gnomAD" ))
    type_data$Database <- factor(type_data$Database, levels=c("LRS15", "DGV", "WGS911", "gnomAD" ))
    type_data$Type <- factor(type_data$Type, levels=c("DEL", "INS", "DUP", "INV"))

    # p <- ggplot(type_data, aes(x=Type, y=Number, fill=Database)) + 
    #     geom_bar(stat="identity", width=0.7, position ="dodge") +

    p <- ggplot(type_data,aes(x=Type, y=Number, fill=Database)) +
        geom_bar(position ="dodge", stat="identity", width=0.8) +
        geom_text(aes(label= 10^ Number), position=position_dodge(width=0.9), hjust = 0, vjust=0, size = 3, angle=60) +
        xlab("Types") + 
        ylab(expression(paste("Number (log"[10], ")"))) +
        # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
        #       labels = trans_format("log10", math_format(10^.x))) +
        theme(panel.border=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        # theme(axis.text.x=element_blank(), axis.ticks=element_blank()) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        scale_fill_manual(values=colors) +
        theme(legend.position = "right") + # bottom # right
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank()) +
        scale_y_continuous(expand=c(0,0)) +
        # ylim(0, max(type_data$Number)+0.2)
        ylim(0, 5) 
    p



    height <- as.numeric(height)
    width <- as.numeric(width)
    ggsave(pdf_file, width=width, height=height)

}
 

SV_type_freq(argv$input, argv$pdf, argv$height, argv$width)
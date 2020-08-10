#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


#usage: Rscript ~/github/NanoHub/script/SVTypeBarMultiple.R --input Samples_SV_summary.xls --hist hist.pdf --box box.pdf --width 8 --height 4

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--hist', help='output hist file with pdf format.')
arg <- add_argument(arg, '--box', help='output box file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)




boxPlot <- function(melt_dt, pdf_file, colors){

    melt_dt$Value <- log10(melt_dt$Value)
    plot_box <- ggplot(melt_dt, aes(x=reorder(Type, -Value, median), y=Value, fill=Type)) +
        geom_boxplot()  + 
        theme_bw() + 
        scale_fill_manual(values=colors) +
        theme(axis.text.x=element_text(angle = 0, hjust = 1, vjust=0.5)) + 
        # ylim(0, 15000) +
        theme(plot.title = element_blank()) + #, axis.title.x = element_blank(), axis.title.y = element_blank()) + 
        guides(fill=FALSE) + 
        xlab("Type") + ylab("Number") +
        # scale_y_continuous(limits=c(0,0.8)) + 
        # theme(plot.margin = margin(0.3, 0.3, 0.5, 0.5, "cm")) +
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())
    plot_box
    ggsave(pdf_file, width=5, height=4)
}

 
 
multipleStackBar <- function(melt_dt, pdf_file, colors, height, width){
    melt_dt$Type <- factor(melt_dt$Type, levels=c("Singleton", "Rare", "Low",  "Common" ))
    
    p <- ggplot(melt_dt, aes(x=Sample, y=Value, fill=Type)) + 
        theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("Samples") + ylab("Number") +
        theme(panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
        theme(axis.ticks = element_line(size=rel(0.5))) + #,axis.ticks=element_blank()
        theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 0.5, size=rel(0.5))) +
        theme(plot.margin = margin(0.3,0.5,0.5,0.3, "cm")) +
        scale_fill_manual(values=colors) +
        theme(legend.position = "top") + # "right", "bottom"
        theme(legend.text=element_text(size=rel(0.8))) +
        theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())
    p

    height <- as.numeric(height)
    width <- as.numeric(width)
    ggsave(pdf_file, width=width, height=height)
}

 
 


### colours
colors <- c( "aquamarine3", "royalblue2", "purple2", "orangered3")


in_file <- argv$input
summary_data <- read_tsv(in_file)
summary_data <- summary_data[, 1:5]
melt_dt <- reshape2::melt(summary_data)
colnames(melt_dt) <- c("Sample", "Type", "Value")



boxPlot(melt_dt, argv$box, colors)
multipleStackBar(melt_dt, argv$hist, colors, argv$height, argv$width)

 
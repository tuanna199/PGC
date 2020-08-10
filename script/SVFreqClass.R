#!/usr/bin/Rscript
library(ggplot2)
library(readr)
library(argparser)
library(scales)
library(plyr)

#usage: Rscript ~/github/NanoHub/script/SVFreqClass.R --input Sample_type_freq.xls --pdf Sample_type_freq.pdf --height 4 --width 4  --pdf2  Sample_type_freq_percentage.pdf


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, '--pdf2', help='output hist file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


SV_type_freq <- function(in_file, pdf_file, pdf_file2, height, width){


    colors <- c( "mediumseagreen", "dodgerblue3", "gold2", "tomato2")
    # colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")

	type_data <- read_tsv(in_file)
	melt_dt <- reshape2::melt(type_data)
	colnames(melt_dt) <- c('Type', 'Content', 'Value')

	# melt_dt$Type <- factor(melt_dt$Type, levels=c( "DEL", "INS", "TRA", "INV",  "DUP"))

    # melt_dt$Type <- factor(melt_dt$Type, levels=c("DEL", "INS", "DUP", "INV", "TRA", "INVDUP"))

    melt_dt$Content <- factor(melt_dt$Content, levels=c("Singleton", "Rare", "Low",  "Common" ))
    melt_dt$Type <- factor(melt_dt$Type, levels=c("DEL", "INS",  "DUP", "INV" ))

	### trans to log10 for values
    melt_dt$Value <- log10(melt_dt$Value)

    # p <- ggplot(melt_dt, aes(x=Type, y=Value, fill=Content)) + 
    #     geom_bar(stat="identity", width=0.7, position ="dodge") +

    p <- ggplot(melt_dt, aes(x=Type, y=Value, fill=Content)) +
        geom_bar(position ="dodge", stat="identity", width=0.8) +
        geom_text(aes(label= 10^ Value), position=position_dodge(width=0.9), hjust = 0, vjust=0, size = 3, angle=60) +
        xlab("Types") + ylab(expression(paste("Number (log"[10], ")"))) +
        # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
        #       labels = trans_format("log10", math_format(10^.x))) +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
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
        ylim(0, 5)
    p



   	height <- as.numeric(height)
   	width <- as.numeric(width)
  	ggsave(pdf_file, width=width, height=height)





    melt_dt2 <- ddply(melt_dt, "Type", transform, percent_Value = Value / sum(Value) * 100)

    p2 <- ggplot(melt_dt2, aes(x=Type, y=percent_Value, fill=Content)) + 
        geom_bar(stat="identity",  width=0.7) +
        # coord_flip() +
        theme_bw() +
        xlab("Types") + ylab("Percentage (%)") +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        # theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 0.5, size=rel(0.5))) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1,  "cm")) +
        scale_fill_manual(values=colors) +
        theme(legend.position = "right") + # bottom
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm"), legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank()) +
        scale_y_continuous(expand=c(0,0))

    p2
    height <- as.numeric(height)
    width <- as.numeric(width)
    ggsave(pdf_file2, width=width, height=height)
}


SV_type_freq(argv$input, argv$pdf, argv$pdf2, argv$height, argv$width)
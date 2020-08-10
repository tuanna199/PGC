#!/usr/bin/Rscript
library(ggplot2)
library(readr)
library(argparser)
library(scales)


#usage: Rscript ~/github/NanoHub/script/SVFreqClass.R --input Sample_type_freq.xls --pdf Sample_type_freq.pdf --height 4 --width 4  --pdf2  Sample_type_freq_percentage.pdf


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)





    
    
disease_category <- function(in_file, pdf_file, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)


    colors <- c("mediumseagreen", "dodgerblue3", "gold2", "tomato2")

    type_data <- read_tsv(in_file)
    type_data$Number <- log2(type_data$Number)
    type_data$Number <- ifelse(type_data$Number == -Inf, 0, type_data$Number)
    print(type_data)

    type_data$Category <- factor(type_data$Category, levels=c("Singleton", "Rare", "Low",  "Common" ))
    type_data$Database <- factor(type_data$Database, levels=c("GWAS", "OMIM", "COSMIC")) #"Clinvar", "Biobank", 

    # p <- ggplot(type_data, aes(x=Database, y=Number, fill=Category)) + 
    #     geom_bar(stat="identity", width=0.7, position ="dodge") +

    p <- ggplot(type_data,aes(x=Database, y=Number, fill=Category)) +
        geom_bar(position ="dodge", stat="identity", width=0.8) +
        geom_text(aes(label= 2^ Number), position=position_dodge(width=0.9), hjust = 0, vjust=0, size = 3, angle=60) +
        xlab("") + ylab(expression(paste("Number (log"[2], ")"))) +
        # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
        #       labels = trans_format("log10", math_format(10^.x))) +
        theme(panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank(), axis.line = element_line(colour = "black")) +
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
        # ylim(0, max(type_data$Number)+0.2) +
        ylim(0, 11)
    p

    ggsave(pdf_file, width=width, height=height)

}

disease_category(argv$input, argv$pdf, argv$width, argv$height)


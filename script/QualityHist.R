#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(readr)

arg <- arg_parser('hist of summary for quality control.')
arg <- add_argument(arg, '--input', help='In file with read quality.')
arg <- add_argument(arg, '--outPrefix', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


#usage: Rscript ~/github/NanoHub/script/QualityHist.R --input Samples_quality_summary.xls --pdf1 temp1.pdf --pdf2 temp2.pdf --width 6 --height 4


Quality_hist <- function(infile, outPrefix, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)
    data$Clean_read_length_N50 <- data$Clean_read_length_N50 / 1000
    data$Clean_mean_read_length <- data$Clean_mean_read_length / 1000

    data$Clean_total_base <- data$Clean_total_base / 3000000000

    ggplot(data, aes(x=Clean_mean_read_length)) +  #Clean_read_length_N50
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3",  binwidth=0.5) +
      xlab("Mean length (kb)") + ylab("Number") +  
    theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text = element_text( size=rel(1.2 ))) + 
      theme(axis.title = element_text( size=rel(1.5 ))) 

    pdf1 <- paste(outPrefix, "_length_hist.pdf", sep="")
    ggsave(pdf1, width=width, height=height)


    ggplot(data, aes(x=Clean_read_quality)) + 
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3", binwidth=0.1) +
        xlab("Mean read quality") + ylab("Number") +  
    theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text = element_text( size=rel(1.2 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))

    pdf2 <- paste(outPrefix, "_quality_hist.pdf", sep="")
    ggsave(pdf2, width=width, height=height)


    ggplot(data, aes(x=Clean_total_base)) + 
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3", binwidth=1) +
        xlab("Read depth") + ylab("Number") +  
    theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text = element_text( size=rel(1.2 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))

    pdf3 <- paste(outPrefix, "_depth_hist.pdf", sep="")
    ggsave(pdf3, width=width, height=height)



}


Quality_hist(argv$input, argv$outPrefix, argv$width, argv$height)



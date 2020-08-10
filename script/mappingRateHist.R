#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(readr)


#usage: Rscript ~/github/NanoHub/script/mappingRateHist.R --input Samples_mapping_summary.xls --pdf mapping_rate.pdf --width 5 --height 4

arg <- arg_parser('hist of summary for quality control.')
arg <- add_argument(arg, '--input', help='In file with read quality.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



Quality_hist <- function(infile, pdf_file, width, height){
    ### old infile
    # Sample  Total_reads     Mapped_reads    Total_bases     Mapped_bases    Mapped_base_rate(%)     Error_rate(%)   Mapped_read_rate(%)
    # CN001   5533458 5529814 43989664936     43983153607     99.99   11.83   99.93
    # CN002   6553456 6548325 47377543973     47368725630     99.98   11.82   99.92




    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)
    colnames(data) <- c("Sample", "Total_reads", "Mapped_reads", "Total_bases", "Mapped_bases", "Mapped_base_rate", "Error_rate", "Mapped_read_rate")



    ggplot(data, aes(x=Mapped_read_rate)) +  #Clean_read_length_N50
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3",  binwidth=0.02) +
      xlab("Read mapping rate (%)") + ylab("Number") +  
      theme_bw()  + 
      theme(axis.text = element_text( size=rel(1.2))) + 
      theme(axis.title = element_text( size=rel(1.5 ))) 

    ggsave(pdf_file, width=width, height=height)





}


Quality_hist(argv$input, argv$pdf, argv$width, argv$height)


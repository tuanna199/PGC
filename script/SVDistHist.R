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
    # Sample  centromere      length  original        read3
    # CN001   18463   18638   19970   18654
    # CN002   18628   18812   19951   18831
    # CN003   18737   18916   19975   18935





    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)




    ggplot(data, aes(x=centromere)) +  #Clean_read_length_N50
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3",  binwidth=200) +
      xlab("SV number") + ylab("Individuals") +  
      theme_bw()  + 
      theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
      xlim(15000, 24000) +
      theme(axis.text = element_text( size=rel(1.2))) + 
      theme(axis.title = element_text( size=rel(1.5 ))) 

    ggsave(pdf_file, width=width, height=height)





}


Quality_hist(argv$input, argv$pdf, argv$width, argv$height)


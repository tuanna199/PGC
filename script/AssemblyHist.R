#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(readr)

arg <- arg_parser('hist of summary for quality control.')
arg <- add_argument(arg, '--input', help='In file with read quality.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


#usage: Rscript ~/github/NanoHub/script/QualityHist.R --input Samples_quality_summary.xls --pdf1 temp1.pdf --pdf2 temp2.pdf --width 6 --height 4


Quality_hist <- function(infile, pdf_file, width, height){
    #infile
    # Sample  No. of contigs  Contig N50 (bp) Total bases (bp)    Genome sequences coverage (%) Protein-coding genes coverage (%)
    # CN001   3,279   6,827,280   2,737,667,863   93.6    90.9 
    # CN002   3,307   6,567,130   2,736,956,236   93.5    91.1 

    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)
    colnames(data) <- c("Sample", "ContigNum", "ContigN50", "AssemblyTotal", "GenomeCov", "ProteinCov")
    data$AssemblyTotal <- data$AssemblyTotal / 1000000000

    ggplot(data, aes(x=AssemblyTotal)) +  #Clean_read_length_N50
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3",  binwidth=0.005) +
      xlab("Assembly length (Gb)") + ylab("Sample count") +  
      theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
      xlim(2.65, 2.80) +
      theme(axis.text = element_text( size=rel(1.1 ))) + 
      theme(axis.title = element_text( size=rel(1.3 ))) 

    ggsave(pdf_file, width=width, height=height)



}


Quality_hist(argv$input, argv$pdf, argv$width, argv$height)



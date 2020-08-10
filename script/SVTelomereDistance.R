#!/usr/bin/Rscript
library(ggplot2)
library(readr) 
library(argparser)


#usage: Rscript ~/github/NanoHub/script/SVTelomereDistance.R --input Sample_common_to_terminal_distance.xls --pdf temp.pdf --width 6 --height 4

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



TelomereDistance <- function(infile, pdf_file, width, height){
    height <- as.numeric(height)
    width <- as.numeric(width)


    data <- read_tsv(infile)

    data$Start <- data$Pos 


    # p<-ggplot(data,aes(x=Pos, y=Number, fill="blue")) + 
      # geom_point() +
    #   theme_bw() +
    #   ylim(0, 100) + # 50, 300
    #   xlab("Meta-chromosome (Mean of Autosomes)") +
    #   ylab("SV counts per 100 kbp") +
    #   geom_smooth(method = loess)
    # p

  p<-ggplot(data, aes(x=Pos)) + 
  #  geom_freqpoly(binwidth = 0.01, fill="dodgerblue3")
    theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
    geom_histogram(binwidth = 0.01, fill="dodgerblue3") +
    xlab("Meta-chromosome") + #(Mean of Autosomes)
    ylab("SV count") +
    theme(axis.text = element_text( size=rel(1.5 ))) + 
    theme(axis.title = element_text( size=rel(1.8 )))  +
    theme(plot.title = element_text(hjust = 0.5, size=rel(1.5)))

  p
    

  ggsave(pdf_file, width=width, height=height)
}

TelomereDistance(argv$input, argv$pdf, argv$width, argv$height)



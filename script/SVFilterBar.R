#!/usr/bin/Rscript
library(argparser)
library(ggplot2)
library(readr)
library(reshape2)

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


#usage: Rscript ~/github/NanoHub/script/SVFilterBar.R --input SVCall/FiltStats/Samples_types_summary.xls --pdf SVCall/FiltStats/Samples_types_summary.pdf --width 4 --height 4


SV_Filter_Bar <- function(in_file, pdf_file, width, height) {

  width <- as.numeric(width)
  height <- as.numeric(height)


  data <- read_tsv(in_file)
  dataMean <- apply(data[, 2:ncol(data)], 2, mean)
  dataSD <- apply(data[, 2:ncol(data)], 2, sd)

  newData <- data.frame(Category=c("Region", "Length", "Original", "Depth"), Mean = dataMean, Sd = dataSD)
  newData$Category <- factor(newData$Category, order=TRUE, levels=c("Original", "Depth", "Length", "Region"))


  p <- ggplot(newData, aes(x=Category, y=Mean)) + 
    geom_bar(stat="identity", fill="cornflowerblue", position=position_dodge(), width=0.7) +
    geom_errorbar(aes(ymin= Mean - Sd, ymax=Mean + Sd), width=.2, position=position_dodge(.7), colour="black", size=0.5) +
    # theme_bw()  + 
    theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Filtering steps") + ylab("SV number") +
    coord_cartesian(ylim=c(15000, 21000)) +
    theme(axis.text = element_text( size=rel(1.2 ))) + 
    theme(axis.title = element_text( size=rel(1.5 )))


  p

  ggsave(pdf_file, width=width, height=height)

}


SV_Filter_Bar(argv$input, argv$pdf, argv$width, argv$height)


#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(readr)
library(reshape2)


arg <- arg_parser('Bar plot for enrichment.')
arg <- add_argument(arg, '--input', help='The input file with enriment.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
arg <- add_argument(arg, '--termNum', help="Selected term number")
argv <- parse_args(arg)


#usage: Rscript ~/github/NanoHub/script/enrichBarPlot.R --input GWAS_Catalog_2019..enrichr.reports.txt --pdf GWAS_Catalog_2019..enrichr.reports.pdf --width 8 --height 4 --termNum 12


enrichment_bar_plot <- function(in_file, pdf_file, width, height, termNum){
    width <- as.numeric(width)
    height <- as.numeric(height)
    termNum <- as.numeric(termNum)

    data <- read_tsv(in_file)
    subdata <- data[1:termNum, 1:5]
    colnames(subdata) <- c("GeneSets", "Term", "Overlap", "Pvalue", "AdjustedPvalue")

    subdata$AdjustedPvalue <- -log10(subdata$AdjustedPvalue)

    print(subdata)

    barPlot <- ggplot(subdata, aes(y=AdjustedPvalue, x=reorder(Term, AdjustedPvalue, reverse=TRUE), fill="orange")) +
        geom_bar(stat="identity") +
        theme_bw() +
        scale_x_discrete(labels = function(y) lapply(strwrap(y, width = 40, simplify = FALSE), paste, collapse="\n")) +
        ylab(expression(-log[10]("Adjusted P-value"))) + xlab("Term") +
        coord_flip() +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +  ##panel.border=element_blank(),
        theme(axis.text.x=element_blank()) + #axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        theme(axis.text = element_text( size=rel(1.3 ))) +
        theme(axis.title = element_text( size=rel(1.4 ))) +
        guides(fill = FALSE)

      barPlot
      
      
      ggsave(pdf_file, width=width, height=height)

}


enrichment_bar_plot(argv$input, argv$pdf, argv$width, argv$height, argv$termNum)

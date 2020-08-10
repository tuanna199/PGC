#!/usr/bin/Rscript
library(argparser)
library(ggplot2)
library(readr)

arg <- arg_parser('Line plot of SV number for shuffled sample number.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


length_box <- function(in_file, pdf_file, width, height) {
    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)
  
    colors <- c("mediumseagreen", "dodgerblue3", "gold2", "tomato2")
    data$Category <- factor(data$Category, order=TRUE, levels=c("Singleton", "Rare",  "Low", "Common"))

    data$Distance <- log2(data$Distance)
    plot_box <- ggplot(data, aes(x=Category, y=Distance, fill=Category)) +
        geom_boxplot()  + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        scale_fill_manual(values=colors) +
        theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0.5)) + 
        theme(plot.title = element_blank()) + #, axis.title.x = element_blank(), axis.title.y = element_blank()) + 
        guides(fill=FALSE) + 
        # xlab("") + ylab(expression(paste("SV length (log"[2], ")"))) +
        xlab("") + ylab(expression(paste("LoF SV length (log"[2], ")"))) +
        theme(plot.margin = margin(0.3, 0.3, 0.5, 0.5, "cm")) + 
        theme(axis.text = element_text( size=rel(1.5 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))
    
    plot_box

    ggsave(pdf_file, width=width, height=height)

}

length_box(argv$input, argv$pdf, argv$width, argv$height)




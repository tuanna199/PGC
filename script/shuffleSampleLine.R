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


shuffle_sample_line <- function(in_file, pdf_file, width, height){
    # in_file:
    # Sample  Category    Number
    # 100 Singleton   14039
    # 100 Rare    16525
    # 100 Low 13868
    # 100 Common  32781
    # 100 All 77214

    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)
    data$Number <- data$Number / 1000

    # colors <- c("black", "tomato2", "gold2", "dodgerblue3", "mediumseagreen")
    # data$Category <- factor(data$Category, order=TRUE, levels=c("All", "Common", "Low", "Rare", "Singleton"))

    colors <- c("tomato2", "gold2", "dodgerblue3", "mediumseagreen")
    data$Category <- factor(data$Category, order=TRUE, levels=c("Common", "Low", "Rare", "Singleton"))


    plot <- ggplot(data, aes(x=Sample, y=Number, color=Category))+
        geom_line(size=rel(1.2)) +
        theme_bw() +
        xlab("Samples") + ylab("SV number (k)") +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x=element_blank()) + #axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_color_manual(values=colors)  + 
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,0.8,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())

    plot
    ggsave(pdf_file, width=width, height=height)

}

shuffle_sample_line(argv$input, argv$pdf, argv$width, argv$height)

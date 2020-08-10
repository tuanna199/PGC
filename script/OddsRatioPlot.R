#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


arg <- arg_parser('Bar plot.')
arg <- add_argument(arg, '--input', help='The file with odds ratio of annotation.')
arg <- add_argument(arg, '--pdf', help='output box file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


odds_ratio_barplot <- function(in_file, pdf_file, width, height){
    # in_file:
    # Category        CDS     Enhancer        Intron  NC_Exon Promoter        UTR5/3  Up/DownStream   Other   All
    # All     3147    42657   63852   2653    2972    2552    7996    45887   132356
    # Common  620     10366   15478   461     538     434     1742    11634   32782
    # Low     290     4420    6581    242     250     223     819     5298    14362
    # Rare    597     9032    14182   542     534     502     1726    10029   28951
    # Singleton       1640    18839   27611   1408    1650    1393    3709    18926   56261




    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)
    ### delete column "Up/DownStream"
    #data1 <- data[ , -which(names(data) %in% c("Up/DownStream", "Enhancer", "NC_Exon"))]
    colnames(data) <- c("Category", "CDS", "Intron", "Promoter", "UTR", "Intergenic") 
    melt_dt <- reshape2::melt(data)
    colnames(melt_dt) <- c("Category", "Type", "Value")
    melt_dt$Value <- log2(melt_dt$Value)

    colors <- c("mediumseagreen", "dodgerblue3", "gold2", "tomato2")
    melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=c("Singleton", "Rare", "Low", "Common"))

    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("CDS", "Exon", "UTR5", "UTR3", "UpStream", "DownStream", "Intron", "Intergenic"))

    melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("Promoter", "UTR", "CDS",  "Intron", "Intergenic"))
 
    barPlot <- ggplot(melt_dt, aes(x=Type, y=Value, fill=Category)) +
        # theme_bw() +
        geom_bar(stat="identity", position=position_dodge()) +
        xlab("") + ylab(expression(log[2]("Odds ratio"))) +
        ylim(-1.0, 1.0) +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(angle=45, hjust= 1, vjust=1)) + #axis.ticks.x=element_blank()
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(axis.text = element_text( size=rel(1.1 ))) +
        theme(axis.title = element_text( size=rel(1.3 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())
        # theme(legend.position="none")
      barPlot
      ggsave(pdf_file, width=width, height=height)

}

odds_ratio_barplot(argv$input, argv$pdf, argv$width, argv$height)

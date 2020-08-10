#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


#usage: Rscript /home/wuzhikun/github/NanoHub/script/CategoryBox.R --input /home/wuzhikun/Project/Population/population/genotype/Sample_SV_type_heter2homo.txt --pdf /home/wuzhikun/Project/Population/population/genotype/Sample_SV_type_heter2homo.pdf --width 4 --height 4

arg <- arg_parser('Box plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output box file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


category_heter2homo_box <- function(in_file, pdf_file, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)

    # ## for type 
    # melt_dt <- reshape2::melt(data[1:nrow(data)-1, ])
    # colnames(melt_dt) <- c("Category", "Sample", "Value")
    # melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=c("DEL", "INS", "DUP", "INV"))
    # colors <-c("limegreen", "royalblue1", "gold2", "tomato")

    ### for category
    melt_dt <- reshape2::melt(data)
    colnames(melt_dt) <- c("Category", "Sample", "Value")
    colors <- c( "mediumseagreen",  "dodgerblue3","gold2",  "tomato2")
    melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=c("Singleton", "Rare", "Low", "Common"))


    # ### for feature
    # melt_dt <- reshape2::melt(data)
    # colnames(melt_dt) <- c("Category", "Sample", "Value")
    # colors <- c("aquamarine3", "royalblue2", "purple2", "orangered3", "orange", "black", "blue")
    # melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=c("Enhancer", "Promoter", "UTR5/3", "CDS", "NC_Exon", "Intron", "Up/DownStream"))


    BoxPlot <- ggplot(melt_dt, aes(x=Category, y=Value, fill=Category)) +
            geom_boxplot() + 
            # theme_bw() +
            ylim(0, 15) + #12  #15
            xlab("") + ylab("Heter/Homo ratio") +
            theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
            # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
            scale_fill_manual(values=colors)  + 
            theme(axis.text = element_text( size=rel(1.2))) +
            theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
            theme(axis.title = element_text( size=rel(1.4 )))  +
            # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1)) + #angle = 90, hjust = 1
            # theme(axis.title = element_text( size=rel(0.8)))  +
            theme(plot.margin = margin(1,1,1,1, "cm")) +
            guides(fill=FALSE)

    BoxPlot
    ggsave(pdf_file, width=width, height=height)

}

category_heter2homo_box(argv$input, argv$pdf, argv$width, argv$height)

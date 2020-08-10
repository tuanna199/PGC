#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/TrioWGS/script/SVTypeBar.R --input M625-0_summary.xls --pdf M625-0_summary.pdf --width 6 --height 4


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)

stackBar <- function(in_file, pdf_file, height, width) {
    height <- as.numeric(height)
    width <- as.numeric(width)

    SV_summary <- read_tsv(in_file)
    ### delete last column with total number of SV type
    # SV_summary <- SV_summary[, -ncol(SV_summary)]

    melt_dt <- reshape2::melt(SV_summary)
    colnames(melt_dt) <- c("Sample", "Type", "Number")


    # melt_dt <- melt_dt[which(melt_dt$Type %in% c("DEL", "INS", "TRA", "DUP", "INV", "INVDUP")), ]
    # ### order of  x axis 
    # # Type_order <- factor(melt_dt$Type, order=TRUE, levels=c("DEL", "INS", "DUP", "INV", "TRA", "INVDUP"))
    # # colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")

    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("DEL", "INS", "TRA", "DUP", "INV",  "INVDUP"))
    # colors <-c("limegreen", "royalblue1", "purple2", "tomato", "gold2",  "gray60")

    melt_dt <- melt_dt[which(melt_dt$Type %in% c("DEL", "INS", "DUP", "INV")), ]
    melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("DEL", "INS", "DUP", "INV"))
    colors <-c("limegreen", "royalblue1", "gold2", "tomato")
    
    melt_dt$Number <- log10(melt_dt$Number)

    print(melt_dt)

    barPlot <- ggplot(melt_dt, aes(x=Type, y=Number, fill=Type)) +
        theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("SV type") + ylab(expression(paste("Number (log"[10], ")"))) +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x=element_blank()) + #axis.ticks=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(axis.text = element_text( size=rel(1.3 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        guides(fill = FALSE)
        # theme(legend.position = "bottom") +
        # theme(legend.text=element_text(size=rel(0.5))) +
        # theme(legend.key.width = unit(0.5, "cm"))
        # theme(legend.title = element_blank()) +
        # theme(legend.position="none")

    barPlot
    ggsave(pdf_file, width=width, height=height)
}

stackBar(argv$input, argv$pdf, argv$height, argv$width)




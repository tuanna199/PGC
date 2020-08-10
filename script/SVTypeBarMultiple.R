#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


#usage: Rscript ~/github/NanoHub/script/SVTypeBarMultiple.R --input Samples_SV_summary.xls --hist hist.pdf --box box.pdf --width 8 --height 4

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--hist', help='output hist file with pdf format.')
arg <- add_argument(arg, '--box', help='output box file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
arg <- add_argument(arg, "--maxNum", default=5, help="The max number for box plot.")
argv <- parse_args(arg)




boxPlot <- function(melt_dt, pdf_file, colors, maxNum){
    # colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")
    colors <-c("limegreen", "royalblue1", "gold2", "tomato",  "purple2",  "black")

    maxNum <- as.numeric(maxNum)

    melt_dt$Value <- log10(melt_dt$Value)
    # melt_dt$Type <- factor(melt_dt$Type, levels=c("DEL", "INS", "DUP", "INV", "TRA", "INVDUP"))
    melt_dt$Type <- factor(melt_dt$Type, levels=c("DEL", "INS", "DUP", "INV"))
    plot_box <- ggplot(melt_dt, aes(x=Type, y=Value, fill=Type)) + # x=reorder(Type, -Value, median),
        geom_boxplot()  + 
        # theme_bw() + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        scale_fill_manual(values=colors) +
        theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0)) + 
        ylim(1, max(maxNum)) + ### maxNum = 15000
        theme(plot.title = element_blank()) + #, axis.title.x = element_blank(), axis.title.y = element_blank()) + 
        guides(fill=FALSE) + 
        xlab("Type") + ylab(expression(paste("SV Number (log"[10], ")"))) +
        # scale_y_continuous(limits=c(0,0.8)) + 
        # theme(plot.margin = margin(0.3, 0.3, 0.5, 0.5, "cm")) +
        theme(axis.text = element_text( size=rel(1.1 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())

    plot_box
    ggsave(pdf_file, width=4, height=4)
}


 
multipleStackBar <- function(melt_dt, pdf_file, colors, height, width){
    # melt_dt$Type <- factor(melt_dt$Type, levels=c("DEL", "INS", "DUP", "INV", "TRA", "INVDUP"))
    melt_dt$Type <- factor(melt_dt$Type, levels=c("DEL", "INS", "DUP", "INV"))
    
    colors <-c("limegreen", "royalblue1", "gold2", "tomato",  "purple2",  "black")

    p <- ggplot(melt_dt, aes(x=reorder(Sample, -Value), y=Value, fill=Type)) + 
        # theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("Samples") + ylab("Number") +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        # theme(axis.ticks = element_line(size=rel(0.5))) + 
        # theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 0.5, size=rel(0.5))) +
        theme(axis.ticks.x = element_blank()) + 
        theme(axis.text.x = element_blank()) +
        theme(plot.margin = margin(0.3,0.5,0.5,0.3, "cm")) +
        scale_fill_manual(values=colors) +
        theme(axis.text.y = element_text( size=12 )) + #size=rel(1.5)
        theme(axis.title= element_text(size=18)) +
        theme(legend.box = "horizontal", legend.position = "top") + #legend.position = c(0.50, 0.85)
        # theme(legend.position = "right") + # bottom
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.5, "cm")) +
        theme(legend.title = element_blank()) 

    p

    height <- as.numeric(height)
    width <- as.numeric(width)
    ggsave(pdf_file, width=width, height=height)
}

 

# colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")



in_file <- argv$input
summary_data <- read_tsv(in_file)
# summary_data <- summary_data[, 1:7]
summary_data <- summary_data[, 1:5]
melt_dt <- reshape2::melt(summary_data)
colnames(melt_dt) <- c("Sample", "Type", "Value")



boxPlot(melt_dt, argv$box, colors, argv$maxNum)
multipleStackBar(melt_dt, argv$hist, colors, argv$height, argv$width)

 

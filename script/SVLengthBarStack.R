#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/TrioWGS/script/SVTypeBarStack.R --input M628-2_dist_summary.xls --pdf M628-2_dist_summary.pdf --width 7 --height 5 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
arg <- add_argument(arg, '--chromosome', default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X", help="The target chromosome name with order.")
argv <- parse_args(arg)

stackBar <- function(in_file, pdf_file, height, width, chromosome="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X") {
    ### in_file 
    # Chr     DEL     INS     DUP     INV     TRA     Total
    # 1       1693    843     25      28      214     2803
    # 2       1769    677     23      25      127     2621


    height <- as.numeric(height)
    width <- as.numeric(width)

    SV_summary <- read_tsv(in_file)

    ### delete the last line with summary
    # SV_summary <- SV_summary[1:nrow(SV_summary)-1, ]

    SV_Len <- apply(SV_summary[, 2:ncol(SV_summary)], 2, function(x) x/1000000)
    newSummary <- cbind(SV_summary[, 1], SV_Len)

    # ### delete last column with total number of SV type
    # newSummary <- newSummary[, -ncol(newSummary)]




    melt_dt <- reshape2::melt(newSummary)
    colnames(melt_dt) <- c("Chr", "Type", "Number")

    ### select types
    melt_dt <- melt_dt[which(melt_dt$Type %in% c("DEL", "INS", "DUP", "INV")), ]

    ### select chromosomes
    chrs <- strsplit(chromosome, ",")
    chrs <- as.vector(chrs[[1]])
    melt_dt <- melt_dt[which(melt_dt$Chr %in% chrs), ]

    melt_dt$Chr=factor(melt_dt$Chr,levels=chrs)

    ### order of  x axis 
    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("INV", "DUP", "DEL", "INS")) 
    
    # ### color of types
    # colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "darkolivegreen3", 
    #     "mediumorchid3",  "gold3", "deeppink3", "green4", "lightblue3", 
    #     "royalblue1",  "thistle1", "skyblue", "purple2", "orangered")
        
    colors <-c("tomato",  "gold2", "royalblue1", "limegreen") 
    melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("INV", "DUP", "INS", "DEL")) 

    barPlot <- ggplot(melt_dt, aes(x=Chr, y=Number, fill=Type)) +
        # theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("Chromosome") + ylab("SV length (Mb)") +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x=element_blank()) + #, axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + # hjust = 1
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
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

stackBar(argv$input, argv$pdf, argv$height, argv$width, chromosome=argv$chromosome)




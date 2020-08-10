#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)
library(grid)

#usage: Rscript ~/github/NanoHub/script/TypeFreqBar.R --input /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype_frequency_stat.xls --pdf /home/wuzhikun/Project/Population/population/genotype/Sample_SV_genotype_frequency_stat.png  --width 6 --height 4


arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)

stackBar <- function(in_file, pdf_file, height, width) {
    ### in_file 
    # Chr     DEL     INS     DUP     INV     TRA     Total
    # 1       1693    843     25      28      214     2803
    # 2       1769    677     23      25      127     2621


    height <- as.numeric(height)
    width <- as.numeric(width)

    SV_summary <- read_tsv(in_file)
    ### delete last column with total number of SV type
    # SV_summary <- SV_summary[, -ncol(SV_summary)]
    melt_dt <- reshape2::melt(SV_summary)
    colnames(melt_dt) <- c("Type", "Region", "Number")

    ### order of  x axis 
    # colors <-c("black","purple2", "gold2", "tomato", "royalblue1", "limegreen")
    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("INVDUP","TRA", "INV", "DUP", "INS", "DEL")) 


    # melt_dt <- melt_dt[which(melt_dt$Type %in% c("DEL", "INS", "DUP", "INV")), ]

    # colors <-c( "gold2", "tomato", "purple2", "royalblue1", "limegreen") 
    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("INV", "DUP", "TRA", "INS", "DEL")) 
    colors <- c( "slategray3", "mediumseagreen",  "dodgerblue3", "gold2", "orchid2",  "tomato2")

    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("DEL", "INS", "DUP", "INV")) 
    melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels= c("dbVar", "DGV", "LRS15", "WGS911",  "WGS17795", "gnomAD")) 
    melt_dt$Region <- factor(melt_dt$Region, order=TRUE, levels=c("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0", "1.0"))

    barPlot1 <- ggplot(melt_dt, aes(x=Region, y=Number, fill=Type)) +
        theme_bw() +
        geom_bar(stat="identity", position="dodge", width=0.7,) +
        # xlab("Chromosome") + ylab("Number") +
        # theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank(), panel.border=element_blank()) +  
        theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
        # theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) + 
        # theme(axis.text.y = element_text( size=rel(1.2 ))) +
        # theme(axis.title.y = element_text( size=rel(1.4 )))  +
        scale_fill_manual(values=colors) +
        # theme(axis.text = element_text( size=rel(1.2 ))) +
        # theme(axis.title = element_text( size=rel(1.4 )))  +
        # theme(plot.margin = margin(1,1,1,1, "cm")) 
        # theme(legend.position = "right") +
        # theme(legend.text=element_text(size=rel(0.7))) +
        # theme(legend.key.width = unit(0.4, "cm")) +
        # theme(legend.key.height = unit(0.4, "cm")) +
        # theme(legend.title = element_blank())
        # theme(legend.position="none")
        xlab("Chromosome") + ylab("Number") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1)) + 
        theme(axis.text = element_text( size=rel(1.0))) +
        theme(axis.title = element_text( size=rel(1.3 )))  +
        # theme(axis.text.x=element_blank(), axis.title=element_blank(), axis.ticks.x = element_blank()) +
        theme(axis.text.y = element_text( size=rel(1.0 ))) +
        theme(legend.position = c(0.85, 0.7)) + #legend.position = "right"
        theme(legend.text=element_text(size=rel(0.9))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())

        barPlot1
        ggsave(pdf_file, width=width, height=height)



    # barPlot <- ggplot(melt_dt, aes(x=Region, y=Number, fill=Type)) +
    #     theme_bw() +
    #     geom_bar(stat="identity", position="dodge", width=0.7,) +
    #     # xlab("Chromosome") + ylab("Number") +
    #     # theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
    #     theme(panel.background=element_blank(),panel.grid.major=element_blank(),
    #           panel.grid.minor=element_blank(),plot.background=element_blank(), panel.border=element_blank()) +  
    #     theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +

    #     scale_fill_manual(values=colors) 



    # # ggsave(pdf_file, width=width, height=height)


    # #使用 coord_cartesian() 分割作图结果
    # split1 <- barPlot + coord_cartesian(ylim = c(0, 8000)) + 
    #     xlab("Chromosome") + ylab("Number") +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1)) + 
    #     theme(axis.text = element_text( size=rel(1.2 ))) +
    #     theme(axis.title = element_text( size=rel(1.4 )))  +
    #     theme(legend.position='none')
    # split2 <- barPlot + coord_cartesian(ylim = c(20000, 50000)) + 
    #     theme(axis.text.x=element_blank(), axis.title=element_blank(), axis.ticks.x = element_blank()) +
    #     theme(axis.text.y = element_text( size=rel(1.4 ))) +
    #     theme(legend.position = c(0.85, 0.7)) + #legend.position = "right"
    #     theme(legend.text=element_text(size=rel(0.9))) +
    #     theme(legend.key.width = unit(0.4, "cm")) +
    #     theme(legend.key.height = unit(0.4, "cm")) +
    #     theme(legend.title = element_blank())

    # # png(pdf_file, width = 2500, height = 1700, res = 300, units = "px")
    # pdf(pdf_file, width=width, height=height)
    # grid.newpage()
    # plot_site1 <- viewport(x = 0, y = 0, width = 1, height = 0.61, just = c("left", "bottom"))
    # plot_site2 <- viewport(x = 0.028, y = 0.61, width = 1, height = 0.37, just = c("left", "bottom"))
    # print(split1, vp = plot_site1)
    # print(split2, vp = plot_site2)
    # dev.off()

}

stackBar(argv$input, argv$pdf, argv$height, argv$width)




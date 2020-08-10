#!/usr/bin/Rscript
library(ggplot2)
library(reshape2)
library(readr) 
library(argparser)

arg <- arg_parser('Hist and pie plot for SV repeats.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--outPrefix', help='output prefix for pdf file.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



SVType_length_repeat <- function(infile, outPrefix){
    # ### infile:
    # Tag SVType  SVLength
    # 1_66288-1_66527-239.0-DEL   DEL 239.0
    # 1_67910-1_68341-431.0-DEL   DEL 431.0
    # 1_83968-1_84057-89.0-INS    INS 89.0
    # 1_88684-1_88831-147.0-DEL   DEL 147.0
    # 1_88889-1_88956-66.0-INS    INS 66.0


    # color2 <- c("gold3", "purple2")
    color2 <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")

    data <- read_tsv(infile)

    data$SVType <- factor(data$SVType, levels=c("DEL", "INS", "DUP", "INV"))

    hist1 <- ggplot(data, aes(x=SVLength)) + 
        geom_histogram(aes(y=..count.., fill=data$SVType), binwidth=10) +
        xlab("Variant Length (bp)") + ylab("Number") +  
        # theme_bw()  + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text = element_text( size=rel(1.5 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))  + 
        xlim(0, 1000) + 
        ggtitle("SV Length <= 1000 bp") + 
        theme(plot.title = element_text(hjust = 0.5, size=rel(1.5))) +
        scale_fill_manual(values=color2)  + 
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(1.2))) +
        theme(legend.key.width = unit(0.8, "cm")) +
        theme(legend.key.height = unit(0.8, "cm")) +
        theme(legend.title = element_blank())

    hist1
    file1 <-  paste(outPrefix, "_length_1000.pdf", sep="")
    ggsave(file1, width=6, height=4)






    hist2 <- ggplot(data, aes(x=SVLength)) + 
        geom_histogram(aes(y=..count.., fill=data$SVType), binwidth=100) +
        xlab("Variant Length (bp)") + ylab("Number") +  
        # theme_bw()  + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text = element_text( size=rel(1.5 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))  + 
        xlim(1000, 10000) + 
        ggtitle("SV Length > 1000 bp") + 
        theme(plot.title = element_text(hjust = 0.5, size=rel(1.5))) +
        scale_fill_manual(values=color2)  + 
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(1.2))) +
        theme(legend.key.width = unit(0.8, "cm")) +
        theme(legend.key.height = unit(0.8, "cm")) +
        theme(legend.title = element_blank())

    hist2
    file2 <-  paste(outPrefix, "_length_10000.pdf", sep="")
    ggsave(file2, width=6, height=4)




}



SVType_length_repeat(argv$input, argv$outPrefix)
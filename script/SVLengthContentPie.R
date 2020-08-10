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
### infile:
#     SVID    SVLength        SVType  DfamName        DfamType
# 10_100093729-10_100093883-120.918-INS:M546-1:10_100093777-10_100093926-149-INS  120.918 INS     SVA_C   Retroposon
# 10_100277459-10_100277518-90.25-DEL:M452-1:10_100277872-10_100278000-128-DEL    90.25   DEL     hAT-N35B_DR     DNA
# 10_100277920-10_100278046-117.875-DEL:M534-2:10_100277895-10_100278026-131-DEL  117.875 DEL     hAT-N35B_DR     DNA
# 10_100278002-10_100278154-148.5-INS:M534-2:10_100278002-10_100278154-157-INS    148.5   INS     SVA_A   Retroposon



    colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "darkolivegreen3", "mediumorchid3",  "gold3", "deeppink3", "green4", "lightblue3", "royalblue1",  "thistle1", "skyblue", "purple2", "orangered")

    # color2 <- c("gold3", "purple2")
    color2 <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "black")

    data <- read_tsv(infile)

    data$SVType <- factor(data$SVType, levels=c("DEL", "INS"))

    hist1 <- ggplot(data, aes(x=SVLength)) + 
        geom_histogram(aes(y=..count.., fill=data$SVType), binwidth=10) +
        xlab("Variant Length (bp)") + ylab("Number") +  
        theme_bw()  + 
        theme(axis.text = element_text( size=rel(1.5 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))  + 
        xlim(0, 1000) + 
        ggtitle("Variant Length <= 1000 bp") + 
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
        theme_bw()  + 
        theme(axis.text = element_text( size=rel(1.5 ))) + 
        theme(axis.title = element_text( size=rel(1.5 )))  + 
        xlim(1000, 10000) + 
        ggtitle("Variant Length > 1000 bp") + 
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





  
    data$Number <- rep(1, nrow(data))

    #agg <- aggregate(data$Number, by=list(data$DfamType, data$SVType), FUN=sum)
    #colnames(agg) <- c("DfamType", "SVType", "Number")
    
    ### INS
    subs <- subset(data, SVType=="INS")
    agg <- aggregate(subs$Number, by=list(subs$DfamType), FUN=sum)
    colnames(agg) <- c("DfamType",  "Number")
    
    ### get top 8
    subdata <-  head(agg[order(agg$Number, decreasing= T),], n = 8)
    subdata$Percentage <- subdata$Number / sum(subdata$Number)

    # Create a basic bar
    INS_pie <- ggplot(subdata, aes(x="", y= Percentage, fill=DfamType)) + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) + 
        geom_text(aes(label = paste0(round(Percentage*100), "%")), position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values=colors) +
        labs(x = NULL, y = NULL, fill = NULL, title = "") +
        theme_classic() + 
        ggtitle("Repeat contents for INS") +
        theme(plot.title = element_text(hjust = 0.5, size=rel(1.5))) +
        theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, color = "black"))


    INS_pie
    file3 <-  paste(outPrefix, "_INS_pie.pdf", sep="")
    ggsave(file3, width=6, height=4)


    ### DEL
    subs <- subset(data, SVType=="DEL")
    agg <- aggregate(subs$Number, by=list(subs$DfamType), FUN=sum)
    colnames(agg) <- c("DfamType",  "Number")
    
    ### get top 8
    subdata <-  head(agg[order(agg$Number, decreasing= T),], n = 8)
    subdata$Percentage <- subdata$Number / sum(subdata$Number)

    # Create a basic bar
    DEL_pie <- ggplot(subdata, aes(x="", y= Percentage, fill=DfamType)) + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) + 
        geom_text(aes(label = paste0(round(Percentage*100), "%")), position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values=colors) +
        labs(x = NULL, y = NULL, fill = NULL, title = "") +
        theme_classic() + 
        ggtitle("Repeat contents for DEL") +
        theme(plot.title = element_text(hjust = 0.5, size=rel(1.5))) +
        theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, color = "black"))



    DEL_pie
    file4 <-  paste(outPrefix, "_DEL_pie.pdf", sep="")
    ggsave(file4, width=6, height=4)

}



SVType_length_repeat(argv$input, argv$outPrefix)
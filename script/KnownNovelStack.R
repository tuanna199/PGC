library(argparser)
library(ggplot2)
library(readr)
library(reshape2)

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


known_novel_stack <- function(infile, pdf_file, width, height){
    #infile:
    # Category        KnownTags       NovelTags       KnownRatio      NovelRatio
    # Common  27990   4792    0.854   0.146
    # Low     7436    6926    0.518   0.482
    # Rare    10517   18434   0.363   0.637
    # Singleton       13262   42999   0.236   0.764

    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)




    # ## for different type SV
    # subdata <- data[, c(1,2,3)]
    # colnames(subdata) <- c("Category", "Known", "Novel")
    # melt_dt <- reshape2::melt(subdata)
    # melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=c("DEL","INS", "DUP", "INV"))
    # colnames(melt_dt) <- c("Category", "Type", "Value")
    # melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("Novel","Known"))
    # colors <-c( "tomato2", "dodgerblue3")


    ### for different category
    colnames(data) <- c("Category", "Known", "Novel", "KnownRatio", "NovelRatio")
    melt_dt <- reshape2::melt(data[, 1:3])
    colnames(melt_dt) <- c("Category", "Type", "Value")
    melt_dt$Value <- melt_dt$Value / 1000

    melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=c("Singleton","Rare", "Low", "Common"))
    melt_dt$Type <- factor(melt_dt$Type, order=TRUE, levels=c("Novel","Known"))
    colors <-c( "tomato2", "dodgerblue3")





    barPlot <- ggplot(melt_dt, aes(x=Category, y=Value, fill=Type)) +
        geom_bar(stat="identity",width=0.7,) +
        xlab("") + ylab("SV number (k)") +
        # ylim(0, 100) +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x=element_blank()) + #axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
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


known_novel_stack(argv$input, argv$pdf, argv$width, argv$height)

library(argparser)
library(readr)
library(ggplot2)

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


individuals_per_SV <- function(in_file, pdf_file, width, height){
    # in_file:
    # SVType  SVNumber        Individuals     Stdev
    # DEL     59735   23.263  46.283
    # DUP     3592    9.697   24.970
    # INS     52364   30.905  53.627
    # INV     2259    8.336   28.225


    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)

    colors <-c( "limegreen", "royalblue1", "gold2", "tomato" ) 
    data$SVType <- factor(data$SVType, order=TRUE, levels=c("DEL", "INS", "DUP",  "INV" )) 

    
    barPlot <- ggplot(data, aes(x=SVType, y=Individuals, fill=SVType)) +
        # theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("SV type") + ylab("Individuals per SV") +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x=element_blank()) + #axis.ticks.x=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(plot.margin = margin(1,1,1,1, "cm")) + # +
        # theme(legend.position = "right") +
        # theme(legend.text=element_text(size=rel(0.7))) +
        # theme(legend.key.width = unit(0.4, "cm")) +
        # theme(legend.key.height = unit(0.4, "cm")) +
        # theme(legend.title = element_blank())
        guides(fill=FALSE)
    barPlot
    ggsave(pdf_file, width=width, height=height)
}


individuals_per_SV(argv$input, argv$pdf, argv$width, argv$height)

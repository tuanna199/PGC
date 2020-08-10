#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(readr)

#usage: Rscript ~/github/TrioWGS/script/SVDensity.R --input M625-0_record.xls --pdf M625-0_record.pdf --width 7 --height 5


arg <- arg_parser('Density plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


SVDensity <- function(input_file, pdf_file, width, height){
    ### input_file
    # Type    SVLEN   RE
    # DEL     106     2
    # DEL     64      2
    # DEL     64      2

    ### color of types
    # colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "darkolivegreen3", 
    #     "mediumorchid3",  "gold3", "deeppink3", "green4", "lightblue3", 
    #     "royalblue1",  "thistle1", "skyblue", "purple2", "orangered")




    width <- as.numeric(width)
    height <- as.numeric(height)

    SV_density <- read_tsv(input_file)

    New_density <-  SV_density[!grepl("TRA", SV_density$Type), ]

    # New_density$logSVLEN <- log10(New_density$SVLEN)

    ### color of types
    # colors <-c("limegreen", "royalblue1","tomato", "gold2", "purple2",  "darkolivegreen3", 
    #     "mediumorchid3",  "gold3", "deeppink3", "green4", "lightblue3", 
    #     "royalblue1",  "thistle1", "skyblue", "purple2", "orangered")

    ### DEL, DUP, INS, INV, INVDUP
	colors <-c( "tomato", "gold2", "royalblue1", "limegreen")
    New_density$Type <- factor(New_density$Type, order=TRUE, levels=c( "INV", "DUP", "INS", "DEL"))

    # colors <-c("gray60", "gold2", "tomato", "royalblue1", "limegreen")
    # New_density$Type <- factor(New_density$Type, order=TRUE, levels=c("INVDUP", "INV", "DUP", "INS", "DEL"))


    density <- ggplot(New_density, aes(x=log10(SVLEN), fill=Type)) +
        # theme_bw() +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        # scale_x_continuous(trans='log10') +
        geom_density(alpha=.5) +
        xlab(expression(log[10]("bp"))) + ylab("Density") +
        xlim(1.5, 4) +
        theme(panel.background=element_blank(),panel.grid.major=element_blank(),  panel.grid.minor=element_blank(),plot.background=element_blank()) +  
        theme(axis.text = element_text( size=rel(1.2 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        scale_fill_manual(values=colors) +
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(legend.position = "right") +
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.4, "cm")) +
        theme(legend.key.height = unit(0.4, "cm")) +
        theme(legend.title = element_blank())


    density
    ggsave(pdf_file, width=width, height=height)
}

SVDensity(argv$input, argv$pdf, argv$width, argv$height)





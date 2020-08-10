#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


#usage: Rscript ~/github/NanoHub/script/LengthHeterBox.R --input temp_del.txt --pdf temp_del.pdf --wdith 4 --height 4

arg <- arg_parser('Box plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output box file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


Length_heter2homo_box <- function(in_file, pdf_file, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(in_file)
    melt_dt <- reshape2::melt(data[1:nrow(data)-1, ])
    colnames(melt_dt) <- c("Category", "Sample", "Value")

    regions <- c("50-100", "100-200", "200-500", "500-1000", "1000-2000", "2000-5000", ">5000")

    melt_dt$Category <- factor(melt_dt$Category, order=TRUE, levels=regions)
    print(head(melt_dt))

    BoxPlot <- ggplot(melt_dt, aes(x=Category, y=Value, fill=Category)) +
            geom_boxplot() + 
            ylim(0, 12) +
            xlab("SV length") + ylab("Heter/Homo ratio") +
            theme(panel.background=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +  ##panel.border=element_blank(),
            theme(axis.text = element_text( size=rel(1.2 ))) +
            theme(axis.title = element_text( size=rel(1.4 )))  +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + #angle = 90, hjust = 1
            # theme(plot.margin = margin(1,1,0.5,0, "cm")) +
            # scale_fill_manual(values=colors)  + 
            theme(plot.margin = margin(1,1,1,1, "cm")) +
            guides(fill=FALSE) +
            theme_bw()
    BoxPlot
    ggsave(pdf_file, width=width, height=height)

}

Length_heter2homo_box(argv$input, argv$pdf, argv$width, argv$height)

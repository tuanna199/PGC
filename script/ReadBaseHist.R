#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(readr)


#usage: Rscript ~/github/NanoHub/script/ReadBaseHist.R --input Samples_quality_summary.xls --pdf Samples_quality_summary_hist.pdf  --width 8 --height 4


arg <- arg_parser('Stack bar plot for total read bases.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, "--column", help='The target column.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


histPlot <- function(in_file, pdf_file, width, height){ #column, 
    read_summary <- read_tsv(in_file)

    ### bases to Gb
    read_summary$Clean_total_base <- read_summary$Clean_total_base / 1000000000

    # read_summary$Clean_total_base <- read_summary[, column] / 1000000000

    # print(head(read_summary))

    histPlot <- ggplot(read_summary, aes(x=reorder(Sample, -Clean_total_base), y=Clean_total_base)) +
        theme_bw() +
        geom_bar(stat="identity",width=0.6, color="royalblue3", fill="royalblue3") +
        # geom_hline(yintercept = 40, color="gray60", size=0.2) + #40
        xlab("Samples") + ylab("Read bases (Gbp)") +
        theme(panel.background=element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
        theme(axis.ticks = element_line(size=rel(0.5))) + #,axis.ticks=element_blank()
        theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 0.5, size=rel(0.6))) +
        theme(axis.text.y = element_text( size=12 )) + #size=rel(1.5)
        theme(plot.margin = margin(0.3,0.5,0.5,0.3, "cm")) +
        theme(axis.title= element_text(size=18)) +
        guides(fill = FALSE)

    histPlot
    width <- as.numeric(width)
    height <- as.numeric(height)

    ggsave(pdf_file, width=width, height=height)

}

histPlot(argv$input, argv$pdf, argv$width, argv$height)



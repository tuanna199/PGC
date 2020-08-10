library(ggplot2)
library(reshape2)
library(readr) 
library(muStat)
library(argparser)

#usage: Rscript ~/github/NanoHub/script/ErorrRateBar.R --input Samples_error_rate_stats.xls --pdf temp.pdf --width 4  --height 5

arg <- arg_parser('Box plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output pdf file.')
arg <- add_argument(arg, '--pdf2', help='output pdf file 2.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



Error_rate <- function(infile, pdf_file, pdf_file2, width, height){

  # #infile
  # Sample  E_ins(%)        E_del(%)        E_Mis(%)        E_total(%)      Accuracy(%)
  # M416-0  4.07    6.15    5.03    15.25   91.90
  # M416-1  3.92    6.35    5.04    15.32   90.94
  # M416-2  3.48    6.83    4.80    15.12   90.87

  width <- as.numeric(width)
  height <- as.numeric(height)

  data <- read_tsv(infile)

  colnames(data) <- c("Sample", "E_ins", "E_del", "E_Mis", "E_total", "Accuracy", "BP_align")

  mean(data$E_ins)

  subdata <- data[, 1:(ncol(data)-2)]

  Mean_data <- apply(data[,2:ncol(subdata)], 2, mean)
  Stdev_data <- apply(data[,2:ncol(subdata)], 2, function(x) stdev(x, na.rm=TRUE, unbiased=TRUE))

  newData <- data.frame(Mean = as.vector(Mean_data), Stdev = as.vector(Stdev_data), Type=c( "Insertion", "Deletion", "Mismatch", "Total"))


  ### error rate
  p<- ggplot(newData, aes(x=Type, y=Mean, fill=Type)) + 
    geom_bar(stat="identity",
    position=position_dodge()) +
    geom_errorbar(aes(ymin=Mean-Stdev, ymax=Mean+Stdev), width=.3, position=position_dodge(.6)) +
    # theme_bw() +
    theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Error type") +
    ylab("Percentage (%)") +
    theme(axis.text = element_text( size=rel(1.2 ))) + 
    theme(axis.title = element_text( size=rel(1.5 )))  +
    theme(plot.title = element_text(hjust = 0.5, size=rel(1.2))) +
    guides(fill=FALSE) 
    # theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))

  p

  ggsave(pdf_file, width=width, height=height)



  ### mapping rate

  p2 <-  ggplot(data, aes(x=BP_align)) +  #Clean_read_length_N50
        geom_histogram(aes(y=..count..), fill="cornflowerblue", color="dodgerblue3",  binwidth=0.5) +
      xlab("Base mapping rate (%)") + ylab("Number") +  
      # theme_bw()  + 
      theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
      xlim(88.5, 97.5) + 
      theme(axis.text = element_text( size=rel(1.2 ))) + 
      theme(axis.title = element_text( size=rel(1.5 ))) 
      # theme(plot.margin = margin(0.2, 0.3, 0.2, 0.2, "cm"))

  p2
  ggsave(pdf_file2, width=width, height=height)

}


Error_rate(argv$input, argv$pdf, argv$pdf2, argv$width, argv$height)




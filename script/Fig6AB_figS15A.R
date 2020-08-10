

########   Figure6AB & FigureS15A   ########

SV_len_hist_plot <- function(dat_path,pdf_file,xlab="Length (Mb)"){

  library(ggplot2)
  
  df <- read.csv(dat_path,sep="\t",check.names=FALSE)
  df$Length = df$sum_len/1000000

  p1 <- ggplot(df, aes(x = Length)) +
    geom_histogram(binwidth = 0.2, fill="#65A0D7") + 
    labs(x=xlab,y="Individuals") +
    theme(axis.text = element_text(size=rel(1)),axis.title = element_text( size=rel(1.2))) +
    theme_classic() +
    theme(plot.margin = unit(c(5.5,10.5,5.5,5.5),"pt"))

  ggsave(pdf_file, width=4, height=4)

}

SV_len_hist_plot("D:/SYSU/cnSV/FigDat/6A.txt","D:/SYSU/cnSV/Fig/Fig6A.pdf",xlab="Assembly length (Gb)")
SV_len_hist_plot("D:/SYSU/cnSV/FigDat/6B.txt","D:/SYSU/cnSV/Fig/Fig6B.pdf")
SV_len_hist_plot("D:/SYSU/cnSV/FigDat/S15A_similar_seqkit_405_hsat23_alpha.txt","D:/SYSU/cnSV/Fig/FigS15A.pdf")


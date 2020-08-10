

###########   Figure6C & FigureS15 -- histogram with subplot    #############

NRS_histogram <- function(dat1,dat2,xrange_main=c(0,12000),xrange_sub=c(10,80),pdf_file){

  library(ggplot2)
  library(grid)

  theme_set(theme_classic())

  p1 <- ggplot(dat1, aes(x = Length)) +
    geom_histogram(binwidth = 100,fill="#65A0D7") + 
    labs(x="NRS length (bp)",y="Sequences count") +
    xlim(800,12000) +
    theme(axis.text = element_text(size=rel(1.4)),axis.title = element_text( size=rel(1.4))) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(plot.margin = unit(c(5.5,12.5,5.5,5.5),"pt"))

  p2 <- ggplot(dat2,aes(x = Length_kb)) +
    geom_histogram(binwidth = 0.1,fill="#65A0D7") +
    labs(x="NRS length (kb)",y="Sequences count") +
    xlim(10,60) +
    theme_bw() +
    theme(axis.text = element_text(size=rel(1.2)),axis.title = element_text( size=rel(1.2))) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

  pdf(pdf_file,width=4,height=4) 
  print(p1)
  subplot <- viewport(width = 0.6, height = 0.6, x = 0.95,y = 0.9,just=c("right","top"))
  print(p2,vp=subplot)
  dev.off()

}

#######    Figure6C    #######
df <- read.csv("D:/SYSU/cnSV/FigDat/6C_final_3type.NRS.length.txt",
				header=FALSE,sep="\t",check.names=FALSE)

colnames(df) = c("Name","Length")
df$Length_kb = df$Length / 1000

df_main = df[df$Length_kb < 12,]
df_sub = df[df$Length_kb >= 10,]

NRS_histogram(df_main,df_sub,xrange_main=c(800,12000),xrange_sub=c(10,60),pdf_file="D:/SYSU/cnSV/Fig/fig6c.pdf")

#######    FigureS15    #######
df_sup <- read.csv("D:/SYSU/cnSV/FigDat/S15_supplementary_round4.minihit.length",
				header=FALSE,sep="\t",check.names=FALSE)

colnames(df_sup) = c("Name","Length")
df_sup$Length_kb = df_sup$Length / 1000

df_main = df_sup[df_sup$Length_kb < 12,]
df_sub = df_sup[df_sup$Length_kb >= 10,]

NRS_histogram(df_main,df_sub,xrange_main=c(800,12000),xrange_sub=c(10,60),pdf_file="D:/SYSU/cnSV/Fig/figS15.pdf")









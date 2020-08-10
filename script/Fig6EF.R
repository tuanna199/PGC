

##### Figure 6E&F #####
library(ggplot2)
theme_set(theme_classic())

dat_6e = read.csv("D:/SYSU/cnSV/FigDat/fig6e_0801.csv",check.names=F)
dat_6f = read.csv("D:/SYSU/cnSV/FigDat/fig6f_0801.csv",check.names=F)

dat_6e$Mb = dat_6e[,3] / 1000000
dat_6f$Mb = dat_6f[,3] / 1000000

dat_6e$Name = factor(dat_6e$Name,levels=dat_6e$Name)
dat_6f$Name = factor(dat_6f$Name,levels=dat_6f$Name)


plot_NRS <- function(input_dat, pdf_file, height, width){
  
  g1 <- ggplot(data = in_dat, aes(x = Name, y = Contigs)) +
        geom_bar(stat = "identity",width=.5, fill="#00ba38") + 
  	    labs(x="NRS number",y=NULL) +
        theme(plot.margin = unit(c(1,0,1,1), "mm")) +
  	    scale_y_reverse() +
        coord_flip()
  
  g2 <- ggplot(data = in_dat, aes(x = Name, y = Mb)) +
        geom_bar(stat = "identity",width=.5, fill="#f8766d") + 
  	    labs(x="NRS length (Mb)",y=NULL) +
        theme(plot.margin = unit(c(1,0,1,1), "mm")) +
        coord_flip()

  library(gridExtra)
  pdf(pdf_file,width=width,height=height)
  grid.arrange(g1,g2,ncol=2, nrow=1)
  dev.off()

}

plot_NRS(dat_6e,"fig6e.pdf",height=5,width=5)
plot_NRS(dat_6f,"fig6f.pdf",height=5,width=5)




##### Figure 4C down #####
library(ggplot2)
dat <- read.csv("D:/SYSU/cnSV/FigDat/fig4C_down.csv")

pdf("D:/SYSU/cnSV/Fig/fig4c_down.pdf",width=4,height=4)
ggplot(dat,aes(x=Genotype,y=Value, fill=Genotype)) +
  geom_boxplot() +
  xlab(NULL) + 
  scale_fill_manual(breaks = c("0/0", "0/1"),values=c("#0073C299", "#EFC00099")) +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),    ##remove grid
        panel.grid.minor=element_blank(),
		legend.position="top") +
  facet_wrap(~Trait,scales="free_y")
dev.off()

#####  Figure4D down #####

library(ggpubr)

dat <- read.csv("D:/SYSU/cnSV/FigDat/fig4D_down.csv")

dat$Sex = factor(dat$Sex,levels=c("Male","Female")) 

my_comparisons <- list( c("0/0", "0/1"), c("0/1", "1/1"), c("0/0", "1/1") )
p <- ggboxplot(dat, x = "Genotype", y = "Height",xlab = "",
          color = "black",fill="Genotype", palette = "jco",
          facet.by = "Sex", short.panel.labs = FALSE)

pdf("D:/SYSU/cnSV/Fig/fig4d_down.pdf",width=5,height=5)
p +  stat_compare_means(comparisons = my_comparisons,method = "t.test", label.y = c(190, 195, 202))
dev.off()





#####  FigureS12A  #####
library(ggplot2)
library(gridExtra)

dat <- read.csv("D:/SYSU/cnSV/FigDat/Assembly_stats.txt",sep='\t',check.names=F)

colnames(dat)[5:6] = c("Gene_sequences_coverage","Protein_coding_genes_coverage")

theme_set(theme_classic())

# remove samples with low coverage
dat = dat[-which(dat[,6]<80),]

g1 <- ggplot(dat, aes(Gene_sequences_coverage, Protein_coding_genes_coverage)) + 
  labs(x = "Genome sequences coverage (%)", y = "Protein-coding genes coverage (%)") + 
  geom_point(shape=21, size=3,alpha=0.4,fill = "dodgerblue3" ,color = "black",stroke = 1) +
  theme(axis.text = element_text( size = rel(1.3 ))) +
  theme(axis.title = element_text( size = rel(1.5 )))

hist_top <- ggplot(dat,aes(Gene_sequences_coverage)) + 
  labs(y = "Count") +
  geom_histogram(colour = 'black',fill = "grey") + 
  theme(axis.text.x = element_blank(),
    axis.text = element_text(size = rel(1.3 )),
    axis.ticks.x = element_blank(),
	axis.title.x = element_blank(),
	axis.title = element_text(size=rel(1.5)))

empty <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(axis.ticks = element_blank(), 
    panel.background = element_blank(), 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),       
    axis.title.x = element_blank(), 
    axis.title.y = element_blank())    

hist_right <- ggplot(dat,aes(Protein_coding_genes_coverage)) + 
  labs(y = "Count",x = "") +
  geom_histogram(colour = 'black',fill="grey") +
  scale_y_continuous(breaks=seq(0, 120, 60),limits = c(0,120)) + 
  theme(axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
	axis.text = element_text(size=rel(1.3 )),
	axis.title = element_text( size=rel(1.5 )),
	plot.margin = unit(c(0, 0.3, 0.15, -0.7), "cm")) +
  coord_flip()

pdf("D:/SYSU/cnSV/Fig/FigureS12A.pdf",width=6,height=6)
grid.arrange(hist_top, empty, g1, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()

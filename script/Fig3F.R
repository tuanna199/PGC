
##### Figure3F scattorplot with sub figure #####

library(ggplot2)
library(grid)

theme_set(theme_classic())

df <- read.csv("D:/SYSU/cnSV/FigDat/Sample_common_SV_genepred_overlap_promoter_all_lof_SV_dist.txt",
				header=TRUE,sep="\t",check.names=FALSE)
				
colnames(df) = c("Name","Length","Count")
df$Length_kb = df$Length / 1000

p1 <- ggplot(df, aes(x = Length_kb,y=Count)) +
  geom_point(fill="darkblue",alpha = 0.4,color="darkblue") +  #color="black",
  labs(x="SV length (kb)",y="Individuals") +
  theme(axis.text = element_text(size=rel(1.4)),axis.title = element_text( size=rel(1.4))) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(plot.margin = unit(c(5.5,12.5,5.5,5.5),"pt"))

p2 <- ggplot(df, aes(x = Length_kb,y=Count)) +
  geom_point(fill="darkorange",alpha = 0.5,color="darkorange") +
  labs(x="SV length (kb)",y="Individuals") +
  xlim(0,20) +
  theme(axis.text = element_text(size=rel(1.1)),axis.title = element_text( size=rel(1.1))) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

pdf("D:/SYSU/cnSV/Fig/Figure3F.pdf",width=4,height=4) 
p1
subfig <- viewport(width = 0.55, height = 0.55, x = 1,y = 0.95,just=c("right","top"))
print(p2,vp=subfig)
dev.off()














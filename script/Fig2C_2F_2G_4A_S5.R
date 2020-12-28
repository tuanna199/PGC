

####### Fig2C
library(ggplot2)

data <- read.csv("D:/SYSU/cnSV/FigDat/fig2C.csv")

data$Type = factor(data$Type,levels=c("DEL","INS","DUP","INV")) 
data$Database = factor(data$Database,levels=c("LRS15","DGV","WGS911","gnomAD"))

pdf("D:/SYSU/cnSV/Figure/Figure2C.pdf",width=6,height=4) 
ggplot(data, aes(fill=Database, y=Number, x=Type)) + 
    geom_bar(position=position_dodge(), stat="identity") +
	scale_fill_manual(values=c("#3CB371", "#1874CD", "#EEC900", "#EE5C42")) +
	geom_text(aes(label=Number), position = position_dodge(width=1),hjust = -0.2) +
	theme_classic() +
	xlab("") +
	ylab("") +
	#theme(legend.position="top",plot.margin = margin(1, 1, 1, -2)) +
	coord_flip()
dev.off()


data2 = data.frame(Database=c("LRS15","DGV","WGS911","gnomAD"),Number=c(97,246,141,36))
data2$Database = factor(data2$Database,levels=c("LRS15","DGV","WGS911","gnomAD"))

pdf("D:/SYSU/cnSV/Figure/Figure2C_2.pdf",width=3,height=2) 
ggplot(data2, aes(x=Database, y=Number, fill=Database)) +
  scale_fill_manual(values=c("#3CB371", "#1874CD", "#EEC900", "#EE5C42"))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Number), position = position_dodge(width=1),hjust = -0.2) +
  theme_classic() +
  xlab("") +
  ylab("") + 
  theme(legend.position="none")+
  coord_flip()
dev.off()

####### Fig2F

library(ggpubr)
library(ggplot2)

dat1 = read.csv("D:/SYSU/cnSV/FigDat/2E_Sample_common_SV_INS_DEL.txt",sep="\t")
dat2 = read.csv("D:/SYSU/cnSV/FigDat/2E_Sample_common_SV_INV_DUP.txt",sep="\t")

dat1$SVType <- factor(dat1$SVType,levels=c("DEL","INS"))
dat2$SVType <- factor(dat2$SVType,levels=c("DUP","INV"))

me_del = median(dat1$SVLength[dat1$SVType=="DEL"])
me_ins = median(dat1$SVLength[dat1$SVType=="INS"])
me_dup = median(dat2$SVLength[dat2$SVType=="DUP"])
me_inv = median(dat2$SVLength[dat2$SVType=="INV"])

p1 <- ggplot(dat1, aes(x = SVLength,fill = SVType)) +
		geom_density(position = "stack", alpha = 0.6) +
		scale_x_continuous(limits=c(50, 3000)) +
		geom_vline(xintercept=me_del) +
		geom_vline(xintercept=me_ins) +
		scale_fill_manual(values=c("#3CB371", "#1874CD"))+
		theme_classic() +
		facet_grid(SVType ~. ,scales="free_x")

p2 <- ggplot(dat2, aes(x = SVLength,fill = SVType)) +
		geom_density(position = "stack", alpha = 0.6) +
		scale_x_continuous(limits=c(50, 9000)) +
		geom_vline(xintercept=me_dup) +
		geom_vline(xintercept=me_inv) +
		scale_fill_manual(values=c("#EEC900", "#EE5C42"))+
		theme_classic() +
		facet_grid(SVType ~. ,scales="free_x")

pdf("D:/SYSU/cnSV/Figure/Fig2E_1125.pdf",width=5,height=4)
ggarrange(p1,p2,ncol=1, nrow=2)
dev.off()


#######   Fig2G  

library(ggplot2)

data = data.frame(
  svtype = c("INV", "DEL", "DUP", "INS"),
  svlength = c( 145.2, 125.7, 104.8, 19.8))

data$fraction <- data$svlength / sum(data$svlength)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2

data$label <- paste0(data$svtype, "\n" ,data$svlength, " Mb")

pdf("D:/SYSU/cnSV/Figure/Fig2G_1228.pdf",width=4,height=4)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=svtype)) +
  geom_rect() +
  geom_label( x=2.5, aes(y=labelPosition, label=label), size=4) +
  scale_fill_manual(values=c("#3CB371", "#1874CD","#EEC900", "#EE5C42"))+
  #scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()


####### Fig4A
library(ggplot2)

data <- read.csv("D:/SYSU/cnSV/FigDat/fig4A.csv")

data$Type = factor(data$Type,levels=c("Singleton","Rare","Low","Common")) 
data$Database = factor(data$Database,levels=c("GWAS","OMIM","COSMIC"))

pdf("D:/SYSU/cnSV/Figure/Figure4A.pdf",width=5,height=4) 
ggplot(data, aes(fill=Type, y=Number, x=Database)) + 
    geom_bar(position="fill", stat="identity", width=0.5) +
	scale_fill_manual(values=c("#41AB6D", "#216EB7", "#E5C31D", "#E15B42")) +
	theme_classic() +
	xlab("") +
	ylab("") +
	theme(legend.position="top") +
	coord_flip()
dev.off()


#######   Fig S5
library(ggplot2)

data <- read.csv("D:/SYSU/cnSV/FigDat/figS5.csv")

data$Type = factor(data$Type,levels=c("DEL","INS","DUP","INV")) 
data$Database = factor(data$Database,levels=c("Singleton","Rare","Low","Common"))

pdf("D:/SYSU/cnSV/Figure/FigureS5.pdf",width=6,height=4) 
ggplot(data, aes(fill=Database, y=Number, x=Type)) + 
    geom_bar(position=position_dodge(), stat="identity") +
	scale_fill_manual(values=c("#3CB371", "#1874CD", "#EEC900", "#EE5C42")) +
	geom_text(aes(label=Number), position = position_dodge(width=1),hjust = -0.2) +
	theme_classic() +
	xlab("") +
	ylab("") +
	#theme(legend.position="top",plot.margin = margin(1, 1, 1, -2)) +
	coord_flip()
dev.off()


data2 = data.frame(Database=c("Singleton","Rare","Low","Common"),Number=c(377,167,101,151))
data2$Database = factor(data2$Database,levels=c("Singleton","Rare","Low","Common"))

pdf("D:/SYSU/cnSV/Figure/FigureS5_2.pdf",width=3,height=2) 
ggplot(data2, aes(x=Database, y=Number, fill=Database)) +
  scale_fill_manual(values=c("#3CB371", "#1874CD", "#EEC900", "#EE5C42"))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Number), position = position_dodge(width=1),hjust = -0.2) +
  theme_classic() +
  xlab("") +
  ylab("") + 
  theme(legend.position="none")+
  coord_flip()
dev.off()


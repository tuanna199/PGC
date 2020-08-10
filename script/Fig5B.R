
#####  Figure5B-1-Fst #####
Fst = read.csv("D:/SYSU/cnSV/FigDat/chr14_region_fst.txt",sep="\t",header=FALSE)
Fst = na.omit(Fst)
colnames(Fst) = c("Chr_start","start","Chr_end","end","SVlen","Type","Fst")

p1 <- ggplot(Fst, aes(x=start)) + 
  geom_point(aes(y=Fst),color='darkblue') +
  xlim(c(105500000,105800000)) +
  labs(x="") + 
  theme_bw() +
  theme(plot.margin = unit(c(50,5.5,0,5.5),"pt"),
		axis.ticks.x = element_blank(),
		axis.title.x = element_blank())

#####  Figure5B-2-Structure #####
library(ggplot2)
library(gggenes)

gene_list = read.csv("D:/SYSU/cnSV/FigDat/chr14_region_gene.txt",sep="\t",header=F)
colnames(gene_list) = c("chr","start","end","gene")
gene_list$length = gene_list$end - gene_list$start
gene_list = aggregate(. ~ gene ,data=gene_list,max)
gene_list = gene_list[order(gene_list$start),]
gene_list$Chr = "Chr14"
gene_list$gene = factor(gene_list$gene,order=TRUE,levels=gene_list$gene)
ref_loc = read.csv("D:/SYSU/cnSV/FigDat/Chr14_genes.txt",sep="\t",header=F)[,-1]
colnames(ref_loc) = c("gene","anno_start","anno_end")

dat = cbind(gene_list,ref_loc[,2:3])
dat = na.omit(dat)
p2 <- ggplot(dat, aes(xmin = anno_start, xmax = anno_end, y = "0.050",fill = gene,label = gene)) +
      geom_gene_arrow(arrowhead_height = unit(5, "mm"), 
                      arrow_body_height = unit(6, "mm"),
    				  arrowhead_width = unit(4, "mm")) +
      xlim(c(105500000,105800000)) +
      labs(y = "Chr14") + 
      theme_genes() +
      theme(legend.position = "bottom") +
      theme(plot.margin = unit(c(80,5.5,0,5.5),"pt")) +
      theme(legend.text=element_text(face="italic"))


#####  Figure5B-3-Repeats #####
Rep = read.csv("D:/SYSU/cnSV/FigDat/ref_chr14_105500000_105800000.txt",sep="\t",header=T)

p3 <- ggplot(Rep,
      aes(x=genoStart, xend=genoEnd,y="0.050",yend="0.050")) + 
	  geom_segment(size=5) +
	  labs(x="",y="Tandem repeat") + theme_bw() +
	  theme(plot.margin = unit(c(0,5.5,0,5.5),"pt"),
	  axis.ticks.x = element_blank(),
	  axis.title.x = element_blank())


#####  Figure5B-4-SV #####
library(ggplot2)
data = read.csv("D:/SYSU/cnSV/FigDat/chr14_region_genotype_select.csv")
data = na.omit(data)
data$SV = paste0("SV",1:7)

p4 <- ggplot(data,
     aes(x=Pos1, xend=Pos2,y="0.050",yend="0.050",color=SV,fill=SV)) + 
   	 geom_segment(size=10) +
   	 scale_colour_brewer(palette = "Dark2") +  
   	 scale_fill_brewer(palette = "Dark2") +
   	 xlim(c(105500000,105800000)) +
   	 labs(x="",y="SV name") + 
	 theme_bw() +
   	 theme(legend.position = "bottom") +
   	 theme(plot.margin = unit(c(0,5.5,0,5.5),"pt"))


#####  Figure5B-5-heatmap #####
library(ComplexHeatmap)
library(RColorBrewer)

data = read.csv("D:/SYSU/cnSV/FigDat/chr14_region_genotype_select.csv")
Area = t(data[1,7:ncol(data)])
colnames(Area) = "Area"
data = data[-1,]

dat = data[,7:ncol(data)]
rownames(dat) = paste0(data$Type,":",data$Pos1,"-",data$Pos2)

col2 <- setNames(RColorBrewer::brewer.pal(name = "Dark2", n = 7),rownames(dat))
ha = rowAnnotation(Area = Area, 
                   col = list(Area = c("South" = "blueviolet","North" = "brown")),
                   show_annotation_name = FALSE)

column_ha = HeatmapAnnotation(SV = rownames(dat),col=list(SV = col2),show_annotation_name=FALSE)

Dat = t(dat)
Dat = apply(Dat,2,function(x) as.numeric(x))

colors = list("0" = "snow2","1" = "steelblue", "2" = "brown1")
pdf("D:/SYSU/cnSV/Fig/chr14_down.pdf",width=10,height=8)
Heatmap(Dat,name = "Genotype",
		width = unit(11, "cm"), height = unit(9, "cm"),
		show_row_names = FALSE, 
		cluster_columns=FALSE,
		cluster_rows=FALSE,
        show_column_names = FALSE,
		top_annotation = column_ha,
		left_annotation = ha,
		heatmap_height = unit(0.1, "npc"),		 
		col=colors) 
dev.off()


#####    merge all figures   #####
library(gridExtra)
pdf("D:/SYSU/cnSV/Fig/Figure5_chr14_upper.pdf")
grid.arrange(p1,p2, p3, p4,ncol=1, nrow=4,heights=c(0.5,0.3,0.1,0.3))
dev.off()











################  Gene structure  ################
bed = read.table("D:/SYSU/cnSV/FigDat/hg38_ensGene_gtf95_sort.bed")

colnames(bed) = c("chromosome","start","end","gene")

### 4c_1 chr16 165396-184704
chr16_bed = subset(bed,chromosome=="16" & start>=165396 & end <= 184704)

### 4c_2 chr11 5201674-5229059
chr11_bed = subset(bed,chromosome=="11" & start>=5201674 & end <= 5229059)

### 4d 5_42628460-5_42630888-2433-DEL
bed2 = read.table("D:/SYSU/cnSV/FigDat/hg38_ensGene_dechr.txt")
chr5_bed = subset(bed2,V2 == "ENST00000230882.8")
start = as.numeric(unlist(strsplit(chr5_bed$V10,",")))
end = as.numeric(unlist(strsplit(chr5_bed$V11,",")))
fig4d_bed = data.frame(gene="GHR",exon=paste0("exon",1:10),start,end)
fig4d_bed = rbind(fig4d_bed,c("GHR","SV",42628460,42630888))
fig4d_bed$start = as.numeric(fig4d_bed$start)
fig4d_bed$end = as.numeric(fig4d_bed$end)


################  Figure 4C&4d gene structure  ################
library(ggplot2)
library(gggenes)

p1 <- ggplot(chr11_bed, aes(xmin = start, xmax = end, y = chromosome, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(3, "mm")) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

p2 <- ggplot(chr16_bed, aes(xmin = start, xmax = end, y = chromosome, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(3, "mm")) +
  #geom_gene_label(align = "left") +
  #geom_blank(data = dummies) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

pdf("D:/SYSU/cnSV/Fig/fig4d_gene.pdf")
ggplot(fig4d_bed, aes(xmin = start, xmax = end, y = gene, fill = exon, label = exon)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(3, "mm")) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()
dev.off()

library(gridExtra)
pdf("D:/SYSU/cnSV/Fig/fig4c_gene.pdf")
grid.arrange(p1,p2,ncol=1, nrow=2)
dev.off()


#################  Figure S9 gene structure  #################
library(ggplot2)
library(gggenes)

bed <- read.csv("D:/SYSU/cnSV/FigDat/structure_716.csv",header=FALSE)

gggene_input <- function(bed,gene,SV_start,SV_end){
  bed_gene = subset(bed,V13==as.character(gene))
  start = as.numeric(unlist(strsplit(bed_gene$V10,",")))
  end = as.numeric(unlist(strsplit(bed_gene$V11,",")))
  bed1_dat = data.frame(gene=as.character(gene),exon=paste0("exon",1:length(start)),start,end)
  bed1_dat = rbind(bed1_dat,c(as.character(gene),"SV",SV_start,SV_end))
  bed1_dat$start = as.numeric(bed1_dat$start)
  bed1_dat$end = as.numeric(bed1_dat$end)
  return(bed1_dat)
}

### MC1R 16_89916622-16_89919621-2999-DEL	
MC1R = gggene_input(bed,"MC1R",89916622,89919621)

### BAX	19_48935540-19_48959076-23536-DUP
BAX = gggene_input(bed,"BAX",48935540,48959076)


p1 <- ggplot(MC1R, aes(xmin = start, xmax = end, y = gene, fill = exon, label = exon)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(3, "mm")) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()


p2 <- ggplot(BAX, aes(xmin = start, xmax = end, y = gene, fill = exon, label = exon)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(3, "mm")) +
  theme_genes()

library(gridExtra)
pdf("D:/SYSU/cnSV/Fig/Fig_s9_upper.pdf",width=10,height=10)
grid.arrange(p1,p2,ncol=1,nrow=2)
dev.off() 

ggsave("D:/SYSU/cnSV/Fig/structures_2_MC1R.pdf",plot=p2)
ggsave("D:/SYSU/cnSV/Fig/structures_4_BAX.pdf",plot=p4,width=15,height=10)


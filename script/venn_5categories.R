#!/usr/bin/env R

args <- commandArgs(TRUE)

library(grid)
library(futile.logger)
library(VennDiagram)



venn_plot <- function(data_file, out_name, out_type, out_height, out_width, text_cex){
    A = levels(data_file$LRS15)
    B = levels(data_file$DGV)
    C = levels(data_file$dbVar)
    D = levels(data_file$WGS17795)
    E = levels(data_file$gnomAD)
    print(length(A))

    # venn.diagram(
    #   x = list(A = A, B = B, C = C, D = D),
    #   category.names = c("Cell2019","DGV", "dbVar","gnomAD"),
    #   filename = out_name,
    #   imagetype = out_type, #"svg", "png"
    #   height = out_height, #15, 3000
    #   width = out_width, #15, 3000
    #   col = "transparent",
    #   fill = c("cornflowerblue", "seagreen3", "orange", "darkorchid1"),
    #   alpha = 0.50,
    #   label.col = c("white", "white", "white", "white", 
    #       "white", "white", "white", "white", "white", "white", 
    #       "white", "white", "white", "white", "white"),
    #   cex = text_cex, #2.5, 1
    #   fontfamily = "serif",
    #   fontface = "bold",
    #   cat.col = c("cornflowerblue", "darkgreen", "orange", "darkorchid1"),
    #   cat.pos = c(-20,20,0,0),
    #   cat.dist = c(0.21,0.21,0.10,0.09),
    #   cat.fontfamily = "serif", 
    #   cat.cex = text_cex, #2.5, 1
    #   rotation.degree = 0,
    #   margin = 0.2 )

    labelColors <- rep("black", 31)
    catColors <- rep("black", 5)

    venn.diagram(
        x = list(A = A, B = B, C = C, D = D, E = E),
        category.names = c("dbVar", "DGV", "LRS15",  "WGS17795", "gnomAD" ),
        filename = out_name,
        imagetype = out_type, #"svg", "png"
        height = out_height, #15, 3000
        width = out_width, #15, 3000
        col = "transparent",
        # fill = c("cornflowerblue", "seagreen3", "orange", "darkorchid1"),
        fill = c( "aquamarine3", "deepskyblue3", "royalblue3", "purple3", "orangered3"),
        alpha = 0.30,
        label.col = labelColors,
        cex = text_cex, #2.5, 1
        fontfamily = "serif",
        # fontface = "bold",
        cat.col = catColors,
        # cat.col = c( "aquamarine1", "deepskyblue", "royalblue2", "purple2", "orangered3"),
        cat.pos = c(-72, 0, 72, 144, -136),
        # cat.dist = c(0.26, 0.21, -0.23, 0.23, -0.23),
        cat.dist = c(-0.11, -0.05, -0.09, 0.06, -0.08),
        cat.fontfamily = "serif", 
        cat.cex = text_cex * 1.1, #2.5, 1
        cat.fontface = "bold",
        rotation.degree = 0,
        margin = 0.2 )

}



data_file <- read.table(args[1],header=TRUE,sep="\t",quote="")
print(head(data_file))
#venn_plot(data_file, "noncoding.svg", "svg", 15, 15, 2.5)
#venn_plot(data_file, "noncoding.png", "png", 3000, 3000, 1)

venn_plot(data_file, args[2], args[3], as.numeric(args[4]), as.numeric(args[5]), as.numeric(args[6]))

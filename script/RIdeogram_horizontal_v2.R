# The global coordinates of the chromosomes
karyotype <- read.table("GCF_000001405.39_GRCh38.p13_genomic.karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
karyotype <- karyotype[, 1:5]
# A heat map of the interior of a chromosome
heterochromatin <- read.table("intersect_two_end.forR_heatmap.txt", sep = "\t", header = T, stringsAsFactors = F)
# Some annotation on the outside of a chromosome
two_end <- read.table("intersect_two_end.forR_feature.txt", sep = "\t", header = T, stringsAsFactors = F)
# The final output is to an SVG file
output <- "RIdeogram_horizontal_v2_test1.svg"


mpx<-3.543307
width <- 170
chr_width <- width / (2.6*nrow(karyotype)) * mpx

karyotype$y9 <- karyotype$y1 <- karyotype$y8 <- karyotype$y4 <- karyotype$y5 <-
      karyotype$y12 <- apply(data.frame(1:nrow(karyotype)),1,function(x)(20*mpx+(x[1]-1)*2.6*chr_width))

maxchrlen<-200

karyotype$x1 <- karyotype$x2 <-
      apply(data.frame(karyotype$End),1,
            function(x) ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 
              ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

karyotype$y10 <- karyotype$y2 <- karyotype$y3 <- karyotype$y7 <-
      karyotype$y6 <- karyotype$y11 <- karyotype$y1+chr_width

karyotype$x3 <- karyotype$x8 <-
        apply(data.frame(karyotype$End,karyotype$CE_start),1,
              function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx - 
                ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

karyotype$x4 <- karyotype$x7 <-
        apply(data.frame(karyotype$End,karyotype$CE_end),1,
              function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx - 
                ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

karyotype$x5 <- karyotype$x6 <- 
        apply(data.frame(karyotype$End),1,
              function(x)((25+maxchrlen)*mpx-chr_width/2 - 
                ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

karyotype$x9 <- karyotype$x10 <-
      apply(data.frame(karyotype$End),1,
            function(x)((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx - 
              ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

karyotype$x11 <- karyotype$x12 <- 
      apply(data.frame(karyotype$End),1,
            function(x)((25+maxchrlen)*mpx - 
              ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

karyotype$path = paste("<path d=\"M", karyotype$x1, ",", karyotype$y1,
                             " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x2, ",", karyotype$y2,
                             " L", karyotype$x3, ",", karyotype$y3,
                             " L", karyotype$x4, ",", karyotype$y4,
                             " L", karyotype$x5, ",", karyotype$y5,
                             " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x6, ",", karyotype$y6,
                             " L", karyotype$x7, ",", karyotype$y7,
                             " L", karyotype$x8, ",", karyotype$y8,
                             " Z" ,"\" style=\"fill:none; stroke:grey; stroke-width:1\"/>", sep = "")

karyotype$hat = paste("<path d=\"M", karyotype$x1, ",", karyotype$y1,
                          " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x2, ",", karyotype$y2,
                          " L", karyotype$x10, ",", karyotype$y10,
                          " L", karyotype$x9, ",", karyotype$y9,
                          " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

karyotype$shoe = paste("<path d=\"M", karyotype$x5, ",", karyotype$y5,
                           " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x6, ",", karyotype$y6,
                           " L", karyotype$x11, ",", karyotype$y11,
                           " L", karyotype$x12, ",", karyotype$y12,
                           " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

karyotype$bow = paste("<path d=\"M", karyotype$x8, ",", karyotype$y8,
                          " L", karyotype$x7, ",", karyotype$y7,
                          " L", karyotype$x3, ",", karyotype$y3,
                          " L", karyotype$x4, ",", karyotype$y4,
                          " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

karyotype$text = paste("<text x=\"", 615.2521 - (150 + 25) * mpx - 15 + 82, "\" y=\"",
                           (karyotype$y1 + karyotype$y2)/2 + nchar(karyotype$Chr) * 2.2 + 2,
                           "\" font-size=\"13\" font-family=\"Arial\" fill=\"black\" >",
                           karyotype$Chr, "</text>", sep = "")


# A heat map of the interior of a chromosome
# cnum<-10000
mydata <- heterochromatin[, 1:6]
# colorset1 <- c("#4575b4", "#ffffbf", "#d73027")
# mydata$color <- colorRampPalette(colorset1)(cnum)[round(scales::rescale(mydata$Value,to=c(1,cnum)))]

mydata<-merge(mydata,data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x1=karyotype$y1,x2=karyotype$y2),by="Chr")

mydata$y1 <- apply(data.frame(mydata$ChrEnd,mydata$Start),1,
  function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx - 
    ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

mydata$y2 <- apply(data.frame(mydata$ChrEnd,mydata$End),1,
  function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx - 
    ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

mydata$rect = paste("<path d=\"M", mydata$y1, ",", mydata$x1,
                          " L", mydata$y2, ",", mydata$x1,
                          " L", mydata$y2, ",", mydata$x2,
                          " L", mydata$y1, ",", mydata$x2,
                          " Z" ,"\" style=\"fill:#", mydata$color, "; stroke:#", 
                          mydata$color, "; stroke-width:0.25\"/>", sep = "")


Lx <- 160
Ly <- 70

# mydata_legend <- data.frame(color=colorRampPalette(colorset1)(cnum))
# mydata_legend$x1<-apply(data.frame(1:cnum),1,function(x)(Lx * mpx + (x[1] - 1) * ((20 * mpx) / cnum)))
# mydata_legend$y1<-Ly * mpx
# mydata_legend$legend<-paste("<rect x=\"", mydata_legend$x1, "\" y=\"", mydata_legend$y1,
#                           "\" width=\"", (20 * mpx) / cnum,
#                           "\" height=\"", 4 * mpx,
                          # "\" style=\"fill:", mydata_legend$color, ";stroke:none\"/>", sep = "")

# legend_text <- c(paste("<text x=\"", min(mydata_legend$x1),
#                      "\" y=\"", min(mydata_legend$y1) + 8 * mpx - 3,
#                      "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >Low</text>", sep = ""),
#                paste("<text x=\"", max(mydata_legend$x1) - 20,
#                      "\" y=\"", min(mydata_legend$y1) + 8 * mpx - 3,
#                      "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >High</text>", sep = "")
# )

mydata_legend <- mydata[!duplicated(mydata$Type), 1:6]
mydata_legend <- mydata_legend[order(mydata_legend$color),]
mydata_legend$x <- Lx * mpx
mydata_legend$y<-apply(data.frame(1:nrow(mydata_legend)),1,
                                function(x)(Ly * mpx + 4 + (12 + (x[1] - 1) * 4 ) * mpx - 44))

mydata_legend$x1<-mydata_legend$x + 4
mydata_legend$y1<-mydata_legend$y- 4 * mpx / 2

mydata_legend$legend <- paste("<rect x=\"",
                              mydata_legend$x1 - 20, "\" y=\"",
                              mydata_legend$y1 - 2, "\" width=\"", 24,
                              "\" height=\"", 4, "\" style=\"fill:#",
                              mydata_legend$color, ";stroke:none\"/>", sep = "")

mydata_legend$legend_text <- paste("<text x=\"",
                              mydata_legend$x + 15, "\" y=\"",
                              mydata_legend$y - (4 * mpx / 2 - 4),
                              "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >",
                              mydata_legend$Type,"</text>", sep = "")



# Some annotation on the outside of a chromosome
mydata_interval <- two_end[, 1:6]
mydata_interval<-merge(mydata_interval,
  data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x2=karyotype$y2),by="Chr")

mydata_interval$x <- mydata_interval$x2 + chr_width / 2
mydata_interval$y0 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,
  function(x)((25+maxchrlen*(1-(x[1]-(x[2]+x[3])/2)/max(karyotype$End))) * mpx - 
    ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2 - 93.4093)))

mydata_interval<-mydata_interval[order(mydata_interval$Chr,mydata_interval$y0),]

repel<- function(mydata,myforce,tag){
          if(min(diff(mydata)) >= myforce | tag==3000){return(mydata)}
          tag<-tag+1
          sp<-which.min(diff(mydata))
          ep<-sp+1
          mydata[sp]<-mydata[sp]-myforce
          mydata[ep]<-mydata[ep]+myforce
          mydepo<-sort(mydata)
          return(repel(mydepo,myforce,tag))
        }

mydata_interval$y <- NA
for (chr in unique(mydata_interval$Chr)){
  if(nrow(mydata_interval[mydata_interval$Chr==chr,])>1){
    mydata_interval[mydata_interval$Chr == chr,]$y <-
    repel(mydata_interval[mydata_interval$Chr == chr,]$y0, chr_width/3, 1)
  }else{mydata_interval[mydata_interval$Chr == chr,]$y <- mydata_interval[mydata_interval$Chr == chr,]$y0}
}

mydata_interval_triangle<- mydata_interval[mydata_interval$Shape == "triangle",]
if ("triangle" %in% mydata_interval$Shape){
  mydata_interval_triangle$interval <- 
  paste("<path d=\"M", mydata_interval_triangle$y + chr_width / 4, ",", mydata_interval_triangle$x - chr_width / 4,
        " L", mydata_interval_triangle$y + chr_width / 4, ",", mydata_interval_triangle$x + chr_width / 4,
        " L", mydata_interval_triangle$y - chr_width / 4, ",", mydata_interval_triangle$x,
        " Z" ,"\" style=\"fill:#", mydata_interval_triangle$color, ";stroke:none\"/>", sep = "")
}

mydata_interval_box<- mydata_interval[mydata_interval$Shape == "box",]
if ("box" %in% mydata_interval$Shape){
  mydata_interval_box$interval <- 
  paste("<rect x=\"", mydata_interval_box$y - chr_width / 4,
        "\" y=\"", mydata_interval_box$x - chr_width / 4,
        "\" width=\"", chr_width / 2,
        "\" height=\"", chr_width / 2,
        "\" style=\"fill:#", mydata_interval_box$color, "; stroke:none\"/>", sep = "")
}

mydata_interval_circle<- mydata_interval[mydata_interval$Shape == "circle",]
if ("circle" %in% mydata_interval$Shape) {
  mydata_interval_circle$interval <- 
  paste("<circle cx=\"", mydata_interval_circle$y,
       "\" cy=\"", mydata_interval_circle$x,
       "\" r=\"", chr_width / 4,
       "\" style=\"fill:#", mydata_interval_circle$color, "; stroke:none\"/>", sep = "")
}

mydata_interval<-rbind(mydata_interval_box,mydata_interval_circle,mydata_interval_triangle)
mydata_interval$line <- paste("<line x1=\"", mydata_interval$y0,
                              "\" y1=\"", mydata_interval$x2,
                              "\" x2=\"", mydata_interval$y,
                              "\" y2=\"", mydata_interval$x,
                              "\" style=\"stroke:#", mydata_interval$color, "; stroke-width:0.25\"/>", sep = "")


mydata2_legend <- mydata_interval[!duplicated(mydata_interval$Type), 1:6]
mydata2_legend <- mydata2_legend[order(mydata2_legend$Shape, mydata2_legend$color),]
mydata2_legend$x <- Lx * mpx - 8
mydata2_legend$y<-apply(data.frame(1:nrow(mydata2_legend)),1,
                                function(x)(Ly * mpx + 4 + (12 + (x[1] - 1) * 4 ) * mpx))

mydata2_legend$x1<-mydata2_legend$x + 4
mydata2_legend$y1<-mydata2_legend$y- 4 * mpx / 2

for (i in 1:nrow(mydata2_legend)){
  if (mydata2_legend[i, 3] == "triangle") {
    mydata2_legend[i,11] = paste("<path d=\"M", mydata2_legend[i,9] - 4, ",",
                                 mydata2_legend[i,10] + 4, " L", mydata2_legend[i,9] + 4, ",",
                                 mydata2_legend[i,10] + 4, " L", mydata2_legend[i,9], ",",
                                 mydata2_legend[i,10] - 4, " Z" ,"\" style=\"fill:#",
                                 mydata2_legend[i, 6], "; stroke:none\"/>", sep = "")
  } else if (mydata2_legend[i, 3] == "box") {
    mydata2_legend[i,11] <- paste("<rect x=\"",
                                  mydata2_legend[i,9] - 4, "\" y=\"",
                                  mydata2_legend[i,10] - 4, "\" width=\"", 8,
                                  "\" height=\"", 8, "\" style=\"fill:#",
                                  mydata2_legend[i, 6], ";stroke:none\"/>", sep = "")
  } else if (mydata2_legend[i, 3] == "circle") {
    mydata2_legend[i,11] <- paste("<circle cx=\"",
                                  mydata2_legend[i,9], "\" cy=\"",
                                  mydata2_legend[i,10], "\" r=\"", 4, "\" style=\"fill:#",
                                  mydata2_legend[i, 6], ";stroke:none\"/>", sep = "")
  }
}

names(mydata2_legend)[11] <- "shape"

for (i in 1:nrow(mydata2_legend)){
  mydata2_legend[i,12] <- paste("<text x=\"",
                                mydata2_legend[i,7] + 15 + 8, "\" y=\"",
                                mydata2_legend[i,8] - (4 * mpx / 2 - 4),
                                "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >",
                                mydata2_legend[i,2],"</text>", sep = "")
}
names(mydata2_legend)[12] <- "name"



# The final output is to an SVG file
first_line <- c("<?xml version=\"1.0\" standalone=\"no\"?>",
                  "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">",
                  "",
                  paste("<svg id=\"svg\" width=\"1116.14175\" height=\"1052.362\">", "\t")
  )


cat(first_line, file = output)


# heat map
cat(mydata$rect, file = output, append = TRUE)
cat(mydata_legend$legend_text, file = output, append = TRUE)
cat(mydata_legend$legend, file = output, append = TRUE)


# annotated map
cat(mydata_interval$interval, file = output, append = TRUE)
cat(mydata_interval$rect, file = output, append = TRUE)
cat(mydata_interval$line, file = output, append = TRUE)
cat(mydata2_legend$shape, file = output, append = TRUE)
cat(mydata2_legend$name, file = output, append = TRUE)


# chromosome map
cat(karyotype$hat, file = output, append = TRUE)
cat(karyotype$shoe, file = output, append = TRUE)
cat(karyotype$bow, file = output, append = TRUE)
cat(karyotype$path, file = output, append = TRUE)
cat(karyotype$text, file = output, append = TRUE)


cat("</svg>", file = output, append = TRUE)

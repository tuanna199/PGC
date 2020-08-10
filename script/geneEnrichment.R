#!/usr/bin/Rscript
library(argparser)
library(ggplot2)
library(readr)


#usage: Rscript ~/github/NanoHub/script/geneEnrichment.R --input /home/wuzhikun/Project/NanoTrio/AnnotSV/Enrichment/KEGG_2019_Human..enrichr.reports.txt --pdf test.pdf --width 8 --height 4 --selectNum 20 --requestNum 355

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output pdf file.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
arg <- add_argument(arg, "--requestNum", default=0, help="Request gene number.")
arg <- add_argument(arg, '--background', default=19041, help='The background gene number.')
arg <- add_argument(arg, '--selectNum', help='The selected top number.')
argv <- parse_args(arg)



gene_set_enrichment <- function(file, pdf_file, width, height, selectNum, requestNum, background){
    backgroundLength <- as.numeric(background)

    width <- as.numeric(width)
    height <- as.numeric(height)
    selectNum <- as.numeric(selectNum)


    if (!file.exists(file)) {
      print(paste0(file, " did not exists"))
      quit()
    }

    if (file.info(file)$size == 0){
            system(paste('touch', args[2], sep=' '))
            quit()
    }

    raw_data <- read_tsv(file)


    if (!"Odds Ratio" %in% colnames(raw_data)) {

        questLength <- as.numeric(requestNum)

        data <- subset(raw_data, select=c("Term", "P-value", "Genes"))
        colnames(data) <- c("Term", "P_Value", "Genes")

        Overlap2 <- strsplit(raw_data$Overlap, "/")
        gene_numbers <-  matrix(as.numeric(unlist(Overlap2)), ncol=2, byrow=TRUE)
        data$Rich_factor <- gene_numbers[,1] * backgroundLength / (gene_numbers[,2] * questLength)
        # print(data$Rich_factor)


    } else{

        data <- subset(raw_data, select=c("Term", "P-value", "Odds Ratio", "Genes"))
        colnames(data) <- c("Term", "P_Value", "Rich_factor", "Genes")

    }


    Number <- c()
    for (i in strsplit(data$Genes, ";")){
      Number <- c(Number, length(i))
    }

    data$Input_number <- Number
      

    #if (!"Rich_factor" %in% colnames(data)) {
    #        data$Rich_factor <- data$Input_number/data$Background_number
    #}

    data <- data[order(data$P_Value),]

    if(nrow(data) < selectNum){data <- data}else{data <- data[1:selectNum,]}
    
    if(max(nchar(as.vector(data$Term)))>90){t_size=7}else{t_size=9}


    p <- ggplot(data, aes(-log10(P_Value), reorder(Term, -log10(P_Value)) )) +
        theme_bw() +
        geom_point(aes(size=Input_number, color=Rich_factor)) +
        #xlim(150, 280) +
        xlab(bquote("-" ~log[10]~ "(P value)")) +
        ylab("Terms")+
        scale_color_gradient(low = "blue",high="red")+
        theme(axis.text=element_text(size=t_size))+
        scale_y_discrete(labels = function(y) lapply(strwrap(y, width = 55, simplify = FALSE), paste, collapse="\n")) +
        labs(colour = "Odds ratio", size="Gene number") +
        theme(legend.position = "right") + # bottom
        theme(legend.text=element_text(size=rel(0.5))) +
        theme(legend.title = element_text(size = 7)) +
        theme(legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3, "cm"), legend.key.size=unit(0.3, "cm"))
    p

    ggsave(pdf_file, width=width, height=height)

}


gene_set_enrichment(argv$input, argv$pdf, argv$width, argv$height, argv$selectNum, argv$requestNum, argv$background)





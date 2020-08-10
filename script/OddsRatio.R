#!/usr/bin/Rscript
library(readr)
library(argparser)

arg <- arg_parser('Calculate the odds ratio and pvalue.')
arg <- add_argument(arg, '--input', help='The file with annotation statistics.')
arg <- add_argument(arg, '--odds', help='The output file with odds ratio.')
arg <- add_argument(arg, '--pvalue', help='The output file with pvalue.')
argv <- parse_args(arg)


#usage: Rscript ~/github/NanoHub/script/OddsRatio.R --input /home/wuzhikun/Project/Population/population/Annotation/Sample_common_SV_annotation_stats.xls --odds temp_odds.xls --pvalue temp_pavlue.xls


annotation_odds_ratio <- function(in_file, odds_file, pvalue_file){
    #in_file :
    # Category    CDS DownStream  Exon    Intron  UTR3    UTR5    UpStream    Intergenic  All
    # All 3147    4686    2653    63852   1581    1628    4804    63767   132356
    # Common  620 979 461 15478   243 270 968 16222   32782
    # Low 290 465 242 6581    154 143 496 7248    14362
    # Rare    597 999 542 14182   300 319 1038    13689   28951
    # Singleton   1640    2243    1408    27611   884 896 2302    26608   56261

    data <- read_tsv(in_file)

    rowNum <- nrow(data)
    colNum <- ncol(data)

    rnames <- data[2:nrow(data), 1]

    cnames <- colnames(data[,1:(ncol(data)-1)])


    Estimate <- c()
    Pvalue <- c()

    for (i in 2:5){
      for (j in 2:6){
        target <- c(data[i, j], data[i, colNum], data[1, j], data[1, colNum])
        target <- as.numeric(target)
        dataM <- matrix(data=target, nrow=2, byrow=TRUE)
        test <- fisher.test(dataM)
        Estimate <- c(Estimate, test$estimate)
        Pvalue <- c(Pvalue, test$p.value)
      }
    }


    EstimateMatrix <- matrix(data=Estimate, nrow=4, byrow=TRUE)
    PvalueMatrix <- matrix(data=Pvalue, nrow=4, byrow=TRUE)

    EstimateData <- as.data.frame(EstimateMatrix)
    EstimateData <- cbind(as.vector(unlist(rnames)), EstimateData)
    colnames(EstimateData) <- cnames


    PvalueData <- as.data.frame(PvalueMatrix)
    PvalueData <- cbind(as.vector(unlist(rnames)), PvalueData)
    colnames(PvalueData) <- cnames

    write.table(EstimateData, odds_file, sep="\t", quote=FALSE, row.names = FALSE)
    write.table(PvalueData, pvalue_file, sep="\t", quote=FALSE, row.names = FALSE)
    
}


annotation_odds_ratio(argv$input, argv$odds, argv$pvalue)


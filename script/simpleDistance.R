#!/usr/bin/env Rscript
library(argparser)
library(scrime)


arg <- arg_parser('Calculate simple matching coefficients.')
arg <- add_argument(arg, '--input', help='Input file.')
arg <- add_argument(arg, '--out', help='output file with distances of sample pairs.')
argv <- parse_args(arg)


simple_matching_coefficients <- function(in_file, out_file){

    mat <- read.table(in_file, header=FALSE, row.names=1)
    mat <- as.matrix(mat)
    distance <- smc(mat, dist=TRUE)
    dist_frame <-  as.data.frame(distance)
    distance <- round(distance, 4)
    write.table(dist_frame, out_file, sep="\t", quote=FALSE)

}



simple_matching_coefficients(argv$input, argv$out)

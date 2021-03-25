#!/usr/bin/Rscript
library.path <- .libPaths()
library("tm", lib.loc = library.path)
library("SnowballC", lib.loc = library.path)
library("RColorBrewer", lib.loc = library.path)
library("XML", lib.loc = library.path)
library("argparse", lib.loc = library.path)


args <- commandArgs(trailingOnly = TRUE)

# Reading in output CSV file
words = read.csv(args[1], header = FALSE)

#outfile = paste(args[2], args[3], sep = "/", collapse = NULL)
tiff(args[2], units="in", width=20, height=10, res=100)
wordcloud(words$V1, scale = c(3,0.25), fixed.asp = TRUE)
dev.off()

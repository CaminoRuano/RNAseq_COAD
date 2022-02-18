


## RNAseq Analysis on Rstudio

#library import 

library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(DESeq2)


setwd("C://Users/Camino/Documents/EpigeneLabs/")

## Load sample info and RAW counts


sampleinfo <- read.csv("./data/metadata/gdc_sample_sheet.2022-02-18.tsv",sep = "\t", header = TRUE)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleinfo, 
                                       directory = "./data/data",
                                       design= ~ 1)

## Load clinical data

clinical_data = read.csv("./data/metadata/clinical.tsv", sep="\t")



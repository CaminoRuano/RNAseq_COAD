

# RNAseq Analysis on Rstudio



## RNAseq Analysis on Rstudio

setwd("C://Users/Camino/Documents/EpigeneLabs/")

#library import 

library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("DESeq2")
library("TCGAutils")
library('GenomicDataCommons')
library('dplyr')

### Download the data ###

GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

TCGAbiolinks:::getProjectSummary("TCGA-COAD")


  ### No metastasis data
query_TCGA = GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts", sample_type=("Primary Tumor", "Solid Tissue Normal"))


lihc_res = getResults(query_TCGA)
colnames(lihc_res)

GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)
dim(tcga_data)

colnames(colData(tcga_data))
table(tcga_data@colData$tissue_or_organ_of_origin)

head(rowData(tcga_data))

saveRDS(object = tcga_data,
        file = "./tcga_data.RDS",
        compress = FALSE)


## Metadata file = match file names with case name from clinical data txt


# file names

sample_file = read.csv("./metadata/gdc_sample_sheet.2022-02-19.tsv", sep="\t", header=TRUE)
sample_file = sample_file[sample_file$Sample.Type=="Primary Tumor",]
uniq = unique(sample_file$Case.ID)   # check the n° of unique cases - sometimes there are samples with several sequences or samples with tumour and solid normal samples-

sample_file <- as.data.frame(sample_file %>% dplyr::group_by(Case.ID) %>% dplyr::filter (! duplicated(Case.ID)))

# clinical data

clinical_data = read.csv("./metadata/clinical.tsv", sep="\t", header=TRUE)
clinical_data = clinical_data[clinical_data$treatment_type=="Radiation Therapy, NOS",]
clinical_data = dplyr::filter(clinical_data, clinical_data$case_submitter_id %in% sample_file$Case.ID)

clinical_data <- as.data.frame(clinical_data %>% dplyr::group_by(case_submitter_id) %>% dplyr::filter (! duplicated(case_submitter_id)))


sample_file = dplyr::filter(sample_file, sample_file$Case.ID %in% clinical_data$case_submitter_id)

# merge both files together
rownames(sample_file) = sample_file$Case.ID
rownames(clinical_data)=clinical_data$case_submitter_id


metadata = cbind(sample_file, clinical_data)



### Create DESeq file for normalization

data = DESeqDataSetFromHTSeqCount(sampleTable = metadata, directory = "./GDCdata", design ~ Sample.Type)



#### matching Case.UUIDs to File.UUIDs

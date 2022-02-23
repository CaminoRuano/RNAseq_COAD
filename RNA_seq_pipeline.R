


## RNAseq Analysis on Rstudio



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
library('janitor')
library('stringr')

### Download the data colon cancer samples###

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

  #remove empty columns
  clinical_data=clinical_data %>%
   mutate_all(funs(na_if(., "'--"))) %>%
   remove_empty("cols")

clinical_data = dplyr::filter(clinical_data, clinical_data$case_submitter_id %in% sample_file$Case.ID)

clinical_data <- as.data.frame(clinical_data %>% dplyr::group_by(case_submitter_id) %>% dplyr::filter (! duplicated(case_submitter_id)))


sample_file = dplyr::filter(sample_file, sample_file$Case.ID %in% clinical_data$case_submitter_id)

# merge both files together
rownames(sample_file) = sample_file$Case.ID
rownames(clinical_data)=clinical_data$case_submitter_id


metadata = cbind(sample_file, clinical_data)

metadata[is.na(metadata)] = 0


for (i in 1:(ncol(metadata))){
  metadata[,i] = as.factor(metadata[,i])
}
 

class(metadata$case_id)

## I have unziped them, so remove the ".gz" 

metadata <- metadata %>%
  dplyr::mutate_at("File.Name", str_replace, ".gz", "")

### Parallelization 
library("BiocParallel")
register(MulticoreParam(4))

### Create DESeq file for normalization
#### design ~ race + gender + ajcc_pathologic_t 
    ## normalization for gender and race and comparison for pathologic_t

data = DESeqDataSetFromHTSeqCount(sampleTable = metadata, directory = "./GDCdata", design = ~ gender + race + ajjc_pathologic_t)
 
### remove rows with very low number of reads
keep <- rowSums(counts(data)) >= 10
data <- data[keep,]
dim(data)


data$ajcc_pathologic_n <- relevel(data$ajcc_pathologic_t, ref = "T1")



### main DESeq (normalize counts for gene lenght and depth of read)
  ### log2 tables + Wald test for p-value
  ### comparison with the LAST variable of the design : "race"
  ### reference level = 
data = DESeq(data, test=c("Wald"))

########### check results
normCounts <-  counts(data, normalized=T)
View(normCounts)

res <- DESeq2::results(data, alpha=0.05)

summary(res)

res = res[order(res$padj),]

##annotation

annotation= read.


### is ethnicity/sex/age involved in the question


##### DESeq2 normalization


featureData <- data.frame(gene=rownames(data))








dds <- DESeq(data)
res <- results(dds)
res


plotPCA(data, intgroup = "ajcc_pathologic_n", ntop = 500, returnData = FALSE)




data$ajcc_pathologic_stage <- relevel(data$ajcc_pathologic_stage, ref = "untreated")

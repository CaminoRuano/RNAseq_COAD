




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
register(SnowParam(4))

### Create DESeq file for normalization
#### design ~ race + gender + ajcc_pathologic_tm
## normalization for gender and race and comparison for pathologic_m

data = DESeqDataSetFromHTSeqCount(sampleTable = metadata, directory = "./GDCdata", design =~gender + race + treatment_or_therapy)


###  Metadata survey

  ## PCA of the raw data

   library('PCAtools')

    vst = assay(vst(data))
      p = pca(vst, metadata=colData(data), removeVar=0.1)
      pairsplot(p, components= getComponents(p, c(1:6)), colby='treatment_or_therapy"')

      biplot(p, showLoadings=FALSE, colby='ethnicity', lab=paste0(p$metadata$treatment_or_therapy))

      eigencorplot(p, components= getComponents(p, 1:20), 
             metavars=c("treatment_or_therapy","ajcc_pathologic_m", "ajcc_pathologic_t",
                        "ajcc_pathologic_n", "ethnicity", "gender", "morphology", 'icd_10_code',
                        "primary_diagnosis", "ajcc_staging_system_edition", "ajcc_pathologic_stage",
                        "year_of_birth", "prior_treatment", "site_of_resection_or_biopsy"),
             corMultipleTestCorrection = 'BH')

  ## Pearson correlation between metadata information

      library(corrplot)

    metadata<- mutate_all(metadata, function(x) as.numeric(as.factor(x)))


      M = cor(metadata)
      res1 <- cor.mtest(metadata, conf.level = .95)

      corrplot(M, na.label=" ",insig = 'blank')



### remove rows with very low number of reads
keep <- rowSums(counts(data)) >= 10
data <- data[keep,]
dim(data)


data$ajcc_pathologic_m <- relevel(data$treatment_or_therapy, ref = "no")



### main DESeq (normalize counts for gene lenght and depth of read)
### log2 tables + Wald test for p-value
### comparison with the LAST variable of the design : "treatment_or_therapy"
### reference level = no


data_norm = DESeq(data, test=c("Wald"))


### PCA observation


vsp = (vst(data_norm))

DESeq2::plotPCA(vsp, intgroup=c("treatment_or_therapy")) + geom_point(size=5) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

########### check results

normCounts <-  counts(data_norm, normalized=T)
head(normCounts)


### perform contrast

res <- DESeq2::results(data_norm, alpha=0.05)

summary(res)

res = res[order(res$padj),]





## plot PCA with vst (variance stabilizing transformations)

vsp = (vst(data_norm))

DESeq2::plotPCA(vsp, intgroup=c("treatment_or_therapy")) + geom_point(size=5) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


#### visualization  and ranking of DEGs

resultsNames(data_norm)

resApeT <- lfcShrink(data_norm, coef="treatment_or_therapy_yes_vs_no", type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)







data$ajcc_pathologic_stage <- relevel(data$ajcc_pathologic_stage, ref = "untreated")





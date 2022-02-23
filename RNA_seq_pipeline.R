

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

# order by adj.p.value
res = res[order(res$padj),]

##annotation



### is ethnicity/sex/age involved in the question


##### DESeq2 normalization


featureData <- data.frame(gene=rownames(data))





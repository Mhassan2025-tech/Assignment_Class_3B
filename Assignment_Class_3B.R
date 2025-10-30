setwd("C:/Users/AHLp/OneDrive/Desktop/AI & Biotech")
list.files()
getwd()


install.packages("BiocManager")
BiocManager::install(version = "3.21", ask = FALSE, force = TRUE)


BiocManager::install("GEOquery", ask = FALSE, force = TRUE)
BiocManager::install("affy", ask = FALSE, force = TRUE)
BiocManager::install("limma", ask = FALSE, force = TRUE)
BiocManager::install("arrayQualityMetrics", ask = FALSE, force = TRUE)
BiocManager::install("AnnotationDbi", ask = FALSE, force = TRUE)
BiocManager::install("hgu133plus2.db", ask = FALSE, force = TRUE)

library(GEOquery)
library(affy)
library(limma)
library(arrayQualityMetrics)
library(AnnotationDbi)
library(hgu133plus2.db)
library(R.utils)


GSE_Data = getGEO(filename = "GSE252168_series_matrix.txt.gz")

Expression_Data = exprs(GSE_Data)

Feature_Data =  fData(GSE_Data)

Phenotype_Data <-  pData(GSE_Data)

sum(is.na(Phenotype_Data$source_name_ch1))

head(Expression_Data[1:5,1:5])

Gse_data = getGEOSuppFiles("GSE252168", baseDir = "Raw_Data", 
                           makeDirectory = TRUE)
untar("Raw_Data/GSE252168/GSE252168_RAW.tar", exdir = "Raw_Data/GSE252168")


arrayQualityMetrics(expressionset = GSE_Data,
                    outdir = "RESULTS/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)


row_median <- apply(Expression_Data,1,median)
head(row_median)

cutoff = median(row_median)

Probes_to_Keep = row_median > cutoff

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 


indx <- row_median > threshold 
filtered_data <- Phenotype_Data[indx, ] 


Processed_Data <- filtered_data 

class(Phenotype_Data$source_name_ch1) 

Labled_Groups <- factor(Phenotype_Data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("Normal", "Cancer_Patient"))

class(Labled_Groups)
levels(Labled_Groups)





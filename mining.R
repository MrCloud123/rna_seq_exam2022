library(DESeq2)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(stats)

## data pre-processing
sampleNames <- c("SRR7821949","SRR7821950","SRR7821951","SRR7821952","SRR7821953","SRR7821968","SRR7821969","SRR7821970")

# The first line is command information, so skip it
data <- read.table("E:\\tox\\matrix.txt", header=TRUE, quote="\t", skip=1)
# We want the count, so we start with the 7th column, here we need blood data, so start with 15th.
names(data)[15:22] <- sampleNames
countData <- as.matrix(data[15:22])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("Blood_WT_Case","Blood_WT_Case","Blood_WT_Case","Blood_WT_Case","Blood_WT_Case","Blood_WT_Control","Blood_WT_Control","Blood_WT_Control"))
rownames(database) <- sampleNames

## Set the group information and build the dds object
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

## The DESeq function was used to estimate the dispersion, and then the res objects were obtained by difference analysis.
dds <- DESeq(dds)
cdds<-DESeq2::counts(dds,normalized=TRUE)
#Zbp1:ENSMUSG00000027514
#Cxcl9:ENSMUSG00000029417
Zbp1_B=cdds["ENSMUSG00000027514",]
Cxcl9_B=cdds["ENSMUSG00000029417",]


## data pre-processing
sampleNames <- c("SRR7821921","SRR7821922","SRR7821918","SRR7821919","SRR7821920","SRR7821937","SRR7821938","SRR7821939")

# The first line is command information, so skip it
data <- read.table("E:\\tox\\matrix.txt", header=TRUE, quote="\t", skip=1)
# We want the count, so we start with the 7th column, here we need lung's data, so start with 7th.
names(data)[7:14] <- sampleNames
countData <- as.matrix(data[7:14])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Control","Lung_WT_Control","Lung_WT_Control"))
rownames(database) <- sampleNames

## Set the group information and build the dds object
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

## The DESeq function was used to estimate the dispersion, and then the res objects were obtained by difference analysis.
dds <- DESeq(dds)
cdds<-DESeq2::counts(dds,normalized=TRUE)

Zbp1_L=cdds["ENSMUSG00000027514",]
Cxcl9_L=cdds["ENSMUSG00000029417",]

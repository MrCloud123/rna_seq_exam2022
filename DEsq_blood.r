library(DESeq2)
library(pheatmap)  # A package for heat-map
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
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
rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, intgroup = "condition")
res <- results(dds)
write.csv(res, "res_des_output_Blood.csv")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, "all_des_output_Blood.csv", row.names=FALSE)

# res format conversion: Use data.frame to convert to table form
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# Sort by p-value log2FoldChange
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),]      #  Gene with a significant increase in expression
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),]    # Genes with significantly reduced expression
res1_total <- rbind(res1_up,res1_down)
df <- countData[intersect(rownames(countData),rownames(res1_total)),] 
df2<- as.matrix(df)                                                 
pheatmap(df2,
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=T,
         height=10,  
         scale = "row",
         frontsize = 10,
         angle_col=45, 
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
         clustering_method = 'single',
) 




genes<- res1
# Add color information according to up-regulated, down-regulated, and unchanged gene expression
genes$color <- ifelse(genes$padj<0.05 & abs(genes$log2FoldChange)>= 1,ifelse(genes$log2FoldChange > 1,'red','blue'),'gray')
color <- c(red = "red",gray = "gray",blue = "blue")

p <- ggplot(
  genes, aes(log2FoldChange, -log10(padj), col = color)) +  
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  # subline
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6)
(p)



d <- read.csv("res_des_output_Blood.csv")
d<-subset(d,abs(d[,3])>2)
d<-subset(d,abs(d[,6])<0.05)
d<-subset(d,abs(d[,7])<0.05)
d<-subset(d,abs(d[,2])>10)
## 1st column is ID (no duplicated allowed)
## 2nd column is fold change
## feature 1: numeric vector
geneList <- d[,3]
## feature 2: named vector
names(geneList) <- as.character(d[,1])
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)
ls=as.vector(names(geneList))
ego<-gseGO(geneList = geneList,keyType = "ENSEMBL",OrgDb = org.Mm.eg.db,ont="CC",minGSSize = 100,maxGSSize = 500,verbose = FALSE)
eg <- bitr(ls, fromType="ENSEMBL",toType=c("UNIPROT","GENENAME"),OrgDb="org.Mm.eg.db")

#barplot(ego,showCategory=10)
#dotplot(ego,split=".sign",showCategory=15)
#plotGOgraph(ego)
#goplot(ego)
#cnetplot(ego,showCategroy=5)




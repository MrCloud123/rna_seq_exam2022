library(DESeq2)
library(pheatmap)  # 用于作热图的包
library(ggplot2)   # 用于作图的包
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)


## data pre-processing
sampleNames <- c("SRR7821921","SRR7821922","SRR7821918","SRR7821919","SRR7821920","SRR7821937","SRR7821938","SRR7821939")

# The first line is command information, so skip it
data <- read.table("C:\\Users\\11951\\Desktop\\tox\\matrix.txt", header=TRUE, quote="\t", skip=1)
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
res <- results(dds)

write.csv(res, "res_des_output_Lung.csv")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, "all_des_output_Lung.csv", row.names=FALSE)

# res format conversion: Use data.frame to convert to table form
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# Sort by p-value log2FoldChange
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1_up<- res1[which(res1$log2FoldChange >= 1 & res1$pvalue < 0.05),]      # Gene with a significant increase in expression
res1_down<- res1[which(res1$log2FoldChange <= -1 & res1$pvalue < 0.05),]    # Genes with significantly reduced expression
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

library(stats)

d <- read.csv("res_des_output_Lung.csv")
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
dotplot(ego,split=".sign",showCategory=15)
#plotGOgraph(ego)
#goplot(ego)
#cnetplot(ego,showCategroy=5)








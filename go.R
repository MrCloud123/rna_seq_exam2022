library(org.Mm.eg.db)
library(stats)
library(clusterProfiler)

d <- read.csv("res_des_output_lung.csv")
#d <- read.csv("res_des_output_Lung.csv")
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
ego<-gseGO(geneList = geneList,keyType = "ENSEMBL",OrgDb = org.Mm.eg.db,ont="ALL",minGSSize = 100,maxGSSize = 500,verbose = FALSE)
eg <- bitr(ls, fromType="ENSEMBL",toType=c("UNIPROT","GENENAME"),OrgDb="org.Mm.eg.db")


#barplot(ego,showCategory=10)
dotplot(ego,split=".sign",showCategory=10)
#plotGOgraph(ego)
#goplot(ego)
#cnetplot(ego,showCategroy=5)

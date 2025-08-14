#######change gene ID and get matrix ######
rm(list = ls())

setwd("/Users/fraya/Documents/project/NRCM_TA/RNAseq/Gclc_OE/")
counts<-read.table("./01_counts/join.count",header=T)
head(counts)
tail(counts)
#rownames(counts) <- counts
library("dbplyr")
library('biomaRt')
library("curl")
library(ggrepel)
library(dbplyr)
my_mart <-useMart("ensembl")
list<-listDatasets(my_mart)
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ensembl"))
listAttributes(mart)
my_ensembl_gene_id<-counts$gene
head(my_ensembl_gene_id)
options(timeout = 4000000)
#####use ensembl_transcript_id to trans######
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',
                                 'ensembl_transcript_id',"description"),
                    filters = 'ensembl_gene_id', 
                    values= my_ensembl_gene_id, mart = mart)
mms_symbols$gene<-mms_symbols$ensembl_gene_id
counts<-merge(counts,mms_symbols,by="gene")
str(counts)
#counts<-counts[,c(2:7,9)]
counts<-counts[!duplicated(counts$external_gene_name), ]
rownames(counts)<-counts$external_gene_name
head(counts)
counts<-counts[,c(8:10,5:7,2:4)]
colnames(counts)<-c("Nor-rep1","Nor-rep2","Nor-rep3",
                    "DMSO-rep1","DMSO-rep2","DMSO-rep3",
                    "Gclc-rep1","Gclc-rep2","Gclc-rep3")

write.table(counts,"./01_counts/rawcounts_add_genesymbol.txt")
write.csv(counts,"./01_counts/rawcounts_add_genesymbol.csv")

# global raw_counts pca
library(ggpubr)
library(ggthemes)
library(gmodels)
library(export)
counts<- read.table("./01_counts/rawcounts_add_genesymbol.txt")
head(counts)
pca.info <- fast.prcomp(counts)
head(pca.info)
summary(pca.info) 
head(pca.info$rotation) 

pca.data <- data.frame(sample = rownames(pca.info$rotation), 
                       Type=c(rep("Nor",3),
                              rep("DMSO",3),
                              rep("Gclc",3)), pca.info$rotation)


p<-ggscatter(pca.data,x="PC1", y="PC2", color="Type", ellipse=TRUE, ellipse.type="convex")
p


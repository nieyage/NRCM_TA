# compare with hypoxia injury
library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE")
TA_counts<-read.table("../counts/rawcounts_add_genesymbol.txt",header=T,row.names=1)
head(TA_counts)
PPOE_counts<- read.table("./counts/rawcounts_only_analysis.txt",header=T,row.names=1)
head(PPOE_counts)

TA_counts$gene<- rownames(TA_counts)
PPOE_counts$gene<- rownames(PPOE_counts)

count<-merge(TA_counts,PPOE_counts,by="gene")
head(count)
rownames(count)<- count$gene
count<- count[,-1]
count<- count[,-16:-18]
batch <- factor(c(rep("batch1",12),rep("batch2",12)))

treat <- factor(c(rep("CT",3),rep("DMSO",3),rep("H_R10",3),rep("R5",3),
                  rep("Nor",3),rep("DMSO",3),rep("PPOE",3),rep("3PGZ",3)))

colData <- data.frame(row.names=colnames(count),treat=treat,batch=batch)
colData
countData <- count[apply(count,1,sum) > 1 , ] 
head(countData)

dds<-DESeqDataSetFromMatrix(countData,colData, formula(~treat+batch)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
# log trans results######
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
#pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
#hc <- hclust(t(rlogMat) )

library(ggplot2)
vsd <- vst(dds)
head(vsd)
pca_data <- plotPCA(vsd, intgroup=c("treat","batch"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, shape =batch,color=treat)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+theme_bw()

####remove batch####
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), c(colData$batch))
pca <- plotPCA(vsd, intgroup=c("batch","treat"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca, "percentVar"))
p<-ggplot(pca, aes(PC1, PC2,shape =batch,color=treat)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+#geom_text(label=paste(pca_data$treat),colour="black",size=4)+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
adjusted_counts <- sva::ComBat_seq(as.matrix(last), batch = batch, group = type)

head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"RNAseq-PCAcount.csv")


#3D PCA#
library("FactoMineR")
library("factoextra")
library("scatterplot3d")
library("gmodels")
pca_count<- normalized_counts

pca.info <- fast.prcomp(pca_count)
head(pca.info$rotation)
head(pca.info)
pca.data <- data.frame(sample =rownames(pca.info$rotation),
                       condition=treat,
                       pca.info$rotation)

library(ggsci)
library("scales")
colors=pal_npg("nrc")(10)
show_col(pal_npg("nrc")(10))
colors_pal<-colors[c(3,2,5,1,6,7)]
colors <- colors_pal[as.factor(pca.data$condition)]

pVar <- pca.info$sdev^2/sum(pca.info$sdev^2)
pVar = round(pVar,digits = 3)

paste0("PC1 (",as.character(pVar[1] * 100 ),"%)")
paste0("PC2 (",as.character(pVar[2] * 100 ),"%)")
paste0("PC3 (",as.character(pVar[3] * 100 ),"%)")

str(pca.data)
s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],pch=15,color = colors,mar=c(5,5,5,5),
                     angle = 60, type="p",cex.symbols = 1,
                     main = "NRCM_TA 3D PCA plot",
                     xlab="PC1 ",
                     ylab = "PC2",
                     zlab = "PC3 ") 
legend("topleft", legend = c("3PGZ","DMSO","DMSO_up","H_R10","Nor","PPOE","R5"),
       col =colors_pal, pch = 15, bty="n",
       inset = -0.2,xpd = TRUE, horiz = FALSE)

###sample corealation heatmap####
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)[,c(1:12,15:26)]
#rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# pearson correlation
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
library(amap)

condition <- factor(c(rep("WT",12),rep("MEIS1 KO",12)))
timepoint <- factor(c(rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3),
                      rep("Day0",3),rep("Day2",3),rep("Day5",3),rep("Day10",3)))

annotation_col = data.frame(
  Group = condition, 
  Stage = timepoint
)

rownames(annotation_col) = factor(colnames(pearson_cor))

ann_colors= list(
  Stage = c("Day0"="white","Day2"="#B3E3E3","Day5"="#66C3A5","Day10"="#28A45F"),
  Group = c("WT"="#F39B7FB2","MEIS1 KO" ="#8491B4B2"))

library(pheatmap)
pheatmap(pearson_cor,
         cluster_cols = T,cluster_rows = T,
         color = hmcol,
         border_color = NA,
         #annotation_col = annotation_col, 
         #annotation_row = annotation_col,
         #annotation_colors = ann_colors,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)



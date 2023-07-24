library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("/Users/fraya/Documents/project/NRCM_TA")
counts<-read.table("./counts/rawcounts_add_genesymbol.txt",header=T,row.names=1)
head(counts)
str(counts)

condition <- factor(c(rep("CT",3),rep("DMSO",3),rep("H_R10",3),rep("R5",3)))
colData <- data.frame(row.names=colnames(counts),condition=condition)
colData
countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts,"NRCM_TA-normalized_counts.csv")
dds <- DESeq(dds)
res <- results(dds)
# condition in HR global PCA
vsd <- vst(dds)
rld <- rlog(dds)
library(ggplot2)
vsd <- varianceStabilizingTransformation(dds)
head(vsd)

###sample corealation heatmap####
rlogMat <- assay(rld)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
library(amap)
library(heatmap.2)
pheatmap(pearson_cor,
         cluster_cols = T,cluster_rows = T,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)


pca_data <- plotPCA(vsd, intgroup=c("condition"),
                    returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
head(pca_data)
pdf("./condition in HR global PCA.pdf",width=8,height=6)
ggplot(pca_data, aes(PC1, PC2, color =condition)) +
  geom_point(size=3) +geom_text(label=paste(pca_data$name),colour="black",size=4)+
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+theme(panel.grid.major=element_line(colour=NA))
dev.off()


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
                       condition=condition,
                       pca.info$rotation)

library(ggsci)
library("scales")
colors=pal_npg("nrc")(10)
show_col(pal_npg("nrc")(10))
colors_pal<-colors[c(3,2,5,1)]
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
legend("topright", legend = c("CT","DMSO","H_R_10","R5"),
       col =colors_pal, pch = 15, bty="n",
       inset = -0.2,xpd = TRUE, horiz = FALSE)

DMSO_CT <-results(dds,contrast = c("condition","DMSO","CT"))
H_R10_DMSO <-results(dds,contrast = c("condition","H_R10","DMSO"))
R5_DMSO <-results(dds,contrast = c("condition","R5","DMSO"))
# save the DEG results 
write.csv(DMSO_CT,"DMSOvsCT-DEG.csv")
write.csv(H_R10_DMSO,"H_R10_DMSO-DEG.csv")
write.csv(R5_DMSO,"R5_DMSO-DEG.csv")

DEG<- rownames(subset(DMSO_CT, padj < 0.05 & abs(log2FoldChange) >1))

rlogMat <- assay(rld)[DEG,]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
pheatmap(pearson_cor,
         cluster_cols = T,cluster_rows = T,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)

pearson_cor




library(pheatmap)
library(clusterProfiler)
library(org.Rn.eg.db)
#DMSO_CT
#volcano plot#
data<-as.data.frame(DMSO_CT)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'DMSO','Control'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
pdf("DMSO-CT-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("DMSO vs Control DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(DMSO_CT, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),1:3]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("DMSO-CT-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(DMSO_CT, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(DMSO_CT, pvalue < 0.01 & log2FoldChange< -1) )
pdf("DMSOvsCT-up_in_DMSO-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-CC.csv")
dev.off()
######KEGG##########
pdf("DMSOvsCT-up_in_DMSO-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"DMSOvsCT-up_in_DMSO-KEGG.csv")
dev.off()

pdf("DMSOvsCT-up_in_DMSO-GO.pdf")
gene.df <- bitr(down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"DMSOvsCT-up_in_DMSO-CC.csv")
dev.off()
######KEGG##########
pdf("DMSOvsCT-up_in_DMSO-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"DMSOvsCT-up_in_DMSO-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(DMSO_CT),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=DMSO_CT$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(DMSO_CT)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-DMSOvsCT-FSEA.csv")
pdf("GSVA-DMSOvsCT-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
pdf("GSVA-DMSOvsCT-gene.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ?????Ì¨???????????????t????????????????2??????????p????
dev.off()


#H_R10_DMSO
#volcano plot#
data<-as.data.frame(H_R10_DMSO)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'H_R10','DMSO'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
pdf("H_R10_DMSO-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("H_R10 vs DMSO DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(H_R10_DMSO, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),1:3]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("H_R10_DMSO-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(H_R10_DMSO, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(H_R10_DMSO, pvalue < 0.01 & log2FoldChange< -1) )
pdf("H_R10vsDMSO-up_in_H_R10-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R10vsDMSO-up_in_H_R10-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R10vsDMSO-up_in_H_R10-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R10vsDMSO-up_in_H_R10-GO-CC.csv")
dev.off()
######KEGG##########
pdf("H_R10vsDMSO-up_in_H_R10-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"H_R10vsDMSO-up_in_H_R10-KEGG.csv")
dev.off()

pdf("H_R10vsDMSO-up_in_DMSO-GO.pdf")
gene.df <- bitr(down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R10vsDMSO-up_in_DMSO-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R10vsDMSO-up_in_DMSO-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"H_R10vsDMSO-up_in_DMSO-CC.csv")
dev.off()
######KEGG##########
pdf("H_R10vsDMSO-up_in_DMSO-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"H_R10vsDMSO-up_in_DMSO-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(H_R10_DMSO),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=H_R10_DMSO$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(H_R10_DMSO)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-H_R10vsDMSO-FSEA.csv")
pdf("GSVA-H_R10vsDMSO-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
pdf("GSVA-H_R10vsDMSO-gene.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ?????Ì¨???????????????t????????????????2??????????p????
dev.off()

#R5_DMSO
#volcano plot#
data<-as.data.frame(R5_DMSO)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'R5','DMSO'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
pdf("R5_DMSO-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","red", "grey"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("R5 vs DMSO DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(R5_DMSO, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),c(4:6,10:12)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("R5_DMSO-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(R5_DMSO, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(R5_DMSO, pvalue < 0.01 & log2FoldChange< -1) )
pdf("R5vsDMSO-up_in_R5-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"R5vsDMSO-up_in_R5-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"R5vsDMSO-up_in_R5-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"R5vsDMSO-up_in_R5-GO-CC.csv")
dev.off()
######KEGG##########
pdf("R5vsDMSO-up_in_R5-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"R5vsDMSO-up_in_R5-KEGG.csv")
dev.off()

pdf("R5vsDMSO-up_in_DMSO-GO.pdf")
gene.df <- bitr(down, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"R5vsDMSO-up_in_DMSO-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"R5vsDMSO-up_in_DMSO-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"R5vsDMSO-up_in_DMSO-CC.csv")
dev.off()
######KEGG##########
pdf("R5vsDMSO-up_in_DMSO-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"R5vsDMSO-up_in_DMSO-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(R5_DMSO),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=R5_DMSO$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(R5_DMSO)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-R5vsDMSO-FSEA.csv")
pdf("GSVA-R5vsDMSO-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
pdf("GSVA-R5vsDMSO-gene.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # dev.off()
H_R10_DEG<- rownames(subset(H_R10_DMSO, padj < 0.05 & abs(log2FoldChange) >1))
R5_DEG<- rownames(subset(R5_DMSO, padj < 0.05 & abs(log2FoldChange) >1))
DEG<- rownames(subset(DMSO_CT, padj < 0.05 & abs(log2FoldChange) >1))
union_DEG<- union(H_R10_DEG,R5_DEG)

normalized_counts<- read.csv("NRCM_TA-normalized_counts.csv",row.names = 1)
# heatmap 
Group = factor(c(rep("CT",3),rep("DMSO",3),rep("R5",3),rep("H_R10",3)))
anndf <- data.frame(Group)
rownames(anndf) <- colnames(normalized_counts)
#自定义分组颜色条的颜色；
anncol = list(Group=c(CT="#40513B",DMSO="#EDF1D6",H_R10="#609966",R5="#9DC08B"))

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%DEG),]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#pdf("CT-DMSO-H-R10-heatmap_allDEGin_HR10.pdf",width=6,height=10)
count<-na.omit(count)
p<-pheatmap(count,cluster_cols = F,cluster_rows = T,
            color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
            #cellwidth = 10, cellheight = 10,
            cutree_rows = 4,
            annotation_col=anndf,
            annotation_colors=anncol,
            show_rownames=F,show_colnames=T)
row_cluster <- cutree(p$tree_row,k=4)
table(row_cluster)
annotation_row <- data.frame(
  type = paste("Cluster",row_cluster,sep="")
)
rownames(annotation_row) <- names(row_cluster);
pdf("All-heatmap_DEG_injury.pdf",width=5,height=10)
#bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         cutree_rows = 4,
         #legend_breaks=seq(-2,2,1),
         #breaks=bk,
         annotation_col=anndf,
         annotation_row=annotation_row,
         annotation_colors=anncol,
         show_rownames=F,show_colnames=T)
dev.off()

#GO and KEGG 
for (i in unique(annotation_row$type)){
  print(i);
  i=3
  gene_in_cluster<-rownames(annotation_row)[which(annotation_row$type=="Cluster3")];
  
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("heatmap",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    type="BP"
    ego <- enrichGO(gene.df$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Rn.eg.db,
                    ont = type,
                    pAdjustMethod  = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      p_ego_go <- goplot(ego)
      
      print(p);
      print(p_ego_go);
      write.csv(ego,paste0("heatmap_",i,"_GO_",type,".csv"))
    }
  }
  dev.off()
  #KEGG 
  pdf(paste0("heatmap_",i,"_KEGG.pdf",sep=""))
  ego <- enrichKEGG(
    gene = gene.df$ENTREZID,
    keyType = "kegg",
    organism  = "rno",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  data<-as.data.frame(ego)
  if(!nrow(data)==0){
    p<-barplot(ego, showCategory=20);
    print(p);
    edox<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
    write.table(edox,paste0("heatmap_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
}

kegg<-read.csv("DMSOvsCT-up_in_CT-GO-BP.csv")
#kegg <- mutate(kegg, ratio = parse_ratio(GeneRatio))
kegg <- kegg[1:10,]
kegg$Description = factor(kegg$Description,levels = rev(kegg$Description))
ggplot(data = kegg,aes(x = Count, y = reorder(Description,Count)))+
  geom_point(aes(size = Count,color = -log10(pvalue)))+
  theme_bw()+
  scale_colour_gradient(low = "green",high = "red")+
  scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
  labs(x = "Count",y = "",title = "",
       color = expression(-log10(pvalue)),size = "Count")+
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 11),legend.text = element_text(size = 10))+theme_bw()

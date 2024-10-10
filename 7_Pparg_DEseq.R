library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE")
counts<-read.table("./counts/rawcounts_add_genesymbol.txt",header=T,row.names=1)
head(counts)
str(counts)
counts<- counts[,-1]
counts<- counts[,-4:-6]

write.table(counts,"./counts/rawcounts_only_analysis.txt")
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

condition <- factor(c(rep("Nor",3),
                      #rep("DMSO_down",3),
                      rep("DMSO",3),
                      rep("PPOE",3),
                      rep("3PGZ",3)))
colData <- data.frame(row.names=colnames(counts),condition=condition)
colData
countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts,"NRCM_TA_PPOE-normalized_counts.csv")
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
         cluster_cols = F,cluster_rows = F,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)


pca_data <- plotPCA(vsd, intgroup=c("condition"),
                    returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
head(pca_data)
pdf("./PPOE global PCA.pdf",width=8,height=6)
ggplot(pca_data, aes(PC1, PC2, color =condition)) +
  geom_point(size=3) +#geom_text(label=paste(pca_data$name),colour="black",size=4)+
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA))
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
colors_pal<-colors[c(3,2,5,1,6)]
colors <- colors_pal[as.factor(pca.data$condition)]

 pVar <- pca.info$sdev^2/sum(pca.info$sdev^2)
 pVar = round(pVar,digits = 3)

   paste0("PC1 (",as.character(pVar[1] * 100 ),"%)")
   paste0("PC2 (",as.character(pVar[2] * 100 ),"%)")
   paste0("PC3 (",as.character(pVar[3] * 100 ),"%)")

str(pca.data)
s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],pch=15,color = colors,mar=c(5,5,5,5),
                     angle = 60, type="p",cex.symbols = 1,
                     main = "NRCM PPOE 3D PCA plot",
                     xlab="PC1 ",
                     ylab = "PC2",
                     zlab = "PC3 ") 
legend("topleft", legend = c("3PGZ","DMSO+","DMSO-","Nor","PPOE"),
       col =colors_pal, pch = 15, bty="n",
       inset = -0.2,xpd = TRUE, horiz = FALSE)

DMSO_Nor <-results(dds,contrast = c("condition","DMSO","Nor"))
PPOE_DMSO <-results(dds,contrast = c("condition","PPOE","DMSO"))
PGZ_DMSO <-results(dds,contrast = c("condition","3PGZ","DMSO"))

# save the DEG results 
write.csv(DMSO_Nor,"DMSOvsCT-DEG.csv")
write.csv(PPOE_DMSO,"PPOE_DMSO-DEG.csv")
write.csv(PGZ_DMSO,"PGZ_DMSO-DEG.csv")

DEG<- rownames(subset(DMSO_Nor, pvalue < 0.05 & abs(log2FoldChange) >1))

rlogMat <- assay(rld)[DEG,]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
pheatmap(pearson_cor,
         cluster_cols = T,cluster_rows = T,
         color = hmcol,
         border_color = NA,
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)


library(pheatmap)
library(clusterProfiler)
library(org.Rn.eg.db)
#DMSO_Nor
#volcano plot#
data<-as.data.frame(DMSO_Nor)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'upregulated','downregulated'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
data<- na.omit(data)
pdf("DMSO-Nor-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue","grey","red"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("DMSO vs Normal DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(DMSO_Nor, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),1:6]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("DMSO-Nor-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(DMSO_Nor, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(DMSO_Nor, pvalue < 0.01 & log2FoldChange< -1) )
pdf("DMSOvsNor-up_in_DMSO-GO.pdf")
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
write.csv(ego,"DMSOvsNor-up_in_DMSO-GO-BP.csv")
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
write.csv(ego,"DMSOvsNor-up_in_DMSO-GO-MF.csv")
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
write.csv(ego,"DMSOvsNor-up_in_DMSO-GO-CC.csv")
dev.off()
######KEGG##########
pdf("DMSOvsNor-up_in_DMSO-KEGG.pdf")
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
write.table(ego,"DMSOvsNor-up_in_DMSO-KEGG.csv")
dev.off()

pdf("DMSOvsNor-down_in_DMSO-GO.pdf")
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
write.csv(ego,"DMSOvsCT-down_in_DMSO-GO-BP.csv")
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
write.csv(ego,"DMSOvsCT-down_in_DMSO-GO-MF.csv")

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
write.csv(ego,"DMSOvsCT-down_in_DMSO-CC.csv")
dev.off()
######KEGG##########
pdf("DMSOvsNor-down_in_DMSO-KEGG.pdf")

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
write.table(ego,"DMSOvsNor-down_in_DMSO-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(DMSO_Nor),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=DMSO_Nor$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(DMSO_Nor)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("../rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-DMSOvsNor-FSEA.csv")
pdf("GSVA-DMSOvsNor-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
pdf("GSVA-DMSOvsNor-gene.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ?????Ì¨???????????????t????????????????2??????????p????
dev.off()


#PPOE_DMSO
#volcano plot#
data<-as.data.frame(PPOE_DMSO)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 0.5, 
                     ifelse(data$log2FoldChange> 0.5 ,'upregulated','downregualted'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
data<- na.omit(data)
pdf("PPOE_DMSO-volcano.pdf")
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
       title=paste0("PPOE vs DMSO DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
setwd("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE/PPOE_0.5")
diff_gene <-subset(PPOE_DMSO, pvalue < 0.05 & abs(log2FoldChange)> 0.5) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),4:9]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("PPOE_DMSO-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(PPOE_DMSO, pvalue < 0.05 & log2FoldChange> 0.5) )
down<-rownames(subset(PPOE_DMSO, pvalue < 0.05 & log2FoldChange< -0.5) )
pdf("PPOEvsDMSO-up_in_PPOE-GO.pdf")
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
write.csv(ego,"PPOEvsDMSO-up_in_PPOE-GO-BP.csv")
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
write.csv(ego,"PPOEvsDMSO-up_in_PPOE-GO-MF.csv")
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
write.csv(ego,"PPOEvsDMSO-up_in_PPOE-GO-CC.csv")
dev.off()
######KEGG##########
pdf("PPOEvsDMSO-up_in_PPOE-KEGG.pdf")
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
write.table(ego,"PPOEvsDMSO-up_in_PPOE-KEGG.csv")
dev.off()

pdf("PPOEvsDMSO-down_in_PPOE-GO.pdf")
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
write.csv(ego,"PPOEvsDMSO-down_in_PPOE-GO-BP.csv")
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
write.csv(ego,"PPOEvsDMSO-down_in_PPOE-GO-MF.csv")

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
write.csv(ego,"PPOEvsDMSO-down_in_PPOE-CC.csv")
dev.off()
######KEGG##########
pdf("PPOEvsDMSO-down_in_PPOE-KEGG.pdf")

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
write.table(ego,"PPOEvsDMSO-down_in_PPOE-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(PPOE_DMSO),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=PPOE_DMSO$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(PPOE_DMSO)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("../../rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-PPOEvsDMSO-GSEA.csv")
pdf("GSVA-PPOEvsDMSO-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()

pdf("GSVA-PPOEvsDMSO-gene.pdf",width = 10,height = 10)
for (i in 1:78){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ?????Ì¨???????????????t????????????????2??????????p????
dev.off()

#PGZ_DMSO
#volcano plot#
data<-as.data.frame(PGZ_DMSO)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'upregulated','downregulated'),
                     'not sig')
DEG_data<-as.data.frame(table(data$change))
data<- na.omit(data)
pdf("PGZ_DMSO-volcano.pdf")
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
       title=paste0("PGZ vs DMSO DEG:","    ",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[2,1],":",DEG_data[2,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

diff_gene <-subset(PGZ_DMSO, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),c(4:6,10:12)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("3PGZ_DMSO-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(PGZ_DMSO, pvalue < 0.01 & log2FoldChange> 1) )
down<-rownames(subset(PGZ_DMSO, pvalue < 0.01 & log2FoldChange< -1) )
pdf("PGZvsDMSO-up_in_PGZ-GO.pdf")
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
write.csv(ego,"PGZvsDMSO-up_in_PGZ-GO-BP.csv")
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
write.csv(ego,"PGZvsDMSO-up_in_PGZ-GO-MF.csv")
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
write.csv(ego,"PGZvsDMSO-up_in_PGZ-GO-CC.csv")
dev.off()
######KEGG##########
pdf("PGZvsDMSO-up_in_PGZ-KEGG.pdf")
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
write.table(ego,"PGZvsDMSO-up_in_PGZ-KEGG.csv")
dev.off()

pdf("PGZvsDMSO-down_in_PGZ-GO.pdf")
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
write.csv(ego,"PGZvsDMSO-down_in_PGZ-GO-BP.csv")
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
write.csv(ego,"PGZvsDMSO-down_in_PGZ-GO-MF.csv")

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
write.csv(ego,"PGZvsDMSO-down_in_PGZ-CC.csv")
dev.off()
######KEGG##########
pdf("PGZvsDMSO-down_in_PGZ-KEGG.pdf")

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
write.table(ego,"PGZvsDMSO-down_in_PGZ-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(PGZ_DMSO),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=PGZ_DMSO$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(PGZ_DMSO)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
#kegmt<-read.gmt("rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-PGZvsDMSO-GSEA.csv")
pdf("GSVA-PGZvsDMSO-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()
pdf("GSVA-PGZvsDMSO-gene.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # dev.off()
dev.off()





PPOE_DEG<- rownames(subset(PPOE_DMSO, padj < 0.05 & abs(log2FoldChange) >1))
PGZ_DEG<- rownames(subset(PGZ_DMSO, padj < 0.05 & abs(log2FoldChange) >1))
DEG<- rownames(subset(DMSO_Nor, padj < 0.05 & abs(log2FoldChange) >1))
union_DEG<- union(PPOE_DEG,PGZ_DEG)

normalized_counts<- read.csv("NRCM_TA-normalized_counts.csv",row.names = 1)
# heatmap 
Group = factor(c(rep("CT",3),rep("DMSO",3),rep("PGZ",3),rep("PPOE",3)))
anndf <- data.frame(Group)
rownames(anndf) <- colnames(normalized_counts)
#自定义分组颜色条的颜色；
anncol = list(Group=c(CT="#40513B",DMSO="#EDF1D6",PPOE="#609966",PGZ="#9DC08B"))

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

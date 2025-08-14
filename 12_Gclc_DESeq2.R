library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("/Users/fraya/Documents/project/NRCM_TA/RNAseq/Gclc_OE/")
counts<-read.table("./01_counts/rawcounts_add_genesymbol.txt",header=T,row.names=1)
head(counts)
str(counts)

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

condition <- factor(c(rep("Nor",3),rep("DMSO",3),rep("Gclc",3)))
colData <- data.frame(row.names=colnames(counts),condition=condition)
colData
countData <- counts[apply(counts, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts,"./01_counts/NRCM_TA_Gclc-OE-normalized_counts.csv")
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
library(pheatmap)
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
pdf("./Gclc OE global PCA.pdf",width=5,height=3)
ggplot(pca_data, aes(PC1, PC2, color =condition)) +
  geom_point(size=3) +#geom_text(label=paste(pca_data$name),colour="black",size=4)+
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA))
dev.off()

DMSO_Nor <-results(dds,contrast = c("condition","DMSO","Nor"))
Gclc_DMSO <-results(dds,contrast = c("condition","Gclc","DMSO"))

# save the DEG results 
write.csv(DMSO_Nor,"Gclc-OE-DMSOvsNor-DEG.csv")
write.csv(Gclc_DMSO,"Gclc-OE-GclcvsDMSO-DEG.csv")

DEG<- rownames(subset(DMSO_Nor, padj < 0.01 & abs(log2FoldChange) >3))

rlogMat <- assay(rld)[DEG,]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
pheatmap(pearson_cor,
         cluster_cols = F,cluster_rows = F,
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
pdf("Gclc_OE_DMSO-Nor-volcano.pdf")
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
write.csv(KEGG,"GSVA-DMSOvsNor-GSEA.csv")
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


#Gclc_DMSO
#volcano plot#
data<-as.data.frame(Gclc_DMSO)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'upregulated','downregualted'),
                     'Stable')
DEG_data<-as.data.frame(table(data$change))
data<- na.omit(data)
pdf("Gclc_DMSO-volcano.pdf")
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-8,8)) +
  #ylim(c(0,100))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (pvalue)",
       title=paste0("Gclc vs DMSO DEG:","\n",DEG_data[1,1],":",DEG_data[1,2],"    ",DEG_data[3,1],":",DEG_data[3,2]))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p+ theme_bw() + theme(panel.grid=element_blank())
dev.off()
#setwd("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE/Gclc_0.5")
diff_gene <-subset(Gclc_DMSO, pvalue < 0.05 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),4:9]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)

pdf("Gclc_DMSO-heatmap.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
head(data)
up<-rownames(subset(Gclc_DMSO, pvalue < 0.05 & log2FoldChange> 1) )
down<-rownames(subset(Gclc_DMSO, pvalue < 0.05 & log2FoldChange< -1) )
pdf("GclcvsDMSO-up_in_Gclc-GO.pdf")
gene.df <- bitr(up, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.2,
                qvalueCutoff = 0.2,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"GclcvsDMSO-up_in_Gclc-GO-BP.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.2,
                qvalueCutoff = 0.2,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"GclcvsDMSO-up_in_Gclc-GO-MF.csv")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.2,
                qvalueCutoff = 0.2,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"GclcvsDMSO-up_in_Gclc-GO-CC.csv")
dev.off()
######KEGG##########
pdf("GclcvsDMSO-up_in_Gclc-KEGG.pdf")
ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.5)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"GclcvsDMSO-up_in_Gclc-KEGG.csv")
dev.off()

pdf("GclcvsDMSO-down_in_Gclc-GO.pdf")
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
write.csv(ego,"GclcvsDMSO-down_in_Gclc-GO-BP.csv")
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
write.csv(ego,"GclcvsDMSO-down_in_Gclc-GO-MF.csv")

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
write.csv(ego,"GclcvsDMSO-down_in_Gclc-CC.csv")
dev.off()
######KEGG##########
pdf("GclcvsDMSO-down_in_Gclc-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'rno',
  pvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.5)
ego
barplot(ego, showCategory=20)
ego<-setReadable(ego,OrgDb=org.Rn.eg.db,keyType = "ENTREZID")
write.table(ego,"GclcvsDMSO-down_in_Gclc-KEGG.csv")
dev.off()

library(enrichplot)
gene=bitr(rownames(Gclc_DMSO),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") 

gene_df <- data.frame(logFC=Gclc_DMSO$log2FoldChange, #????????????foldchange
                      SYMBOL = rownames(Gclc_DMSO)) #?????????????????????????À¨?????????????
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df$logFC
names(geneList)=gene_df$SYMBOL
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("../rattus_norvegicus_brown-rat_gmt2.gmt") #????gmt??????t
names(geneList)<-toupper(names(geneList))
KEGG<-GSEA(geneList,
           pvalueCutoff = 0.45,
           TERM2GENE = kegmt) #GSEA????????
write.csv(KEGG,"GSVA-GclcvsDMSO-GSEA.csv")
pdf("GSVA-GclcvsDMSO-dotplot.pdf",width = 16,height = 8)
dotplot(KEGG,color="pvalue",split=".sign")+facet_wrap(~.sign,scales = "free") #????????????????????????
dev.off()

pdf("GSVA-GclcvsDMSO-gene.pdf",width = 10,height = 10)
for (i in 1:200){
  print(i)
  p<-gseaplot2(KEGG,i,color="red",pvalue_table = T)
  print(p)
}
gseaplot2(KEGG,1:5,color="red",pvalue_table = T) # ?????Ì¨???????????????t????????????????2??????????p????
dev.off()


 Gclc_DEG<- rownames(subset(Gclc_DMSO, padj < 0.05 & abs(log2FoldChange) >1))
DEG<- rownames(subset(DMSO_Nor, padj < 0.05 & abs(log2FoldChange) >1))
union_DEG<- union(Gclc_DEG,DEG)

normalized_counts<- read.csv("./01_counts/NRCM_TA_Gclc-OE-normalized_counts.csv",row.names = 1)
# heatmap 
Group = factor(c(rep("CT",3),rep("DMSO",3),rep("Gclc",3)))
anndf <- data.frame(Group)
rownames(anndf) <- colnames(normalized_counts)
#自定义分组颜色条的颜色；
anncol = list(Group=c(CT="#40513B",DMSO="#EDF1D6",Gclc="#609966"))

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
  gene_in_cluster<-rownames(annotation_row)[which(annotation_row$type==i)];
  
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

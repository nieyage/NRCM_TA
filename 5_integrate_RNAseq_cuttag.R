annotatePeaks.pl Nor-PPARG_IgG_peaks_homer.tmp rn5  > Nor-PPARG_IgG_peaks_rn5.anno.xls 

library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
library(org.Rn.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db) 
txdb <- TxDb.Rnorvegicus.UCSC.rn4.ensGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, downstream=1000)
peakAnno <- annotatePeak("../4_MACS2/Nor-PPARG_IgG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Rn.eg.db")

gene.df <- bitr(as.data.frame(peakAnno)$geneId, fromType = "ENSEMBL",
                toType = c("SYMBOL","ENTREZID"),
                OrgDb = org.Rn.eg.db)

bindinggene<-data$Gene.Name
DEG<- read.csv("/public/home/nieyg/project/NRCM_TA/RNAseq/5_DESeq2/H_R10_DMSO-DEG.csv")
up_gene<- DEG[DEG$pvalue<0.05&DEG$log2FoldChange>1,1]
down_gene<- DEG[DEG$pvalue<0.05&DEG$log2FoldChange< -1,1]

up_overlap<-intersect(bindinggene,up_gene)
down_overlap<-intersect(bindinggene,down_gene)

gene.df <- bitr(down_overlap, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
pdf("NRCM_TA-HR10-downregulatedtarget-GO.pdf")
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"NRCM_TA-HR10-downregulatedtarget-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"NRCM_TA-HR10-downregulatedtarget-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"NRCM_TA-HR10-downregulatedtarget-CC.csv")
dev.off()

pdf("NRCM_TA-HR10-upregulatedtarget-GO.pdf")
gene.df <- bitr(up_overlap, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"NRCM_TA-HR10-upregulatedtarget-BP.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"NRCM_TA-HR10-upregulatedtarget-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
ego
barplot(ego, showCategory=20)
write.csv(ego,"NRCM_TA-HR10-upregulatedtarget-CC.csv")
dev.off()

library(VennDiagram)
library(RColorBrewer)
vennplot<-venn.diagram(
  x = list(up_gene,down_gene,na.omit(bindinggene)),
  category.names = c("HR10_CT UP" , "HR10_CT DOWN" , "Nor_PPARG BINDING"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:3],
  alpha = 0.50,
  output=TRUE
)

pdf("vennplot.pdf")
grid.draw(vennplot)
dev.off()





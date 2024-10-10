annotatePeaks.pl Nor-PPARG_IgG_peaks_homer.tmp rn5  > Nor-PPARG_IgG_peaks_rn5.anno.xls 

library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
library(org.Rn.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db) 

data<- read.csv("/Users/fraya/Documents/project/NRCM_TA/Cuttag/Nor-PPARG_IgG_peaks_rn5.anno.csv")
head(data)
bindinggene<-data$Gene.Name
# PPOE overlap 
DEG<- read.csv("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE/PPOE_DMSO-DEG.csv")
up_gene<- DEG[DEG$pvalue<0.05&DEG$log2FoldChange>0.5,1]
down_gene<- DEG[DEG$pvalue<0.05&DEG$log2FoldChange< -0.5,1]

up_overlap<-intersect(bindinggene,up_gene)
down_overlap<-intersect(bindinggene,down_gene)

HR10_DEG<- read.csv("/Users/fraya/Documents/project/NRCM_TA/RNAseq/H_R10vsDMSO/H_R10_DMSO-DEG.csv")
HR10_DEGup_gene<- HR10_DEG[HR10_DEG$pvalue<0.05&HR10_DEG$log2FoldChange>0.5,1]

library(VennDiagram)
library(RColorBrewer)
vennplot<-venn.diagram(
  x = list(up_gene,HR10_DEGup_gene),
  category.names = c("PPOE UP" , "HR_10 UP"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:2],
  alpha = 0.50,
  output=TRUE
)

pdf("PPOE_HR10_active_vennplot.pdf")
grid.draw(vennplot)
dev.off()



gene.df <- bitr(down_overlap, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
pdf("NRCM_PPOE-downregulatedtarget-GO.pdf")
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
write.csv(ego,"NRCM_PPOE-downregulatedtarget-BP.csv")


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
write.csv(ego,"NRCM_PPOE-downregulatedtarget-MF.csv")

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
write.csv(ego,"NRCM_PPOE-downregulatedtarget-CC.csv")
dev.off()

pdf("NRCM_PPOE-upregulatedtarget-GO.pdf")
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
write.csv(ego,"NRCM_PPOE-upregulatedtarget-BP.csv")

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
write.csv(ego,"NRCM_PPOE-upregulatedtarget-MF.csv")

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
write.csv(ego,"NRCM_PPOE-upregulatedtarget-CC.csv")
dev.off()

library(VennDiagram)
library(RColorBrewer)
vennplot<-venn.diagram(
  x = list(up_gene,down_gene,na.omit(bindinggene)),
  category.names = c("PPOE UP" , "PPOE DOWN" , "Nor_PPARG BINDING"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:3],
  alpha = 0.50,
  output=TRUE
)

pdf("vennplot.pdf")
grid.draw(vennplot)
dev.off()

# 3PGZ overlap 
DEG<- read.csv("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE/PGZ_DMSO-DEG.csv")
up_gene<- DEG[DEG$pvalue<0.05&DEG$log2FoldChange>1,1]
down_gene<- DEG[DEG$pvalue<0.05&DEG$log2FoldChange< -1,1]

up_overlap<-intersect(bindinggene,up_gene)
down_overlap<-intersect(bindinggene,down_gene)

gene.df <- bitr(down_overlap, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
pdf("NRCM_3PGZ-downregulatedtarget-GO.pdf")
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
write.csv(ego,"NRCM_3PGZ-downregulatedtarget-BP.csv")


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
write.csv(ego,"NRCM_3PGZ-downregulatedtarget-MF.csv")

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
write.csv(ego,"NRCM_3PGZ-downregulatedtarget-CC.csv")
dev.off()

pdf("NRCM_3PGZ-upregulatedtarget-GO.pdf")
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
write.csv(ego,"NRCM_3PGZ-upregulatedtarget-BP.csv")

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
write.csv(ego,"NRCM_3PGZ-upregulatedtarget-MF.csv")

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
write.csv(ego,"NRCM_3PGZ-upregulatedtarget-CC.csv")
dev.off()

library(VennDiagram)
library(RColorBrewer)
vennplot<-venn.diagram(
  x = list(na.omit(up_gene),na.omit(down_gene),na.omit(bindinggene)),
  category.names = c("3PGZ UP" , "3PGZ DOWN" , "Nor_PPARG BINDING"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:3],
  alpha = 0.50,
  output=TRUE
)

pdf("vennplot.pdf")
grid.draw(vennplot)
dev.off()








# create GSEA genelist and plot the GSEA enrich plot

PPOE_res<- read.csv("/Users/fraya/Documents/project/NRCM_TA/RNAseq/PpargOE/PPOE_DMSO-DEG.csv")
HR_10_res<- read.csv("/Users/fraya/Documents/project/NRCM_TA/RNAseq/H_R10vsDMSO/H_R10_DMSO-DEG.csv")

HR_10_up<-subset(HR_10_res, pvalue < 0.05 & log2FoldChange> 0.5) 
HR_10_up<- HR_10_up$X
PPOE_up<-subset(PPOE_res, padj < 0.01 & log2FoldChange> 1) 
PPOE_up<- PPOE_up$X
PPOE_res<- PPOE_res[order(PPOE_res$log2FoldChang,decreasing = T),]

library(enrichplot)
geneList <- PPOE_res$log2FoldChange
names(geneList) <- PPOE_res$X
geneList <- sort(geneList, decreasing = T)
geneList <- geneList[geneList != 0]
head(geneList)

HR_10_up_geneset<- data.frame(term = c(rep("HR_10_up",length(HR_10_up))),
                      gene = c(HR_10_up))
head(HR_10_up_geneset)

library(clusterProfiler)
set.seed(123456)
egmt <- GSEA(geneList, TERM2GENE=HR_10_up_geneset, verbose=FALSE, 
             nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gsea_results <- egmt@result
gsea_results2 <- gsea_results[order(gsea_results$enrichmentScore,decreasing=F),]

library(enrichplot)
gseaplot2(egmt, geneSetID = c('HR_10_up'), 
          base_size = 15,subplots = 1:2,pvalue_table = TRUE,
          title = "PPOE vs DMSO in HR10 up Gene")


geneList <- HR_10_res$log2FoldChange
names(geneList) <- HR_10_res$X
geneList <- sort(geneList, decreasing = T)
geneList <- geneList[geneList != 0]
head(geneList)

PPOE_up_geneset<- data.frame(term = c(rep("PPOE_up",length(PPOE_up))),
                      gene = c(PPOE_up))
head(PPOE_up_geneset)

library(clusterProfiler)
set.seed(123456)
egmt <- GSEA(geneList, TERM2GENE=PPOE_up_geneset, verbose=FALSE, 
             nPerm = 10000, minGSSize = 10, maxGSSize = 10000, pvalueCutoff=1)
gsea_results <- egmt@result
gsea_results2 <- gsea_results[order(gsea_results$enrichmentScore,decreasing=F),]

library(enrichplot)
gseaplot2(egmt, geneSetID = c('PPOE_up'), 
          base_size = 15,subplots = 1:2,pvalue_table = TRUE,
          title = "HR_10 vs DMSO in PPOE up Gene")


# overlap

PPOE_up<-subset(PPOE_res, padj < 0.01) 
PPOE_up<- PPOE_up[order(PPOE_up$log2FoldChang,decreasing = T),]
PPOE_up_top100<- PPOE_up[1:100,]$X;
PPOE_up_top200<- PPOE_up[1:200,]$X;
PPOE_up_top300<- PPOE_up[1:300,]$X;
PPOE_up_top400<- PPOE_up[1:400,]$X;
PPOE_up_top500<- PPOE_up[1:500,]$X;


HR_10_up<-subset(HR_10_res, pvalue < 0.05 & log2FoldChange> 0.5) 
HR_10_up<- HR_10_up[order(HR_10_up$log2FoldChang,decreasing = T),]
HR_10_up_top100<- HR_10_up[1:100,]$X;
HR_10_up_top200<- HR_10_up[1:200,]$X;
HR_10_up_top300<- HR_10_up[1:300,]$X;
HR_10_up_top400<- HR_10_up[1:400,]$X;
HR_10_up_top500<- HR_10_up[1:500,]$X;

library(VennDiagram)
library(RColorBrewer)
vennplot<-venn.diagram(
  x = list(PPOE_up_top100,HR_10_up_top100),
  category.names = c("PPOE UP top100" , "HR_10 UP top100"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:2],
  alpha = 0.50,
  output=TRUE
)

pdf("PPOE_top100_HR10_active_vennplot.pdf")
grid.draw(vennplot)
dev.off()

vennplot<-venn.diagram(
  x = list(PPOE_up_top200,HR_10_up_top200),
  category.names = c("PPOE UP top200" , "HR_10 UP top200"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:2],
  alpha = 0.50,
  output=TRUE
)

pdf("PPOE_top200_HR10_active_vennplot.pdf")
grid.draw(vennplot)
dev.off()


vennplot<-venn.diagram(
  x = list(PPOE_up_top300,HR_10_up_top300),
  category.names = c("PPOE UP top300" , "HR_10 UP top300"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:2],
  alpha = 0.50,
  output=TRUE
)

pdf("PPOE_top300_HR10_active_vennplot.pdf")
grid.draw(vennplot)
dev.off()

vennplot<-venn.diagram(
  x = list(PPOE_up_top400,HR_10_up_top400),
  category.names = c("PPOE UP top400" , "HR_10 UP top400"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:2],
  alpha = 0.50,
  output=TRUE
)

pdf("PPOE_top400_HR10_active_vennplot.pdf")
grid.draw(vennplot)
dev.off()

vennplot<-venn.diagram(
  x = list(PPOE_up_top500,na.omit(HR_10_up_top500)),
  category.names = c("PPOE UP top500" , "HR_10 UP top500"),
  filename =NULL,
  fill = brewer.pal(7, "Set2")[1:2],
  alpha = 0.50,
  output=TRUE
)

pdf("PPOE_top500_HR10_active_vennplot.pdf")
grid.draw(vennplot)
dev.off()



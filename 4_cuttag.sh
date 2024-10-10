nohup macs2 callpeak -t ../2_bam/Hyp-PPARG_sorted_rmDup_mapped.bam \
               -c ../2_bam/IgG-IgG_sorted_rmDup_mapped.bam  \
               -f BAMPE \
               -g 2.7e+9 \
               -n Hyp-PPARG_IgG &
nohup macs2 callpeak -t ../2_bam/Nor-PPARG_sorted_rmDup_mapped.bam \
               -c ../2_bam/IgG-IgG_sorted_rmDup_mapped.bam  \
               -f BAMPE \
               -g 2.7e+9 \
               -n Nor-PPARG_IgG &

nohup macs2 callpeak -t ../2_bam/Hyp-PPARG_sorted_rmDup_mapped.bam \
               -c ../2_bam/Nor-PPARG_sorted_rmDup_mapped.bam  \
               -f BAMPE \
               -g 2.7e+9 \
               -n Hyp-PPARG_Nor-PPARG &
nohup macs2 callpeak -c ../2_bam/Hyp-PPARG_sorted_rmDup_mapped.bam \
               -t ../2_bam/Nor-PPARG_sorted_rmDup_mapped.bam  \
               -f BAMPE \
               -g 2.7e+9 \
               -n Nor-PPARG_Hyp-PPARG &

awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./Hyp-PPARG_IgG_peaks.narrowPeak > Hyp-PPARG_IgG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./Nor-PPARG_IgG_peaks.narrowPeak > Nor-PPARG_IgG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./Hyp-PPARG_Nor-PPARG_peaks.narrowPeak > Hyp-PPARG_Nor-PPARG_peaks_homer.tmp
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' ./Nor-PPARG_Hyp-PPARG_peaks.narrowPeak > Nor-PPARG_Hyp-PPARG_peaks_homer.tmp

nohup findMotifsGenome.pl Hyp-PPARG_IgG_peaks_homer.tmp  rn4 Hyp-PPARG_IgG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl Nor-PPARG_IgG_peaks_homer.tmp  rn4 Nor-PPARG_IgG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl Hyp-PPARG_Nor-PPARG_peaks_homer.tmp  rn4 Hyp-PPARG_Nor-PPARG_motifDir_Homer -len 8,10,12  &
nohup findMotifsGenome.pl Nor-PPARG_Hyp-PPARG_peaks_homer.tmp  rn4 Nor-PPARG_Hyp-PPARG_motifDir_Homer -len 8,10,12  &

conda activate deeptools
computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 1000 -a 1000    \
-R /public/home/nieyg/reference/ref_gene/rn4_RefSeq.bed  \
-S Nor-PPARG_CPM_normalized.bw Hyp-PPARG_CPM_normalized.bw IgG-IgG_CPM_normalized.bw \
--skipZeros  -o PPARG-CPM.gz  \
--outFileSortedRegions A3A-bedtools-out.bed
plotHeatmap -m PPARG-CPM.gz  -out PPARG-CPM-heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m PPARG-CPM.gz  -out PPARG-CPM-profile.pdf --plotFileFormat pdf --perGroup --dpi 720 


library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
library(org.Rn.eg.db)
library(clusterProfiler)
library(org.Rn.eg.db) 
txdb <- TxDb.Rnorvegicus.UCSC.rn4.ensGene
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, downstream=1000)

# Nor-PPARG_IgG
peakAnno <- annotatePeak("Nor-PPARG_IgG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Rn.eg.db")
pdf("Nor-PPARG_IgG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of Nor-PPARG_IgG-binding loci\nrelative to TSS")
dev.off()

# GO 
data<- read.csv("/public/home/nieyg/project/NRCM_TA/cuttag/4_MACS2/Nor-PPARG_IgG_peaks_rn5.anno.csv")
gene<- data$Gene.Name
#gene<-as.data.frame(peakAnno)$SYMBOL

gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego_BP <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
ego

select<- c("GO:0007215","GO:0006536","GO:0010906","GO:0001676","GO:0015909","GO:0006631",
    "GO:0006635","GO:2000191","GO:0042542","GO:0006801","GO:0042743","GO:0070301",
    "GO:0042542","GO:0090322","GO:0036293")

ego<- ego[ego$ID%in%select,]

ego<- ego[order(ego$pvalue,decreasing=T),]
ego$Description<- factor(ego$Description,levels=ego$Description)

pdf("Nor-PPARG_IgG_GO-BP-select.pdf")
ggplot(ego,aes(x=Description,y=Count,fill=-1*log10(pvalue))) + 
      geom_bar(stat="identity",position = "dodge") +
      coord_flip() + 
      scale_fill_gradient(low = "blue",high = "red")+
      theme_bw() 
dev.off()



write.csv(ego,"Nor-PPARG_IgG-GO-BP.csv")
write.csv(as.data.frame(peakAnno),"Nor-PPARG_IgG-cuttag-peak.annotation.csv",row.names = F)


# Hyp-PPARG_IgG
peakAnno <- annotatePeak("Hyp-PPARG_IgG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Rn.eg.db")
pdf("Hyp-PPARG_IgG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of Hyp-PPARG_IgG-binding loci\nrelative to TSS")
dev.off()

# GO 
gene<-as.data.frame(peakAnno)$SYMBOL
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
ego
pdf("Hyp-PPARG_IgG_GO-BP.pdf")
barplot(ego, showCategory=20)
dev.off()
write.csv(ego,"Hyp-PPARG_IgG-GO-BP.csv")
write.csv(as.data.frame(peakAnno),"Hyp-PPARG_IgG-cuttag-peak.annotation.csv",row.names = F)



# Hyp-PPARG_Nor-PPARG
peakAnno <- annotatePeak("Hyp-PPARG_Nor-PPARG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Rn.eg.db")
pdf("Hyp-PPARG_Nor-PPARG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of Hyp-PPARG_Nor-PPARG-binding loci\nrelative to TSS")
dev.off()

# GO 
gene<-as.data.frame(peakAnno)$SYMBOL
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)
ego
pdf("Hyp-PPARG_Nor-PPARG_GO-BP.pdf")
barplot(ego, showCategory=20)
dev.off()
write.csv(ego,"Hyp-PPARG_Nor-PPARG-GO-BP.csv")
write.csv(as.data.frame(peakAnno),"Hyp-PPARG_Nor-PPARG-cuttag-peak.annotation.csv",row.names = F)


# Nor-PPARG_Hyp-PPARG
peakAnno <- annotatePeak("Nor-PPARG_Hyp-PPARG_summits.bed", 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Rn.eg.db")
pdf("Nor-PPARG_Hyp-PPARG_summits.pdf")
p1<-plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of Nor-PPARG_Hyp-PPARG-binding loci\nrelative to TSS")
dev.off()


# GO 
gene<-as.data.frame(peakAnno)$SYMBOL
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Rn.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Rn.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)
ego
pdf("Nor-PPARG_Hyp-PPARG_GO-BP.pdf")
barplot(ego, showCategory=20)
dev.off()
write.csv(ego,"Nor-PPARG_Hyp-PPARG-GO-BP.csv")
write.csv(as.data.frame(peakAnno),"Nor-PPARG_Hyp-PPARG-cuttag-peak.annotation.csv",row.names = F)


nohup annotatePeaks.pl Nor-PPARG_IgG_peaks_homer.tmp  rn4 >Nor-PPARG_IgG_peaks_rn4







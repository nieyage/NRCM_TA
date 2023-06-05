# compare with hypoxia injury
library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")
rm(list = ls())
setwd("D:/project/NRCM_TA")
counts<-read.table("./counts/rawcounts_add_genesymbol.txt",header=T,row.names=1)
head(counts)
str(counts)
DEG_in_DMSO_CT   <-read.csv("DMSOvsCT-DEG.csv",row.names = 1)
DEG_in_H_R10_DMSO<-read.csv("./H_R10_DMSO-DEG.csv",row.names = 1)
DEG_in_R5_DMSO   <-read.csv("./R5_DMSO-DEG.csv",row.names = 1)

# the gene cluster up after hypoxia but down in add drup in H_R10
hypoxia_up<- rownames(DEG_in_DMSO_CT[which(DEG_in_DMSO_CT$pvalue<0.05&DEG_in_DMSO_CT$log2FoldChange>1),])
H_R10_down <- rownames(DEG_in_DMSO_CT[which(DEG_in_H_R10_DMSO$pvalue<0.05&DEG_in_H_R10_DMSO$log2FoldChange< -1),])
up_down_overlap<-intersect(hypoxia_up,H_R10_down)
# the gene cluster down after hypoxia but up in add drup in H_R10
hypoxia_down<- rownames(DEG_in_DMSO_CT[which(DEG_in_DMSO_CT$pvalue<0.05&DEG_in_DMSO_CT$log2FoldChange< -1),])
H_R10_up <- rownames(DEG_in_DMSO_CT[which(DEG_in_H_R10_DMSO$pvalue<0.05&DEG_in_H_R10_DMSO$log2FoldChange>1),])
down_up_overlap<-intersect(hypoxia_down,H_R10_up)

normalized_counts<- read.csv("NRCM_TA-normalized_counts.csv",row.names = 1)
# heatmap 
Group = factor(c(rep("CT",3),rep("DMSO",3),rep("H_R10",3)))
anndf <- data.frame(Group)
rownames(anndf) <- colnames(normalized_counts[,1:9])
#自定义分组颜色条的颜色；
anncol = list(Group=c(CT="#40513B",DMSO="#EDF1D6",H_R10="#609966"))
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%c(up_down_overlap,down_up_overlap)),1:9]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
#bk <- c(seq(-1.5,-0.1,by=0.05),seq(0,1.5,by=0.05))
pdf("CT-DMSO-H-R10-heatmap_signif.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         annotation_col=anndf,
         legend_breaks=seq(-1.5,1.5,1),
         #breaks=bk,
         cutree_rows = 2,
         annotation_colors=anncol,
         show_rownames=T,show_colnames=T)
dev.off()

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%c(H_R10_up,H_R10_down)),1:9]
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
pdf("CT-DMSO-H-R10-heatmap_allDEGin_HR10.pdf",width=5,height=10)
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

library(pheatmap)
library(clusterProfiler)
library(STRINGdb) 
library(tidyverse) 
library(org.Rn.eg.db) 
library(igraph) 
library(ggraph)
string_db <- STRINGdb$new(version="11.5",species=10116)
#GO and KEGG 
for (i in unique(annotation_row$type)){
  print(i);
  gene_in_cluster<-rownames(annotation_row)[which(annotation_row$type==i)];
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  #GO
  pdf(paste0("H_R10",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
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
      print(p);
      write.csv(ego,paste0("H_R10_",i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  #KEGG 
  pdf(paste0("H_R10_",i,"_KEGG.pdf",sep=""))
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
    write.table(edox,paste0("H_R10_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
  # Gene-Concept Network
  input<-DEG_in_DMSO_CT$baseMean
  names(input)<-rownames(DEG_in_DMSO_CT)
  p1 <- cnetplot(edox, foldChange=input, cex_label_category=0.8)
  p3 <- cnetplot(edox, foldChange=input, circular = TRUE, colorEdge = TRUE) 
  pdf(paste0("H_R10_",i,"_Gene-Concept_Network.pdf",sep=""))
  print(p1)
  print(p3)
  dev.off()
  # Stringdb 
  data_mapped <- gene.df %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",removeUnmappedRows = TRUE)  
  string_db$plot_network( data_mapped$STRING_id )
  data_links <- data_mapped$STRING_id %>% string_db$get_interactions() 
  links <- data_links %>%
    mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
    mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
    dplyr::select(from, to , last_col()) %>% 
    dplyr::rename(weight = combined_score)
  nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
  net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
  igraph::V(net)$deg <- igraph::degree(net)
  igraph::V(net)$size <- igraph::degree(net)/5 #
  igraph::E(net)$width <- igraph::E(net)$weight/10
  pdf(paste0("H_R10_",i,"_Stringdb-Network.pdf",sep=""))
  p4<-ggraph(net,layout = "kk")+
    geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
    geom_node_point(aes(size=size), color="orange", alpha=0.7)+
    geom_node_text(aes(filter=deg>1, label=name), size = 5, repel = T)+
    scale_edge_width(range = c(0.2,1))+
    scale_size_continuous(range = c(1,10) )+
    guides(size=F)+
    theme_graph()
  print(p4)
  dev.off()
}

 



# the gene cluster up after hypoxia but down in add drup in R5
hypoxia_up<- rownames(DEG_in_DMSO_CT[which(DEG_in_DMSO_CT$pvalue<0.05&DEG_in_DMSO_CT$log2FoldChange>1),])
R5_down <- rownames(DEG_in_DMSO_CT[which(DEG_in_R5_DMSO$pvalue<0.05&DEG_in_R5_DMSO$log2FoldChange< -1),])
up_down_overlap<-intersect(hypoxia_up,R5_down)
# the gene cluster down after hypoxia but up in add drup in R5
hypoxia_down<- rownames(DEG_in_DMSO_CT[which(DEG_in_DMSO_CT$pvalue<0.05&DEG_in_DMSO_CT$log2FoldChange< -1),])
R5_up <- rownames(DEG_in_DMSO_CT[which(DEG_in_R5_DMSO$pvalue<0.05&DEG_in_R5_DMSO$log2FoldChange>1),])
down_up_overlap<-intersect(hypoxia_down,R5_up)

normalized_counts<- read.csv("NRCM_TA-normalized_counts.csv",row.names = 1)
# heatmap 
Group = factor(c(rep("CT",3),rep("DMSO",3),rep("R5",3)))
anndf <- data.frame(Group)
rownames(anndf) <- colnames(normalized_counts[,1:9])
#自定义分组颜色条的颜色；
anncol = list(Group=c(CT="#40513B",DMSO="#EDF1D6",R5="#9DC08B"))
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%c(up_down_overlap,down_up_overlap)),c(1:6,10:12)]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
pdf("CT-DMSO-R5-heatmap_signif.pdf",width=6,height=10)
count<-na.omit(count)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         annotation_col=anndf,
         legend_breaks=seq(-1.5,1.5,1),
         #breaks=bk,
         cutree_rows = 2,
         annotation_colors=anncol,
         show_rownames=T,show_colnames=T)
dev.off()

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%c(R5_up,R5_down)),1:9]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
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
pdf("CT-DMSO-R5-heatmap_allDEGin_R5.pdf",width=5,height=10)
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
  pdf(paste0("R5",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
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
      print(p);
      write.csv(ego,paste0("R5_",i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  #KEGG 
  pdf(paste0("R5_",i,"_KEGG.pdf",sep=""))
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
    write.table(edox,paste0("R5_",i,"_KEGG_",type,".csv"))
  }
  dev.off()
  # Gene-Concept Network
  input<-DEG_in_DMSO_CT$baseMean
  names(input)<-rownames(DEG_in_DMSO_CT)
  p1 <- cnetplot(edox, foldChange=input, cex_label_category=0.8)
  p3 <- cnetplot(edox, foldChange=input, circular = TRUE, colorEdge = TRUE) 
  pdf(paste0("R5_",i,"_Gene-Concept_Network.pdf",sep=""))
  print(p1)
  print(p3)
  dev.off()
}

string_db <- STRINGdb$new(version="11.5",species=10116,score_threshold=400)
#for (i in unique(annotation_row$type)){
  print(i);
  i="Cluster1"
  gene_in_cluster<-rownames(annotation_row)[which(annotation_row$type==i)];
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Rn.eg.db)
  # Stringdb 
  data_mapped <- gene.df %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",removeUnmappedRows = TRUE)  
  pdf(paste0("HR10_",i,"_Stringdb-Network.pdf",sep=""))
  string_db$plot_network( data_mapped$STRING_id )
  dev.off()
  #data_links <- data_mapped$STRING_id %>% string_db$get_interactions() 
  #links <- data_links %>%
  #  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  #  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  #  dplyr::select(from, to , last_col()) %>% 
  #  dplyr::rename(weight = combined_score)
  #nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
  #net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
  #igraph::V(net)$deg <- igraph::degree(net)
  #igraph::V(net)$size <- igraph::degree(net)/5 #
  #igraph::E(net)$width <- igraph::E(net)$weight/10
 
  #p4<-ggraph(net,layout = "kk")+
  #  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  #  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  #  geom_node_text(aes(filter=deg>0.5, label=name), size = 5, repel = T)+
  #  scale_edge_width(range = c(0.2,1))+
  #  scale_size_continuous(range = c(1,10) )+
  #  guides(size=F)+
  #  theme_graph()
  print(p4)
  dev.off()
}

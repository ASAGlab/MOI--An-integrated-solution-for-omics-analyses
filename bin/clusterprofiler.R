library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pathview)
library(enrichplot)
library(DOSE)
library(GOSemSim)
library(ggplot2)
library(optparse)
library(stringr)
library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
hsGO <- godata('org.Hs.eg.db', ont="MF")

#data(geneList)
#gene <- names(geneList)[1:100]
#gene <- sample(geneList,24)


# 
#args<-c("/home/bianca/Desktop/Dockerfiles/clusterprofiler/defeatures.txt", "rp",0.05)
# genes<-read.table("/home/bianca/Desktop/Dockerfiles/clusterprofiler/short_counts", sep="\t", row.names = 1, header=T)
#genes<-read.table("//home/bianca/Downloads/mcia/mcia_results/integrated.txt", sep="\t", header = T, row.names = 1)
genes<-read.table(args[1], sep="\t", row.names = 1, header=T)

#print("first genes")
#print(head(genes))
#genes<-read.table("/home/bianca/Desktop/Dockerfiles/dea/dea.txt", sep="\t", row.names = 1, header=T)
#hsGO <- godata('org.Hs.eg.db', ont="BP", keytype = "SYMBOL")
hsGO<-godata('org.Hs.eg.db',ont="BP")
geneList<-c()
pval<-as.numeric(args[3])
if(args[2]=="edger"){
  genes<-genes[,c(1,5)]
  colnames(genes)<-c("logFC","pval")
  genes<-genes[genes$pval < pval,]
  geneList<-as.numeric(genes[,2])
  names(geneList) = as.character(rownames(genes))
  geneList = sort(geneList, decreasing = TRUE)
  print(head(geneList))
  print(tail(geneList))
}else if(args[2]=="deseq2"){
  genes<-genes[,c(2,6)]
  colnames(genes)<-c("logFC","pval")
  genes<-genes[genes$pval < pval,]
  geneList<-as.numeric(genes[,2])
  names(geneList) = as.character(rownames(genes))
  geneList = sort(geneList, decreasing = TRUE)
}else if(args[2]=="rp"){
  genes<-genes[,c(3,4)]
  colnames(genes)<-c("logFC","pval")
  genes<-genes[genes$pval < pval,]
  geneList<-as.numeric(genes[,1])
  names(geneList) = as.character(rownames(genes))
  geneList = sort(geneList, decreasing = TRUE)
}else if(args[2]=="mcia"){
  # print("mcia genes")
  # genes<-read.table(args[1], sep=" ", row.names = 1, header=T)
  # print(head(genes))
  geneList<-as.numeric(rev(seq_along(genes$var)))
  names(geneList)<-genes$var
}

#gl<-geneList
#print(head(geneList))
#print(args[2])

#library(biomaRt)

homo.anno<-useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl")


check_mimat <- function(geneList) {
  starts_with_mimat <- grepl("^MIMAT", names(geneList))
  return(starts_with_mimat)
}

# bic<-function(geneList){
#   if (any(check_mimat(geneList))) {
#     geneLista<-geneList
#     mirna_rows <- names(geneList)[check_mimat(geneList)]
#     values=gsub("MAT","", mirna_rows)
#     mirna.df=getBM(c("entrezgene_id","mirbase_accession"),"mirbase_accession",values,homo.anno)
#     colnames(mirna.df)=c("ENTREZID","V2")
#     
#     mirna.df<-mirna.df[!duplicated(mirna.df$V2),]
#     rownames(mirna.df)<-gsub("MI","MIMAT",mirna.df$V2)
#     
#     mirna.df<-na.omit(mirna.df)
#     geneList<-na.omit(geneList)
#     mirna.df<-mirna.df[names(geneList),]
#     gsub("MI","MIMAT",mirna.df$V2)==names(geneList)
#     names(geneList)<-mirna.df$ENTREZID
#     geneListmi<-sort(geneList, decreasing=TRUE)
#     return(geneListmi)
#     print("done annotation miRNA")
#     #print(head(geneList))
#     #print(head(geneLista))
#     if((length(names(geneLista)[!check_mimat(geneLista)]) > 0)){
#       print("here")
#       gene.df <- bitr(names(geneLista), fromType = "SYMBOL",
#                       toType = "ENTREZID",
#                       OrgDb = org.Hs.eg.db)
#       
#       
#       
#       
#       gene.df<-gene.df[!duplicated(gene.df$SYMBOL),]
#       rownames(gene.df)<-gene.df$SYMBOL
#       
#       gene.df<-na.omit(gene.df)
#       geneListb<-na.omit(geneList)
#       gene.df<-gene.df[names(geneListb),]
#       gene.df$SYMBOL==names(geneListb)
#       names(geneListb)<-gene.df$ENTREZID
#       allf<-sort(c(geneListmi,geneListb), decreasing = TRUE)
#       #print(head(allf))
#       geneList<-sort(allf, decreasing=TRUE)
#       return(geneList)
#       print("done annotation both")
#     }}else if(!any(check_mimat(geneList))) {
#       print("doing gene annotation")
#       gene.df <- bitr(names(geneList), fromType = "SYMBOL",
#                       toType =  "ENTREZID",
#                       OrgDb = org.Hs.eg.db)
#       
#       gene.df<-gene.df[!duplicated(gene.df$SYMBOL),]
#       rownames(gene.df)<-gene.df$SYMBOL
#       gene.df<-na.omit(gene.df)
#       geneList<-na.omit(geneList)
#       gene.df<-gene.df[names(geneList),]
#       gene.df$SYMBOL==names(geneList)
#       names(geneList)<-gene.df$ENTREZID
#       geneList<-sort(geneList, decreasing=TRUE)
#       
#       return(geneList)
#       print("done annotation genes")
#     }
# }


gene.df<-data.frame("ENTREZID"="ENTREZID","V2"="V2")
mirna.df<-data.frame("ENTREZID"="ENTREZID","V2"="V2")
trans.df<-data.frame("ENTREZID"="ENTREZID","V2"="V2")
trans<-c()
mirna_rows<-c()
### mirna annotation

mirna_rows <- names(geneList)[check_mimat(geneList)]
values=gsub("MAT","", mirna_rows)
if(length(values)!=0){
  mirna.df=getBM(c("entrezgene_id","mirbase_accession"),"mirbase_accession",values,homo.anno)
  colnames(mirna.df)=c("ENTREZID","V2")

  mirna.df<-mirna.df[!duplicated(mirna.df$V2),]
  rownames(mirna.df)<-gsub("MI","MIMAT",mirna.df$V2)

  mirna.df<-na.omit(mirna.df)
  mirna_rows<-na.omit(mirna_rows)
  mirna.df<-mirna.df[mirna_rows,]
  gsub("MI","MIMAT",mirna.df$V2)==mirna_rows

}

#geneList[mirna_rows]<-mirna.df$ENTREZID

### transcript annotation
trans<-names(geneList)[grepl("^ENST", names(geneList))]


if(length(trans)!=0){
  trans.df=getBM(c("entrezgene_id","ensembl_transcript_id"),"ensembl_transcript_id",trans,homo.anno)
  colnames(trans.df)=c("ENTREZID","V2")
  trans.df<-trans.df[!duplicated(trans.df$V2),]
  rownames(trans.df)<-trans.df$V2
  trans.df<-na.omit(trans.df)
  trans=na.omit(trans)
  trans.df<-trans.df[trans,]
  trans.df$V2==trans
}

#geneList[trans.df$V2]<-trans.df$ENTREZID


############## genes annotation
allothers<-c(trans,mirna_rows)
if(length(allothers!=0)){
  print(head(allothers))
  others<-geneList[!(names(geneList) %in% c(trans,mirna_rows))]
  others<-na.omit(others)
  others<-names(others)
  print(head(others))
}else{
  others<-names(geneList)
}

gene.df<-getBM(c("entrezgene_id","hgnc_symbol"),"hgnc_symbol",others,homo.anno)
colnames(gene.df)=c("ENTREZID","V2")
gene.df<-gene.df[!duplicated(gene.df$V2),]

rownames(gene.df)<-gene.df$V2
#head(gene.df)
gene.df<-na.omit(gene.df)
others=na.omit(others)
gene.df<-gene.df[others,]
gene.df$V2==others
#geneList[gene.df$V2]<-gene.df$ENTREZID

##########################################
# names(geneList)<-geneList
# g1eneList<-geneListback[geneList]
# geneList<-na.omit(geneList)
# trans.df<-trans.df[names(geneList),]


#print(head(gene.df))
#print(head(mirna.df))
#print(head(trans.df))
df<-rbind(gene.df,trans.df,mirna.df)
df<-na.omit(df)
geneList<-geneList[df$V2]
names(geneList)==df$V2

names(geneList)<-df$ENTREZID
print(head(geneList))
#names(geneList)<-trans.df$ENTREZID



#geneList<-bic(geneList)
#print(length(geneList))
#print(head(geneList))
#if(length(geneList>2)){
#    print("You have more than 2")
enrichres <- enrichGO(gene          = unique(names(geneList)),
                      #universe      = unique(allf),
                      #keyType = "SYMBOL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
#}

#if(length(geneList>2)){
# ego3 <- gseGO(geneList     =geneList,
#               OrgDb        = org.Hs.eg.db,
#               #keytype="SYMBOL",
#               ont          = "ALL",
#               minGSSize    = 10,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.05,
#               verbose      = FALSE)
#   enrichres <- enrichGO(gene          = unique(names(geneList)),
#                         #universe      = unique(allf),
#                         #keyType = "SYMBOL",
#                         OrgDb         = org.Hs.eg.db,
#                         ont           = "ALL",
#                         pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05,
#                         qvalueCutoff  = 0.05,
#                         readable      = TRUE)
#}

#enrichedgenes<-geneList[names(geneList) %in% unlist(strsplit(ego3@result$core_enrichment,"/"))]
#gseGO(enrichedgenes)


svg("barplot.svg",width = 10, height = 10)
mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>%
  barplot(x = "qscore",showCategory = 10)
dev.off()
svg("dotplot.svg")
dotplot(enrichres, showCategory = 10)
dev.off()
svg("cnetplot.svg")
cnetplot(enrichres)
dev.off()


svg("upsetplot.svg")
upsetplot(enrichres,showCategory = 10)
dev.off()
svg("heatplot.svg")
heatplot(enrichres, showCategory = 10)
dev.off()

enrichres2 <- try(pairwise_termsim(enrichres),silent = TRUE) # calculate pairwise similarities of the enriched terms using Jaccard’s similarity index

svg("ematplot.svg")
emapplot(enrichres2,showCategory = 10)
dev.off()
svg("treeplot.svg")
treeplot(enrichres2,showCategory = 10)
dev.off()



#gene<-rownames(genes)


# if (startsWith(rownames(genes)[1], "ENSG")){
#   genes$names<-gsub("\\..*$", "", rownames(genes))
#   duplicated_rows <- duplicated(genes$names)
#   genes <- aggregate(. ~ names, data = genes, sum)
#   name<-bitr(genes$names, fromType = "ENSEMBL", toType = "SYMBOL",OrgDb="org.Hs.eg.db")
#   rownames(genes)<-name$SYMBOL
#   genes<-genes[,-1]
# }
# print(head(genes))



#semantic<-mgeneSim(rownames(genes), semData=hsGO, measure="Wang", combine="BMA", verbose=FALSE)
#distance_matrix <- 1 - semantic
#hc <- hclust(as.dist(distance_matrix), method = "complete")
#num_clusters <- 5  
#clusters <- as.data.frame(cutree(hc, k = num_clusters))
#colnames(clusters)="clusters"


## feature 1: numeric vector
# geneList = as.numeric(deg[,2])
## genes will have on second column FC



# # for (cluster_id in unique(clusters$clusters)) {
# #   # Extract the differentially expressed genes for the current cluster
# #   de_genes <- clusters[[cluster_id]]
# #   geneList<- geneList[[de_genes]]
# #   
# #   # Perform pathway enrichment analysis using the 'enrichPathway' function
# #   enrichment_results[[as.character(cluster_id)]] <- gseGO(geneList     = gene_list,
# #                                                           OrgDb        = org.Hs.eg.db,
# #                                                           keyType= "SYMBOL",
# #                                                           ont          = "ALL",
# #                                                           minGSSize    = 100,
# #                                                           maxGSSize    = 500,
# #                                                           pvalueCutoff = 1,
# #                                                           verbose      = FALSE)
# # }
# 
# 
# 
# 
# # Perform pathway enrichment analysis for each cluster
# enrichment_results <- list()
# genelist<-rownames(clusters)
# 
# 
# # View the pathway enrichment results for each cluster
# print(enrichment_results)
# 
# 
# # GO over-representation analysis
# ego <- enrichGO(gene          = gene,
#                 universe      = gene,
#                 OrgDb         = org.Hs.eg.db,
#                 keyType= "SYMBOL",
#                 ont           = "ALL",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 1,
#                 qvalueCutoff  = 0,
#                 readable      = TRUE)
# 
# ego <- enrichGO(gene          = gene,
#                 #universe      = gene,
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "ALL",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 1,
#                 qvalueCutoff  = 0,
#                 readable      = TRUE)
# head(ego)
# 
# gene<-read.table("/home/bianca/Desktop/Dockerfiles/dea/dea.txt")
# geneList=gene[,-c(2:dim(gene)[2])]
# 
# names(geneList) = as.character(rownames(gene))
# geneList = sort(geneList, decreasing = TRUE)
# # GO gene set enrichment analysis
# ego3 <- gseGO(geneList     = geneList,
#               OrgDb        = org.Hs.eg.db,
#               keyType= "SYMBOL",
#               ont          = "ALL",
#               minGSSize    = 100,
#               maxGSSize    = 500,
#               pvalueCutoff = 1,
#               verbose      = FALSE)
# ego3 <- gseGO(geneList     = geneList,
#               OrgDb        = org.Hs.eg.db,
#               ont          = "BP",
#               minGSSize    = 100,
#               maxGSSize    = 500,
#               pvalueCutoff = 0.1,
#               verbose      = FALSE)
# svg("Go_geneset.svg")
# goplot(ego3)
# dev.off()
# 
# go_ids<-as.data.frame(ego3@result$ID)
# go_ids<-as.data.frame(ego@result$ID)
# 
# go_ids<-sample(go_ids$`ego3@result$ID`,  24)
# 
# 
# 
# 
# 
# # similarity_list=list()
# # for (i in 1:length(go_ids)) {
# #   
# #   go_id <- go_ids[i]
# #   go_id2<-go_ids[-i]
# #   sim <- mgoSim(go_id, go_id2, hsGO, combine=NULL)
# #   similarity_list[[i]] <- sim
# # }
# 
# # similarity_list
# # threshold <- 0.1
# # for(i in 1:length(similarity_list)){
# #   filtered_list <- apply(as.data.frame(similarity_list[[i]]), 1,function(df) {
# #   print(df)  
# #   row_means <- rowMeans(df[-1])  # Calculate row means excluding the first column
# #   df[row_means > threshold, ]
# #   })
# # }
# # 
# # filtered_list
# # 
# # 
# 
# 
# 
# 
# 
# similarity_matrix <- matrix(NA, nrow = length(go_ids), ncol = length(go_ids))
# rownames(similarity_matrix) <- go_ids
# colnames(similarity_matrix) <- go_ids
# for (i in 1:length(go_ids)) {
#   
#   go_id <- go_ids[i]
#   go_id2<-go_ids[-i]
#   sim <- mgoSim(go_id, go_id2, hsGO)
#   similarity_matrix[i,] <- sim
#   
# }
# similarity_matrix
# 
# threshold <- 0.1
# 
# row_means <- rowMeans(similarity_matrix)
# 
# filtered_matrix <- similarity_matrix[row_means > threshold, ]
# 
# goOs<-c(colnames(filtered_matrix, rownames(filtered_matrix)))
# 
# ego3@result$ID<-ego3@result$ID %in% goOs
# clusterProfiler::plotGOgraph(ego3)
# svg("gooo")
# ggplot(ego3)
# dev.off()

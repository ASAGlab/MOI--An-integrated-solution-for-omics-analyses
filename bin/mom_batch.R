# batch effect corrections
library(sva)
library(ggplot2)
library(tidyverse)
library(optparse)

args <- commandArgs(trailingOnly=TRUE)

### color palettes
# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

#sampleInfo = DATA FRAME samples condition batch
comsva<- function(samples,sampleInfo,condition,batch){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + condition, as.data.frame(sampleInfo))
  adjustedc <- ComBat(samples, batch=sampleInfo$batch)
  svaob<-adjustedc - min(adjustedc)
  svobj <- svaseq(svaob, mod1, mod0) 
  comsvare <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(comsvare)
}

svacom<- function(svaob,sampleInfo,condition,batch){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + condition, as.data.frame(sampleInfo))
  svobj <- svaseq(as.matrix(svaob), mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv[c(1,2)])) 
  svacomre <- ComBat(cleaned_count, batch=sampleInfo$batch)
  return(svacomre)
}

svAll<- function(svaob,sampleInfo,condition,batch){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + condition, as.data.frame(sampleInfo))
  svobj <- svaseq(svaob, mod1, mod0) 
  svAllre <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(svAllre)
}

com<- function(samples,sampleInfo,batch){
  comre <- ComBat(samples, batch=sampleInfo$batch)
  return(comre)
}

### functions pca

PreparePlotpca<-function(counts,sample_info,var1,var2,values){
  pca_matrix <- counts %>% 
    as.matrix() %>% 
    t()
  sample_pca <- prcomp(pca_matrix)
  tibpca<-as_tibble(pca_matrix, rownames = "sample")
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>%
    mutate(pct = variance/sum(variance)*100) %>% 
    mutate(pct_cum = cumsum(pct))
  pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct), fill="deeppink4") +
    geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")
  pc_scores <- sample_pca$x
  pc_scores <- pc_scores %>% 
    as_tibble(rownames = "sample")
  #wow<-list(sample_pca,pc_eigenvalues)
  p<-sample_pca$x %>% 
    as_tibble(rownames = "sample") %>% 
    full_join(sample_info, by = "sample") %>% 
#    ggplot(aes(x =PC1, y = PC2, col = !!sym(var1), shape = !!sym(var2))) +
    ggplot(aes(x =PC1, y = PC2, col = condition, shape = batch)) +
    scale_colour_manual(values=values) +
    geom_point(size=5) +
    geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
    labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
    ggtitle("PCA of condition and batch") +
    theme_bw()
  p
}


## without manual colour
title="PCA of condition and batch"
PreparePlotpca2<-function(counts,sample_info,var1,var2){
  #title=paste("PCA of ",var1," and ",var2,sep="")
  pca_matrix <- counts %>% 
    as.matrix() %>% 
    t()
  sample_pca <- prcomp(pca_matrix)
  tibpca<-as_tibble(pca_matrix, rownames = "sample")
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>%
    mutate(pct = variance/sum(variance)*100) %>% 
    mutate(pct_cum = cumsum(pct))
  pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct), fill="deeppink4") +
    geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")
  pc_scores <- sample_pca$x
  pc_scores <- pc_scores %>% 
    as_tibble(rownames = "sample")
  #wow<-list(sample_pca,pc_eigenvalues)
  p<-sample_pca$x %>% 
    as_tibble(rownames = "sample") %>% 
    full_join(sample_info, by = "sample") %>% 
    #ggplot(aes(x =PC1, y = PC2, col = !!sym(var1), shape = !!sym(var2))) +
    ggplot(aes(x =PC1, y = PC2, col = condition, shape = batch)) +
    scale_fill_brewer() +
    geom_point(size=5) +
    geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
    labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
    ggtitle(title) +
    theme_bw()
  p
}

##################### runs
samples=read.table(args[1], sep="\t", row.names=1, header=T)
#samples<-read.table("/home/bianca/Desktop/Dockerfiles/preprocess/normalized.txt", sep="\t", row.names = 1, stringsAsFactors = FALSE)
#

samplesinfo=read.table(args[2], sep="\t", row.names=1, header=T)
#samplesinfo<-read.table("/home/bianca/Desktop/Dockerfiles/preprocess/all_sample_info.txt", sep="\t", row.names = 1)
#print(head(samples))
#print(head(samplesinfo))
samples<-read.table("/tmp//normalized.txt", sep="\t", row.names = 1)
samplesinfo<-read.table("/tmp//brca_sampleinfo_pt.txt", sep="\t", row.names = 1)

if(dim(samples)[2]!=dim(samplesinfo)[1]){
  colnames(samplesinfo)=samplesinfo[1,]
  samplesinfo=samplesinfo[-1,]
}
samplesinfo$sample=rownames(samplesinfo)
if(typeof(samples[1,1])=="character"){
    colnames(samples)=samples[1,]
    samples=samples[-1,]
}
#print(head(rownames(samples)))
#samplesr<-rownames(samples)
#print(head(samples))
#print(head(samplesinfo))
#print("make numeric")
samples<-samples[is.finite(rowSums(as.matrix(samples))),]
samplesr<-rownames(samples)
samples<-apply(samples,2,as.numeric)
#print("done apply")
#print(head(samples))
rownames(samples)<-samplesr
#print(head(samples))
#print(head(samplesinfo))
sorder<-samplesinfo$sample
samples=samples[,sorder]
#samples<-samples[,samplesinfo$sample]
#colnames(samples)=samplesinfo$sample
#has_inf <- apply(samples, 1, function(row) any(is.infinite(row)))
#samples <- samples[!has_inf, ]

#print(head(samples))
#print(head(samplesinfo))
print(colnames(samples)==samplesinfo$sample)
#args4<-args[4]
#args5<-args[5]
#samplesinfo$condition<-samplesinfo$args4
#samplesinfo$batch<-samplesinfo$args5
#print(head(samplesinfo))

vec0<-which(colnames(samplesinfo)=="sampleID")
vec1<-which(colnames(samplesinfo)==as.character(args[4]))
vec2<-which(colnames(samplesinfo)==as.character(args[5]))
colnames(samplesinfo)[vec0] <- "sample"
colnames(samplesinfo)[vec1] <- "condition"
colnames(samplesinfo)[vec2] <- "batch"


samplesinfo[] <- lapply(samplesinfo, factor)
print(head(samplesinfo))
#print(samplesinfo$batch)
print("pca of raw")
svg("pca_raw.svg")
PreparePlotpca2(samples,samplesinfo,"condition","batch")
dev.off()


print("cleanining")
if (args[3]=="svacom"){
  print("clean with sva and combat")
  cleaned=svacom(samples,samplesinfo,"condition","batch")
  write.table(cleaned,"cleaned_svacom.txt",sep="\t", quote=F)
  #svg("svacom.svg")
}
if (args[3]=="comsva"){
  print("clean with combat and sva")
  cleaned=comsva(samples,samplesinfo,"condition","batch")
  write.table(cleaned,"cleaned_comsva.txt",sep="\t", quote=F)
  #svg("comsva.svg")
  #title=paste0("PCA of ", as.character(args[4]), " and " ,as.character(args[5]), ".")
  #svg("comsva.svg")
  ##PreparePlotpca2(cleaned_comsva,samplesinfo,args[4],as.factor(args[5]))
  #PreparePlotpca2(cleaned_comsva,samplesinfo,"condition","batch")
  #dev.off()
}
if (args[3]=="com"){
  print("clean with combat")
  cleaned=com(samples,samplesinfo,"batch")
  write.table(cleaned,"cleaned_com.txt",sep="\t", quote=F)
  #svg("com.svg")
  #title=paste0("PCA of ", as.character(args[4]), " and " ,as.character(args[5]), ".")
  #svg("com.svg")
  #PreparePlotpca2(cleaned_comsva,samplesinfo,"condition","batch")
  #PreparePlotpca2(cleaned_com,samplesinfo,as.character(args[4]),as.character(args[5]))
  #dev.off()
}
if (args[3]=="sva"){
  print("clean with sva")
  cleaned=svAll(samples,samplesinfo,args[4],args[5])
  write.table(cleaned,"cleaned_sva.txt",sep="\t", quote=F)
  #title=paste0("PCA of " ,as.character(args[4]), " and " ,as.character(args[5]), ".")
  #svg("sva.svg")
  #PreparePlotpca2(cleaned_comsva,samplesinfo,"condition","batch")
  #PreparePlotpca2(cleaned_sva,samplesinfo,args[4],as.factor(args[5]))
  #dev.off()
}

title=paste0("PCA of ", as.character(args[4]), " and ", as.character(args[5]) ,".")
image<-paste0(args[3],".svg")
svg(image)
PreparePlotpca2(cleaned,samplesinfo,"condition","batch")
dev.off()

if (args[3]=="all"){
  print("clean with sva and combat")
  cleaned1=svacom(samples,samplesinfo,"condition","batch")
  write.table(cleaned1,"cleaned_svacom.txt",sep="\t", quote=F)
  #svg("svacom.svg")
  #title=paste0("PCA of ", as.character(args[4]), " and ", as.character(args[5]) ,".")
  #svg("svacom.svg")
  #PreparePlotpca2(cleaned_svacom,samplesinfo,"condition","batch")
  #dev.off()


  print("clean with combat and sva")
  cleaned2=comsva(samples,samplesinfo,"condition","batch")
  write.table(cleaned2,"cleaned_comsva.txt",sep="\t", quote=F)
  #svg("comsva.svg")
  #title=paste0("PCA of ", as.character(args[4]), " and " ,as.character(args[5]), ".")
  #svg("comsva.svg")
  #PreparePlotpca2(cleaned_comsva,samplesinfo,args[4],as.factor(args[5]))
  #PreparePlotpca2(cleaned_comsva,samplesinfo,"condition","batch")
  #dev.off()

  print("clean with combat")
  cleaned3=com(samples,samplesinfo,"batch")
  write.table(cleaned3,"cleaned_com.txt",sep="\t", quote=F)
  #svg("com.svg")
  #title=paste0("PCA of ", as.character(args[4]), " and " ,as.character(args[5]), ".")
  #svg("com.svg")
  #PreparePlotpca2(cleaned_comsva,samplesinfo,"condition","batch")
  #PreparePlotpca2(cleaned_com,samplesinfo,as.character(args[4]),as.character(args[5]))
  #dev.off()

  print("clean with sva")
  cleaned4=svAll(samples,samplesinfo,args[4],args[5])
  write.table(cleaned4,"cleaned_sva.txt",sep="\t", quote=F)
  #title=paste0("PCA of " ,as.character(args[4]), " and " ,as.character(args[5]), ".")
  #svg("sva.svg")
  #PreparePlotpca2(cleaned_comsva,samplesinfo,"condition","batch")
  ##PreparePlotpca2(cleaned_sva,samplesinfo,args[4],as.factor(args[5]))
  #dev.off()
  
}

title=paste0("PCA of ", as.character(args[4]), " and ", as.character(args[5]) ,".")
svg("svacom.svg")
PreparePlotpca2(cleaned1,samplesinfo,"condition","batch")
dev.off()
svg("comsva.svg")
PreparePlotpca2(cleaned2,samplesinfo,"condition","batch")
dev.off()
svg("com.svg")
PreparePlotpca2(cleaned3,samplesinfo,"condition","batch")
dev.off()
svg("sva.svg")
PreparePlotpca2(cleaned4,samplesinfo,"condition","batch")
dev.off()

if (args[3]=="none"){
  print("no cleaning all is perfect as it is!")
  write.table(samples,"cleaned_samples.txt",sep="\t", quote=F)
}


save.image("cleaned.RData")


library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stats)


#method="pearson"
#cut.off=0.8
#pval_cutoff=0.1
args <- commandArgs(trailingOnly=TRUE)

genenames<-read.table(args[1],sep="\t", row.names=1, header=T)
lncrnsnames<-read.table(args[1],sep="\t", row.names=1, header=T)
genenames$ID_Symbol<-rownames(genenames)
lncrnsnames$ID_Symbol<-rownames(lncrnsnames)

method=args[3]
cut.off=args[4]
pval_cutoff=args[5]


calculate_permutation_correlation <- function(mrna_row, lncrna_row, method = "pearson") {
  lncrna_row_permuted <- sample(lncrna_row)
  correlation_result <- cor.test(mrna_row, lncrna_row_permuted, method = method)
  return(correlation_result$estimate)
}

cal_cor_big <- function(data_1, data_2,  method = "pearson", cut.off, num_permutations = 10, pval_cutoff) {
  corr <- list()
  # Create all combinations of row indices using expand.grid
  combinations <- expand.grid(row_data_1 = 1:nrow(data_1), row_data_2 = 1:nrow(data_2))
  for (k in 1:nrow(combinations)) {
   # mrna<-subset(data_1,select=-c(ID,ID_Symbol,logFC,Pval))
   # lncrna<-subset(data_2,select=-c(ID,ID_Symbol,logFC,Pval))
    mrna<-subset(data_1,select=-c(ID,ID_Symbol))
    lncrna<-subset(data_2,select=-c(ID,ID_Symbol))
    i <- combinations$row_data_1[k]
    j <- combinations$row_data_2[k]
    mrna_row <- as.numeric(mrna[i, ])  # Exclude the first column (e.g., "ID")
    lncrna_row <- as.numeric(lncrna[j, ])  # Exclude the first column (e.g., "ID")
    observed_correlation <- cor.test(mrna_row, lncrna_row, method = method)
    permutation_correlations <- replicate(
      num_permutations,
      calculate_permutation_correlation(mrna_row, lncrna_row, method = method)
    )
    p_value <- sum(permutation_correlations >= abs(observed_correlation$estimate)) / num_permutations
    if (!is.na(observed_correlation$estimate) && !is.na(p_value) && 
        abs(observed_correlation$estimate) > cut.off && p_value < pval_cutoff) {
      #if (!is.na(observed_correlation) && !is.na(p_value) && 
      #    abs(observed_correlation) > cut.off && p_value < pval_cutoff) {
      corr[[length(corr) + 1]] <- c(
        lncRNA = data_2[j, "ID"],
        Gene = data_1[i, "ID"],
        lncRNA_Symbol = data_2[j, "ID_Symbol"],
        Gene_Symbol = data_1[i, "ID_Symbol"],
        Correlation = as.numeric(observed_correlation$estimate)
        #Observed_Correlation = as.numeric(observed_correlation),
       # Pval_Correlation = as.numeric(p_value),
       # logFC_lncRNA = data_2[j, "lncRNA_logFC"],
       # Pval_lncRNA = data_2[j, "lncRNA_Pval"],
       # logFC_gene = data_1[i, "Gene_logFC"],
       # Pval_gene = data_1[i, "Gene_Pval"]
      )
    }
  }
  
  if (length(corr) == 0) {
    return(NULL)
  } else {
    corr <- do.call(rbind, corr)
    #colnames(corr) <- c("lncRNA", "Gene", "lncRNA_Symbol", "Gene_Symbol","Correlation", "Pval_Correlation",
    #                    "logFC_lncRNA", "Pval_lncRNA",
    #                    "logFC_gene", "Pval_gene")
    return(corr)
  }
}

neg_cor_b_big <- function(mrna_data, lncRNA_data, method = "pearson", cut.off = 0.1, num_permutations = 10, pval_cutoff=0.1) {
  corr <- cal_cor_big(mrna_data, lncRNA_data, method, cut.off, num_permutations, pval_cutoff)
  return(corr)
}

cal_cor <- function(data_1, data_2,  method = "pearson", cut.off, num_permutations = 10, pval_cutoff) {
  corr <- list()
  # Create all combinations of row indices using expand.grid
  combinations <- expand.grid(row_data_1 = 1:nrow(data_1), row_data_2 = 1:nrow(data_2))
  for (k in 1:nrow(combinations)) {
    i <- combinations$row_data_1[k]
    j <- combinations$row_data_2[k]
    mrna_row <- as.numeric(data_1[i,-1 ])  # Exclude the first column (e.g., "ID")
    lncrna_row <- as.numeric(data_2[j,-1 ])  # Exclude the first column (e.g., "ID")
    observed_correlation <- cor.test(mrna_row, lncrna_row, method = method)
    permutation_correlations <- replicate(
      num_permutations,
      calculate_permutation_correlation(mrna_row, lncrna_row, method = method)
    )
    p_value <- sum(permutation_correlations >= abs(observed_correlation$estimate)) / num_permutations
    if (!is.na(observed_correlation$estimate) && !is.na(p_value) && 
        abs(observed_correlation$estimate) > cut.off && p_value < pval_cutoff) {
      #if (!is.na(observed_correlation) && !is.na(p_value) && 
      #    abs(observed_correlation) > cut.off && p_value < pval_cutoff) {
      corr[[length(corr) + 1]] <- c(
        lncRNA = data_2[j, "ID_Symbol"],
        Gene = data_1[i, "ID_Symbol"],
        Correlation = as.numeric(observed_correlation$estimate),
        #Observed_Correlation = as.numeric(observed_correlation),
        Pval_Correlation = as.numeric(p_value))
    }
  }
  
  if (length(corr) == 0) {
    return(NULL)
  } else {
    corr <- do.call(rbind, corr)
   # colnames(corr) <- c("lncRNA", "Gene","Correlation", "Pval_Correlation")
    return(corr)
  }
}

neg_cor_b <- function(mrna_data, lncRNA_data, method = "pearson", cut.off = 0.1, num_permutations = 10, pval_cutoff=0.1) {
  corr <- cal_cor(mrna_data, lncRNA_data, method, cut.off, num_permutations, pval_cutoff)
  return(corr)
}



# short gene nemes
matrx<-neg_cor_b(genenames,lncrnsnames,method,cut.off)



dfmat<-as.data.frame(matrx)
head(dfmat)
ftenln<-unique(dfmat$lncRNA)[1:10]
head(ftenln)
df1<-dfmat[dfmat$lncRNA %in% ftenln,]


#df1<-as.data.frame(matrx[1:25,])
df1$Correlation<-as.numeric(df1$Correlation)
print(head(df1))
svg("correlation_between_features.svg")
ggplot(as.data.frame(df1), aes(lncRNA, Gene, fill= Correlation)) + 
  geom_tile() + ggtitle("Correlation of 10 first features and other features") 
dev.off()

result1 <- dfmat %>%
  group_by(lncRNA) %>%
  summarize(targets = paste(unique(Gene), collapse = ","))


targets<-unique(dfmat$Gene)
write.table(targets,"targets.txt", sep="\t", quote = F, row.names = F)
write.table(result1,"correlation_features_and_targets.txt", sep="\t", quote = F, row.names = F)

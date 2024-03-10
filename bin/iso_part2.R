#args1 = path to switchList part1
# args2 = path/to/cpat results
# args3 = path/to/pfam result
# args4 = path/to/signaip result
# args5 = TRUE/FALSE reduce to gene switches
# args6 = alpha (fdr)
# args7 = dif or fold change
# args8 =  run saturn TRUE/FALSE
# args7 = output dexseq
# args8 =  output satur switchlist

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(optparse))

args <- commandArgs(trailingOnly=TRUE)

SwitchList<-readRDS(args[1])

# SwitchList <- isoformSwitchAnalysisPart2(
#   switchAnalyzeRlist        = SwitchList, 
#   n                         = Inf,    # if plotting was enabled, it would only output the top 10 switches
#   removeNoncodinORFs        = TRUE,
#   pathToCPC2resultFile      = system.file("extdata/cpc2_result.txt"         , package = "IsoformSwitchAnalyzeR"),
#   pathToPFAMresultFile      = system.file("extdata/pfam_results.txt"        , package = "IsoformSwitchAnalyzeR"),
#   pathToIUPred2AresultFile  = system.file("extdata/iupred2a_result.txt.gz"  , package = "IsoformSwitchAnalyzeR"),
#   pathToSignalPresultFile   = system.file("extdata/signalP_results.txt"     , package = "IsoformSwitchAnalyzeR"),
#   outputPlots               = FALSE  # keeps the function from outputting the plots from this example
# )

file_info<-file.info(args[4])

if (file_info$size > 0) {
SwitchListDex <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList, 
  codingCutoff =  0.364,
  n                         = Inf,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPATresultFile      = args[2],
  pathToPFAMresultFile      = args[3],
  pathToSignalPresultFile   = args[4],
 # pathToIUPred2AresultFile  = args[5],
  outputPlots               = FALSE ,# keeps the function from outputting the plots from this example
  consequencesToAnalyze = c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified','signal_peptide_identified') 
)}else{
SwitchListDex <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList,
  codingCutoff =  0.364,
  n                         = Inf,    # if plotting was enabled, it would only output the top 10 switches
  removeNoncodinORFs        = TRUE,
  pathToCPATresultFile      = args[2],
  pathToPFAMresultFile      = args[3],
  #pathToSignalPresultFile   = args[4],
 # pathToIUPred2AresultFile  = args[5],
  outputPlots               = FALSE ,# keeps the function from outputting the plots from this example
  consequencesToAnalyze = c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status','domains_identified'))
}

lncList<-subsetSwitchAnalyzeRlist(
    SwitchListDex, 
    SwitchListDex$isoformFeatures$gene_biotype=="lncRNA"
)

library(biomaRt)
values<-gsub("\\..*", "",unique(lncList$isoformFeatures$isoform_id))
print(head(values))
filters="ensembl_transcript_id"
attributes=c("ensembl_transcript_id", "hgnc_symbol","chromosome_name", "start_position", "end_position")

homo.anno = useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")

lncs_cor<-getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
write.table(lncs_cor,"lncrnas.txt", sep="\t", quote=F, row.names=T)
saveRDS(SwitchListDex, args[9])
isoanalysis<-as.data.frame(SwitchListDex$isoformSwitchAnalysis)
write.table(isoanalysis, "isoSwitchAnalysis.txt", sep="\t", quote=F, row.names=T)
orfanalysis<-as.data.frame(SwitchListDex$orfAnalysis)
write.table(orfanalysis, "orfAnalysis.txt", sep="\t", quote=F, row.names=T)
res<-as.data.frame(SwitchListDex$isoformFeatures)
write.table(res, "results.txt", sep="\t", quote=F, row.names=T)

library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stats)
library(tidyverse)


method="pearson"
cut.off=0.8
pval_cutoff=0.1

# Function to calculate correlations for a single permutation
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


# inputs
#all<-read.table("/home/bianca/Desktop/correlation_input_switch_detailed.txt", sep="\t")

#all<-as.data.frame(cbind(SwitchListDex$isoformFeatures$isoform_id,SwitchListDex$isoformFeatures$gene_id,SwitchListDex$isoformFeatures$iso_biotype, SwitchListDex$isoformFeatures$iso_value_1,as.numeric(SwitchListDex$isoformFeatures$iso_value_1 + SwitchListDex$isoformFeatures$iso_stderr_1), as.numeric(SwitchListDex$isoformFeatures$iso_value_1 - SwitchListDex$isoformFeatures$iso_stderr_1), SwitchListDex$isoformFeatures$iso_value_2,as.numeric(SwitchListDex$isoformFeatures$iso_value_2 + SwitchListDex$isoformFeatures$iso_stderr_2), as.numeric(SwitchListDex$isoformFeatures$iso_value_2 - SwitchListDex$isoformFeatures$iso_stderr_2),SwitchListDex$isoformFeatures$iso_log2_fold_change, SwitchListDex$isoformFeatures$iso_q_value ))
#colnames(all)=c("ID","ID_Symbol","biotype","con1","con1","con1","con2","con2","con2","logFC","Pval")
#all<-all[,c("ID","ID_Symbol","logFC","Pval","biotype","con1","con1","con1","con2","con2","con2")]
#all<-subset(all, Pval < 0.05 )

df<-SwitchListDex$isoformRepExpression
write.table(df,"isoforms_expression.txt",sep="\t",quote=F,row.names=F)


difisres<-SwitchListDex$isoformSwitchAnalysis
difisres<-subset(difisres, padj < 0.05)
difisres<-subset(difisres, abs(dIF) > 0.1, select=iso_ref)

write.table(difisres,"de_isoforms.txt",sep="\t",quote=F,row.names=F)

all<-as.data.frame(cbind(SwitchListDex$isoformFeatures$iso_ref,SwitchListDex$isoformFeatures$isoform_id,SwitchListDex$isoformFeatures$gene_id,SwitchListDex$isoformFeatures$iso_biotype))
all<-all[all$V1 %in% difisres$iso_ref,]

isomat<-as.data.frame(SwitchListDex$isoformCountMatrix)

all <- full_join(all, isomat, by = c("V2" = "isoform_id"))

all <- all %>%
  mutate_at(-c(1:4), as.numeric)

colnames(all)[2:3]<-c("ID","ID_Symbol")
#head(all)

genesnames<- subset(all, V4 == "protein_coding", select = -c(1,2,4))
lncrnsnames <- subset(all, V4 == "lncRNA",  select = -c(1,2,4))

genes<- subset(all,V4== "protein_coding",select=-c(1,4))
#print(head(genes))
lncrns <- subset(all, V4== "lncRNA",select=-c(1,4))
#print(head(lncrns))

#all<-cbind(all, all[,6:ncol(all)])

#genesnames<- subset(all, biotype == "protein_coding", select = c(ID_Symbol,  con1, con1, con1, con2, con2, con2))
#lncrnsnames <- subset(all, biotype == "lncRNA", select = c(ID_Symbol,  con1, con1, con1, con2, con2, con2))

#genes<- subset(all, biotype == "protein_coding",select=-biotype)
#lncrns <- subset(all, biotype == "lncRNA",select=-biotype)

# genesnames<-genesnames[1:20,]
# lncrnsnames<-lncrnsnames[1:25,]
# genes<-genes[1:20,]
# lncrns<-lncrns[1:25,]


## detailed analysis

matrx<-neg_cor_b_big(genes,lncrns,method,cut.off)
#head(matrx)
dfmat<-as.data.frame(matrx)
write.table(dfmat,"lncrna_correlation_targets_detailed.txt", sep="\t", quote=F, row.names=F)

# short gene nemes
matrx<-neg_cor_b(genesnames,lncrnsnames,method,cut.off)



dfmat<-as.data.frame(matrx)
head(dfmat)
ftenln<-unique(dfmat$lncRNA)[1:10]
head(ftenln)
df1<-dfmat[dfmat$lncRNA %in% ftenln,]


#df1<-as.data.frame(matrx[1:25,])
df1$Correlation<-as.numeric(df1$Correlation)
print(head(df1))
svg("correlation_lncRNAs_Genes.svg")
ggplot(as.data.frame(df1), aes(lncRNA, Gene, fill= Correlation)) + 
  geom_tile() + ggtitle("Correlation of 10 first lncRNA and Genes") 
dev.off()

result1 <- dfmat %>%
  group_by(lncRNA) %>%
  summarize(targets = paste(unique(Gene), collapse = ","))

#library(dplyr)
#print(head(lncs_cor))
#print(head(result1))
result2<-left_join(result1, lncs_cor, by=c("lncRNA"="hgnc_symbol"))
result2<-result2[,c(1,3,4,5,6,2)]

#head(result1)
write.table(result2,"lncrna_correlation_targets.txt", sep="\t", quote = F, row.names = F)






if(as.logical(args[6])){
  print("now saturn")
  SwitchListSat <- isoformSwitchTestSatuRn(
    ### Core arguments
    switchAnalyzeRlist =SwitchList,
    alpha = args[7],
    dIFcutoff = args[8],
    
    ### Advanced arguments
    reduceToSwitchingGenes = as.logical(args[5]),
    reduceFurtherToGenesWithConsequencePotential = TRUE,
    onlySigIsoforms = FALSE,
    keepIsoformInAllConditions = TRUE,
    diagplots = TRUE,
    showProgress = TRUE,
    quiet = FALSE
  )
  saveRDS(SwitchListSat, args[10])
  }else if(!as.logical(args[6])){
    file_path <- args[10]
    file_conn <- file(file_path, "w")
    text <- "You didn't run satuRn, wuabalabadab dab"
    writeLines(text, file_conn)
    close(file_conn)
}

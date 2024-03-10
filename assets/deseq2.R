#Rscript /r/deseq2.R /deseq2_files/ /sampleTable.txt T "~0 + condition" mkc dmso defeatures.txt
suppressMessages(library("DESeq2"))
suppressMessages(library(sva)) #for combat-seq
suppressMessages(library("stringr"))
suppressMessages(library(dplyr))
library(tidyverse)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)

# Get input file and output file names from command-line arguments
args <- commandArgs(trailingOnly=TRUE)
correctBatch<-as.logical(args[3])
#correctBatch<-as.logical(T)
single_matrix=as.logical(args[7])

#path<-args[1] ## path of files
#path="/home/bianca/Desktop/testdeseq2/"
#sampleTable<-read.table(args[2], sep=",", row.names=1, header=T) 
## path of sampleTable
sampleTable<-read.table(args[2], sep=",", header=T) 
#sampleTable<-read.table("/home/bianca/Desktop/testdeseq2/sampleTable.csv", sep=",")
#num<-colnames(sampleTable)

sampleTable <- sampleTable[!duplicated(sampleTable$sampleID), ]
#colnames(sampleTable)<-num

sampleTable$groups = dense_rank(sampleTable$condition)
if(correctBatch){
  sampleTable$batches<-dense_rank(sampleTable$batch)
}

rownames(sampleTable)<-sampleTable$sampleID

print("Importing formula")
#form<-"~0+ condition"
form<-args[4] ### formula of design matrix
print(form)
formula<-as.formula(form)
design <- model.matrix(formula, data=sampleTable)
#design<-model.matrix(~sampleTable$condition)


if(!single_matrix){
    path<-args[1]
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable =  sampleTable,
                                       directory = path,
                                       design=formula)
                                       #design = ~ condition
                                       #design = ~ batch + condition, #alternative design for 2 conditions, e.g. batch, male/female etc)
}

if(single_matrix){
path=read.table(args[1], sep="\t", row.names=1, header=T)
path <- round(path)
path<-path - min(path)
#print(head(path))
#print(head(sampleTable))
ddsHTSeq <-DESeqDataSetFromMatrix(colData =  sampleTable,
                                       countData = path,
                                       design=formula)
}

if(correctBatch){
  counts=counts(ddsHTSeq)
  corrected_data = ComBat_seq(counts = as.matrix(counts), batch = sampleTable$batches, group = sampleTable$groups)
  corrected_data2 <- apply (corrected_data, c (1, 2), function (x) {
    (as.integer(x))
  })
  assay(ddsHTSeq)= corrected_data2
}

keep <- rowSums(counts(ddsHTSeq)==0) <= round(0.3*ncol(counts(ddsHTSeq)))
ddsHTSeq <- ddsHTSeq[keep, ]

ddsHTSeq <- estimateSizeFactors(ddsHTSeq)



performFirstBlock <- function() {
  dds <- DESeq(ddsHTSeq)
  dds <- vst(dds)
  print("Successfully completed the first block.")
  return(dds)
}

# Function to perform the second block
performSecondBlock <- function() {
  dds <- estimateDispersionsGeneEst(ddsHTSeq)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  print("Successfully completed the second block.")
  return(dds)
}

# Try the first block
aha <- tryCatch({
  performFirstBlock()
}, error = function(e1) {
  cat("Error in the first block:", conditionMessage(e1), "\n")
  
  # Try the second block
  tryCatch({
    performSecondBlock()
  }, error = function(e2) {
    cat("Error in the second block:", conditionMessage(e2), "\n")
    stop("Both attempts failed.")
  })
})

# Check if the result is a DESeqDataSet
if (inherits(aha, "DESeqDataSet")) {
  # Continue with the rest of your analysis using 'result'
} else {
  cat("Downstream analysis cannot be performed due to errors in both attempts.\n")
}

#con1<-"dmso"
#con2<-"mkc"
con1<-args[5]
con2<-args[6]


print("outputing results")  
res= results(aha, contrast=c("condition",con1, con2))

# # Write results to output file
# write.table(res, file=output_file, sep="\t", quote=FALSE)


write.table(res, file="defeatures.txt", sep="\t", quote=FALSE)
de<-as.data.frame(res)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > 0.1 & de$padj < 0.05] <-  "UP"
de$diffexpressed[de$log2FoldChange < -0.1 & de$padj < 0.05] <-  "DOWN"
de<-de[1:100,]
bi<-ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + geom_point() +
  scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4")) +
  geom_text(aes(label = rownames(de)), vjust = -0.5, hjust = 0.5, size = 3) +
  geom_hline(yintercept=-log10(0.05), col="red") +
  labs (y = "-log10(P value)", x = "Log fold change") +
  ggtitle("Significantly differentially expressed with logFC above 0.1")  + theme_minimal()
svg("volcano_plot.svg")
plot(bi)
dev.off()
y<-counts(ddsHTSeq) 
de_features_exp<-y[rownames(y) %in% rownames(de),]
test<-scale(as.data.frame(de_features_exp))
test<-test[1:20,]
col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5
signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(condition=sampleTable[,"condition"])
row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
svg("Heatmap.svg")
hm<-Heatmap(test, name="expression",
            col=col_fun,
            cluster_rows=T,
            cluster_columns = T,
            row_names_side = "left",
            show_row_names = T,
            show_column_names = F,
            show_row_dend = FALSE,
            show_column_dend=FALSE,
            row_names_gp=gpar(cex=fontsize+0.1),
            row_names_max_width = unit(5, "cm"),
            clustering_distance_rows ="euclidean",
            clustering_method_rows = "ward.D",
            clustering_distance_columns =  "euclidean",
            clustering_method_columns = "ward.D",
            row_dend_width = unit(10, "mm"),
            left_annotation=row,
            bottom_annotation = col,
            width = NULL)

draw(hm)
dev.off()

print("DONE")



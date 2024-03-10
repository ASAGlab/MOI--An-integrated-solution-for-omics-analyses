
# Rscript prepare_df.R --rna T --proteins T --integrated T --path "/home/bianca/biotrans/tmp/"
library(optparse)
library(stringr)
library(dplyr)
library(biomaRt)
################################################################################################################################


########################################### make options python wise ###############################################################
option_list = list(
  make_option(
    "--rna",
    action = "store",
    default = TRUE,
    type = 'logical',
    help ='mRNA'
  ),
  make_option(
    "--mirna",
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'miRNA'
  ),
  make_option(
    "--proteins",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Proteins'
  ),
  make_option(
    "--lipids",
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Lipids'
  ),
  make_option(
    "--isoforms",
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Isoforms'
  ),
  make_option(
    "--integrated",
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Isoforms'
  ),
  make_option(
    "--integratedafterlipids",
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Isoforms'
  ),
  make_option(
    "--path",
    action = "store",
    default = "/home/bianca/biotranscompare/tmp/",
    type = 'character',
    help = 'Path to inputs'
  ),
  make_option(
    "--alg_genes",
    action = "store",
    default = "rp",
    type = 'character',
    help = 'Path to samples info'
  ),
  make_option(
    "--alg_mirna",
    action = "store",
    default = "rp",
    type = 'character',
    help = 'Path to samples info'
  ),
  make_option(
    "--alg_proteins",
    action = "store",
    default = "rp",
    type = 'character',
    help = 'Path to samples info'
  ),
  make_option(
    "--pval",
    action = "store",
    default = 0.05,
    type = 'double',
    help = 'Path to samples info'
  )
  
)

opt <- parse_args(OptionParser(option_list=option_list))
pval=opt$pval
path=opt$path
features<-list()
featheat<-data.frame(cbind("Gene", "Category"))
colnames(featheat)<-c("Gene", "Category")
#path<- "/home/bianca/biotranscompare/tmp"

keepfeat<-function(alg,genes){
  if(alg=="edger"){
    genes<-genes[,c(1,4)]
    colnames(genes)<-c("logFC","pval")
    genes<-genes[genes$pval < pval,]
    genes<-rownames(genes)
    return(genes)
  }else if(alg=="deseq2"){
    genes<-genes[,c(2,6)]
    colnames(genes)<-c("logFC","pval")
    genes<-genes[genes$pval < pval,]
    genes<-rownames(genes)
    return(genes)
  }else if(alg=="rp"){
    genes<-genes[,c(3,4)]
    colnames(genes)<-c("logFC","pval")
    genes<-genes[genes$pval < pval,]
    genes<-rownames(genes)
    return(genes)
  }
  
}


#path<- "/home/bianca/biotranscompare/tmp"
if(opt$rna){
  #rnapath <- file.path(opt$path,"genes_defeatures")
  rnapath <- file.path(opt$path,"genes_defeatures.txt")
  pattern <- "*.txt"
  file_rna <- list.files(rnapath, pattern = pattern, full.names = TRUE)
  print(file_rna[1])
  #rna<-read.table(file_rna[1], sep="\t", row.names = 1, header=T)
  rna<-read.table(rnapath, sep="\t", row.names = 1, header=T)
  rna2<-keepfeat(opt$alg_genes,rna)
  features$rna<-as.data.frame(rna2)
  colnames(features$rna)<-"Genes"
  rna2<-as.data.frame(cbind(rna2,rep("Genes",length(rna2))))
  colnames(rna2)<-c("Gene","Category")
  featheat<-rbind(featheat, rna2)
}
if(opt$mirna){
  mirnapath <- file.path(opt$path,"targets.txt")
  #mirnapath <- file.path(opt$path,"mirna_defeatures")
  #pattern <- "*.txt"
  #file_mirna <- list.files(mirnapath, pattern = pattern, full.names = TRUE)
  
  #mirna<-read.table(file_mirna[1], sep="\t", row.names = 1, header=T)
  mirna<-read.table(mirnapath, sep="\t", row.names = 1, header=T)
  #mirna2<-keepfeat(opt$alg_mirna,mirna)
  #features$mirna<-as.data.frame(mirna2)
  #homo.anno<-useEnsembl(biomart = "ensembl",
  #                    dataset = "hsapiens_gene_ensembl")
  #values=mirna2
  #if(length(values)!=0){
  #mirna.df=getBM(c("hgnc_symbol","mirbase_accession"),"mirbase_accession",values,homo.anno)
  #colnames(mirna.df)=c("hgnc_symbol","V2")
  #genesmirna = unique(mirna.df$hgnc_symbol)
  #}
  genesmirna<-mirna
  print(head(genesmirna))
  features$mirna<-as.data.frame(genesmirna)
  colnames(features$mirna)<-"miRNA"
  mirna2<-as.data.frame(cbind(genesmirna,rep("miRNA",length(genesmirna))))
  colnames(mirna2)<-c("Gene","Category")
  featheat<-rbind(featheat, mirna2)
}
if(opt$proteins){
  print("proteins")
  proteinspath <- file.path(opt$path,"proteins_defeatures.txt")
  #proteinspath <- file.path(opt$path,"proteins_defeatures")
  # pattern <- "*.txt"
  # file_protein <- list.files(proteinspath, pattern = pattern, full.names = TRUE)
  # print(file_protein)
  # proteins<-read.table(file_protein[1], sep="\t", row.names = 1, header=T)
  proteins<-read.table(proteinspath, sep="\t", row.names = 1, header=T)
  proteins2<-keepfeat(opt$alg_proteins,proteins)
  features$proteins<-as.data.frame(proteins2)
  colnames(features$proteins)<-"Proteins"
  proteins2<-as.data.frame(cbind(proteins2,rep("Proteins",length(proteins2))))
  colnames(proteins2)<-c("Gene","Category")
  featheat<-rbind(featheat, proteins2)
}
if(opt$integrated){
  integratedpath <- file.path(opt$path,"integrated.txt")
  #proteinspath <- file.path(opt$path,"proteins_defeatures")
  # pattern <- "*.txt"
  # file_protein <- list.files(proteinspath, pattern = pattern, full.names = TRUE)
  # print(file_protein)
  # proteins<-read.table(file_protein[1], sep="\t", row.names = 1, header=T)
  genes1<-read.table(integratedpath, sep="\t", row.names = 1, header=T)
  tryCatch({
    # Attempt to retrieve orthology for the module
    if(any(genes1$lipids==TRUE)){
      integrated_lipids<-as.data.frame(genes1$var[genes1$lipids==TRUE]) 
      integrated_lipids$Category<-rep("Integrated_Lipids",length(integrated_lipids))
      colnames(integrated_lipids)[1]<-"Gene"
      #write.table(integrated_lipids,"integrated_Lipids.txt", sep="\t", row.names = F)
    }
  }, error = function(e) {
    # Handle errors
    cat("NO integrated lipids", conditionMessage(e), "\n")
    
  })
  features$integrated<-as.data.frame(genes1$var[genes1$lipids!=TRUE]) 
  integrated2<-as.data.frame(cbind(features$integrated,rep("Integrated",length(features$integrated))))
  colnames(integrated2)<-c("Gene","Category")
    if (opt$integratedafterlipids){
    integratedLpath <- file.path(opt$path,"integrated_lipids.txt")
    #proteinspath <- file.path(opt$path,"proteins_defeatures")
    # pattern <- "*.txt"
    # file_protein <- list.files(proteinspath, pattern = pattern, full.names = TRUE)
    # print(file_protein)
    # proteins<-read.table(file_protein[1], sep="\t", row.names = 1, header=T)
    genes2<-read.table(integratedLpath, sep="\t",  header=T)
    print(head(genes2))
    integrated2<-rbind(integrated2, genes2)
  }else if(!opt$integratedafterlipids){
    integrated2<-rbind(integrated2, integrated_lipids)
    }
  featheat<-rbind(featheat, integrated2)
}
if(opt$isoforms){
  #isoformspath <- file.path(opt$path,"isoforms")
  isoformspath <- file.path(opt$path,"isoforms.txt")
  #pattern <- "*.txt"
  #file_isoforms <- list.files(isoformspath, pattern = pattern, full.names = TRUE)
  #print(file_isoforms[1])
  #isoforms<-read.table(file_isoforms[1], sep="\t", row.names = 1, header=T)
  isoforms<-read.table(isoformspath, sep="\t", row.names = 1, header=T)
  isoforms<-rownames(isoforms)
  ### MANY TODOS!!!!!!!!!!!!!!!!!!!!!
  }



#name2<-paste(path,"genes_across_omics.txt",sep="/")
write.table(featheat, "genes_across_omics.txt",row.names = F, quote = F, sep="\t")


max_rows <- max(sapply(features, function(x) length(unlist(x))))


#featuresnomi<-features[features !="miRNA"]
features2 <- lapply(features, function(x) {
  if (length(unlist(x)) < max_rows) {
    c(unlist(x), rep(NA, max_rows - length(unlist(x))))
  } else {
    unlist(x)
  }
})
# Combine the columns row-wise
result <- do.call(cbind, features2)
namecom<-paste(path,"bio_comp.txt",sep="/")
write.table(result,"bio_comp.txt", sep="\t",quote=F, na="", row.names = F)


#write.table(result,namecom, sep="\t",quote=F, na="", row.names = F)

# 
# for (i in names(features)){
#   #nam<-featheat[[i]][1,2]
#   name<-paste(colnames(features[[i]]),"bio.txt",sep="_")
#   print(name)
#   name2<-paste(path,name,sep="/")
#   df<-as.data.frame(features[[i]])
#   write.table(df, name2,row.names = F, quote = F, sep="\t")
# }


args<- commandArgs(trailingOnly=TRUE)
gene_symbols<-read.table(args[1], sep="\t", header=T, row.names = 1)
gene_symbols<-gene_symbols[which(gene_symbols$P.Value < 0.05),]
gene_symbols<-rownames(gene_symbols)
head(gene_symbols)
# index <- which(grepl("C\\d+:\\d+", gene_symbols))[1]
# gene_s<-gene_symbols[index:length(gene_symbols)]

gene_symbols_clean <- gsub("C(\\d+):(\\d+) ([A-Z]+)", "\\3(\\1:\\2)", gene_symbols)
print(head(gene_symbols_clean))
gene_symbols<-c(gene_symbols_clean)
library(MetaboAnalystR)


tmp.vec <-gene_symbols
mSet<-InitDataObjects("conc", "msetora", FALSE)

mSet<-Setup.MapData(mSet, tmp.vec);



mSet<-Setup.MapData(mSet, gene_symbols_clean);

mSet<-CrossReferencing(mSet, "name")

mSet<-CreateMappingResultTable(mSet)



mSet <- GetCandidateList(mSet);


# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, F);

# Select metabolite set library, refer to 
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0);

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)


mSet<-PlotORA(mSet, "lipids_pathway_enrichment", "bar", "svg", 72, width=NA)

#args1 = path to salmon output directory
# args2 = path/to/deignMatrix
# args3 = path/to/gtf
# args4 = path/to/fasta
# args4 = path/to/output
# args6 = path/to/outputSwitchList

library(IsoformSwitchAnalyzeR)
library(optparse)

args <- commandArgs(trailingOnly=TRUE)
#args="/home/bianca/R/x86_64-pc-linux-gnu-library/4.2/IsoformSwitchAnalyzeR/extdata"

# salmonFileDataFrame <-prepareSalmonFileDataFrame(
#   ### Core arguments
#  args[1],
#   
#   ### Advanced arguments
#   pattern='',
#   invertPattern=FALSE,
#   ignore.case=FALSE,
#   quiet = FALSE
# )

myDesign<-read.table(args[2],header=T,sep=",")
vec0<-which(colnames(myDesign)=="sampleID") 
vec1<-which(colnames(myDesign)=="condition") 
vec2<-which(colnames(myDesign)=="batch") 


if(length(vec0)!=0){colnames(myDesign)[vec0]<-"sampleID"}
if(length(vec0)==0){myDesign$sampleID<-rownames(myDesign)}

myDesign<-myDesign[,c(vec0,vec1,vec2)]

salmonQuant<-importIsoformExpression(
  ### Core arguments
  args[1],
)

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = args[3],
  isoformNtFasta       =args[4],
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE,
  #removeNonConvensionalChr = TRUE
)

aSwitchList<-preFilter( aSwitchList )
#aSwitchList<-isoformSwitchTestDEXSeq(aSwitchList, alpha = args[5],dIFcutoff = args[6],reduceFurtherToGenesWithConsequencePotential = FALSE)
#aSwitchList<-extractSequence(aSwitchList, alpha=args[5],dIFcutoff=args[6], pathToOutput=args[7] )
alpha=as.numeric(args[5])
dif<-as.numeric(args[6])
print(alpha)
print(dif)
aSwitchList<-isoformSwitchTestDEXSeq(aSwitchList, alpha = alpha,dIFcutoff = dif)
aSwitchList<-extractSequence(aSwitchList, alpha=alpha,dIFcutoff=dif, pathToOutput=args[7] )

#aSwitchList <- isoformSwitchAnalysisPart1(
#  switchAnalyzeRlist   = aSwitchList,
#  pathToOutput =args[5],
#  outputSequences      = TRUE, # change to TRUE whan analyzing your own data 
#  prepareForWebServers = FALSE  # change to TRUE if you will use webservers for external sequence analysis
#)

saveRDS(aSwitchList, args[8])

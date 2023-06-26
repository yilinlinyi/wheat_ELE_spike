library(CAGEr)
library(BSgenome.wheat.IWGSC)
library(dplyr)
args <- commandArgs(T)
sample <- args[1]
ce <- CAGEexp( genomeName     = "BSgenome.wheat.IWGSC",
inputFiles     = paste(sample,".ctss",sep="") ,
inputFilesType = "ctss" ,
sampleLabels   = sample)
getCTSS(ce)
CTSStagCountSE(ce)
CTSScoordinatesGR(ce)
CTSStagCountDF(ce)
###distribution and obtain α
#png(paste(sample,"cumulative_distribution.png",sep="_"))
#plotReverseCumulatives(ce, fitInRange = c(5, 20000), onePlot = TRUE)
#dev.off()
##normalize
a1 <- c("cs_cage_4tissue_merge","cs_cage_embryo_merge","cs_cage_root_merge","cs_cage_seedling_merge","cs_cage_spikelet_I_merge")
a2 <- c(1.07,1.07,1.13,1.08,1.18)
sa <- data.frame(V1=a1,V2=a2)
alp <- filter(sa,V1==sample)[,2]
normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5,20000),alpha = alp ,T=10000000)
####cluster call
clusterCTSS( object = ce
           , threshold = 0.5
           , thresholdIsTpm = TRUE
           , nrPassThreshold = 1
           , method = "distclu"
           , maxDist = 50
           , removeSingletons = TRUE
           , keepSingletonsAbove = 5)
test_cluster <- tagClusters(ce, sample =sample)
write.table(test_cluster,paste(sample,".bed",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
dat <- tagClustersGR(ce, sample,returnInterquantileWidth = TRUE,qLow=0.1,qUp=0.9)
write.table(dat,paste(sample,"_TC_width.txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)

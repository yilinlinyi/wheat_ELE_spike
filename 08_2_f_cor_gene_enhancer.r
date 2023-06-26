library(data.table)
library(dplyr)
args <- commandArgs(T)
file <- args[1]
rpm <- fread("enhancer_promoter_CTC_rep_rpm.txt")
ab <- read.table(file)
ab$V1 <- as.character(ab$V1)
ab$V3 <- as.character(ab$V3)

pair_cor <- function(x){
d <- as.data.frame(x[c(1,3)])
names(d)<-"id"
d1 <- left_join(d,rpm,by="id")
row.names(d1) <- d1$id
d1 <- d1[,-1]
result=cor.test(t(d1)[,1],t(d1)[,2])
cor=result$estimate
p=result$p.value
as.vector(c(p,cor))
}

ab <- filter(ab,V1!=V3)
re <- apply(ab,1,pair_cor)
result <- cbind(ab,as.data.frame(t(re)))
write.table(result,paste(file,"_cage_cor.txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)

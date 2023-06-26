library(data.table)
library(dplyr)
RP_cor_50 <- function(x){
d <- fread(paste(x,".sp.hic_2M.cage_cor0.5.genelist",sep=""),head=F)%>%
 unique()
rp  <- fread(paste(x,"_sp_gene_RP.txt",sep=""),head=F)
n <- quantile(rp$V2,0.5)
rp1 <- filter(rp,V2>n)
over <- semi_join(rp1,d,by="V1")[,1]%>%
unique()
write.table(over,paste(x,".sp.2M_cage_cor0.5_RP_top50.genelist",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
}

RP_cor_50("seedling")
RP_cor_50("embryo")
RP_cor_50("root")
RP_cor_50("spikelet_I")

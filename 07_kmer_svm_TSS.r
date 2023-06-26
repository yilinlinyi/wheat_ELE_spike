library(ggplot2)
library(dplyr)
library(data.table)
library(pROC)
library(randomForest)
library(nnet)
library(multiROC)


d <- fread("en_gene_8020.out")
roc(d$V3,d$V2,plot=TRUE,legacy.axes=T,percent=TRUE,xlab="False positive percentage",ylab="True positive percentage",lwd=4,col="#4885b8",print.auc=TRUE)
###trianed
roc.info=roc(d$V3,d$V2,egacy.axes=TRUE)
auc=round(roc.info$auc[1],4)*100
roc.df <- data.frame(tpp=roc.info$sensitivities*100,
                     fpp=(1-roc.info$specificities)*100,
                     thresholds=roc.info$thresholds)

ggplot(roc.df,aes(fpp,tpp))+geom_line(size=1.5,color="#4885b8")+theme_bw()+
  annotate("text", x = 75, y = 50, label = paste("AUC\n",auc,"%",sep=""),size=6,fontface=2,color="#4e4846")+
  theme(axis.title.x = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20,color="black"),
        legend.position = "none")+
  labs(x="False positive percentage",y="True positive percentage")
ggsave("slop500_svm_enhancer_promoter_auc_trained.png",width=3.54,height=3.28)




###simi predicted 
pre1 <- fread("en_test.out")
pre1$type <- "enhancer"
pre2 <- fread("gene_test.out")
pre2$type <- "gene"
pre <- rbind(pre1,pre2)
pre <- mutate(pre,true=ifelse(type=="enhancer",1,-1))
roc(pre$true,pre$V2,plot=TRUE,legacy.axes=T,percent=TRUE,xlab="False positive percentage",ylab="True positive percentage",lwd=4,col="#4885b8",print.auc=TRUE)
roc.info=roc(pre$true,pre$V2,egacy.axes=TRUE)
auc=round(roc.info$auc[1],4)*100
roc.df <- data.frame(tpp=roc.info$sensitivities*100,
                     fpp=(1-roc.info$specificities)*100,
                     thresholds=roc.info$thresholds)

ggplot(roc.df,aes(fpp,tpp))+geom_line(size=1.5,color="#4885b8")+theme_bw()+
  annotate("text", x = 75, y = 50, label = paste("AUC\n",auc,"%",sep=""),size=6,fontface=2,color="#4e4846")+
  theme(axis.title.x = element_text(size = 20,color="black"),
        axis.text = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 20,color="black"),
        legend.position = "none")+
  labs(x="False positive percentage",y="True positive percentage")
ggsave("slop500_svm_enhancer_promoter_auc_prediction.png",width=3.54,height=3.28)
ggsave("slop500_svm_enhancer_promoter_auc_prediction.pdf",width=3.54,height=3.28)



####not predicted right
pre <- pre[,c(1,2,4)]
names(pre)[3] <-"real"
d <- d[,c(1,2,3)]
names(d)[3] <-"real"
dat <- rbind(d,pre)
dat$real <- as.factor(dat$real)
ggplot(dat,aes(real,V2))+geom_boxplot()

roc(dat$real,dat$V2,plot=TRUE,legacy.axes=T,percent=TRUE,xlab="False positive percentage",ylab="True positive percentage",lwd=4,col="#4885b8",print.auc=TRUE)
###find optimal cutoff

library(cutpointr)
cp <- cutpointr(dat, V2, real, 
                method = maximize_metric, metric = sum_sens_spec)
summary(cp)
dat  <- mutate(dat,pre=ifelse(V2<(-1.6381),"-1","1"))
wr <- filter(dat,real!=pre)
wr_e  <- filter(wr,real==1) ####1159
write.table(wr_e[,1],"enhancer_not_correct.list.txt",quote=F,sep="\t",col.names=F,row.names=F)


###not_corr and simi
d <- fread("enhancer_not_correct.list.txt",head=F)
al <- fread("/Users/xieyilin/Library/CloudStorage/OneDrive-个人/SIPPE_备份/working/CAGE/multi-tissue_202105/TE_TSS_figure/p0_all_CTC/h_p_vs_e/consensus/01_all_blast_all/v4_all_enhancer/promoter_enhancer_simi.txt",head=F)
dat <- left_join(d,al,by=c("V1"="V2"))%>%
  na.omit()
dat <- dat[,1:3]
names(dat)  <- c("e","p","ratio")
pos <- fread("/Users/xieyilin/Library/CloudStorage/OneDrive-个人/SIPPE_备份/working/CAGE/multi-tissue_202105/TE_TSS_figure/p1/support_file/cs_cage_4tissue_merge.dominant.merge.bed")
pos <- pos[,c(1,2,4)]
dat <- left_join(dat,pos,by=c("e"="V4"))%>%
  left_join(pos,by=c("p"="V4"))
names(dat)[4:7] <- c("chre","pose","chrp","posp")
dat <- mutate(dat,chr=ifelse(chre==chrp,"same","diff"))
dat1 <- filter(dat,chr=="same"&chre!="chrUn")
dat1 <- mutate(dat1,dis=abs(pose-posp))
ggplot(dat1,aes(dis))+stat_ecdf(geom="step")+
scale_x_log10(breaks=c(0,10,100,1000,10000,100000,100000000))
dat2 <- filter(dat1,dis<=2000000&ratio>0.1)
length(unique(dat2$e))

simi <- fread("/Users/xieyilin/Library/CloudStorage/OneDrive-个人/mac_work-2022/CAGE/p0/h/func/01/01_cor_expr/pr_en_pair_blast_length_cage_cor.txt")
simi <- mutate(simi,type=ifelse(V3>=30,"high","low"))
simi <- simi[,c(1,5)]
simi1 <- filter(simi,type=="high")
simi2 <- anti_join(simi,simi1,by="V1")
simi <- rbind(simi1,simi2)%>%
  unique()
names(simi)[2]  <- "vs_gene"
dat <- left_join(pre,simi,by="V1")
dat <- na.omit(dat)
ggplot(dat,aes(type_p,fill=vs_gene))+geom_bar(stat="count",position="fill")

write.table(dat,"not_classify_enhancer_blast_promoter.txt",quote=F,sep="\t",col.names = T,row.names = F)


roc(as.factor(pre$type),pre$V2,plot=TRUE,legacy.axes=T,percent=TRUE,xlab="False positive percentage",ylab="True positive percentage",lwd=4,col="#4885b8",print.auc=TRUE)
ggplot(pre,aes(type,V2))+geom_boxplot()




##wheat_ABD_mcscan
##01 prepare blast and gff files
makeblastdb -in IWGSC_v1.1_HC_20170706_pep.longsest.geneid.fa   -dbtype prot
blastp -query  IWGSC_v1.1_HC_20170706_pep.longsest.geneid.fa   -db  IWGSC_v1.1_HC_20170706_pep.longsest.geneid.fa   -out aa_bb_dd.blast   -outfmt  6  -num_threads 16   -num_alignments 10
cat IWGSC_v1.1_HC_20170706.merge.gene.bed |awk '$1~/A/' |sed 's/chr/aa/'|sed 's/A//' |awk 'BEGIN{OFS="\t"}{print $1,$4,$2,$3}' > a.tmp
cat IWGSC_v1.1_HC_20170706.merge.gene.bed |awk '$1~/B/' |sed 's/chr/bb/'|sed 's/B//' |awk 'BEGIN{OFS="\t"}{print $1,$4,$2,$3}'> b.tmp
cat IWGSC_v1.1_HC_20170706.merge.gene.bed |awk '$1~/D/' |sed 's/chr/dd/'|sed 's/D//' |awk 'BEGIN{OFS="\t"}{print $1,$4,$2,$3}'> d.tmp
cat a.tmp  b.tmp  d.tmp   > aa_bb_dd.gff
MCScanX aa_bb_dd -m 5   -b 2 -s 2
grep -v "#" aa_bb_dd.collinearity |sed 's/ //g'|sed 's/://'|awk '$2~/A/&&$3~/B/'|awk '{print $2"\t"$1"\n"$3"\t"$1}' > AB.txt
grep -v "#" aa_bb_dd.collinearity |sed 's/ //g'|sed 's/://'|awk '$2~/A/&&$3~/D/'|awk '{print $2"\t"$1"\n"$3"\t"$1}' > AD.txt
grep -v "#" aa_bb_dd.collinearity |sed 's/ //g'|sed 's/://'|awk '$2~/B/&&$3~/D/'|awk '{print $2"\t"$1"\n"$3"\t"$1}' > BD.txt
 Rscript f_bed.r

##f_bed.r
library(dplyr)
library(data.table)
bed1 <- fread("/public/home/xieyilin/yilin/genome/iwgsc_cs/IWGSC_v1.1_HC_20170706.merge.gene.bed")
bed2 <- fread("/public/home/xieyilin/yilin/genome/wew_AABB/WEW_v2.0.gene.bed")
bed <- rbind(bed1,bed2)
i <- "cs_A_vs_AB_A"
g <-fread(paste(i,".txt",sep=""),head=F)
re <- left_join(bed,g,by=c("V4"="V1"))
re <- na.omit(re)
write.table(re,paste(i,".bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)

###############################
##02 synteny block
for i in `ls ~/mcscan_genome/default/*.bed`
do
name=`basename $i .bed`
awk 'BEGIN{OFS="\t"}{$6="+";print $0}' $i > $name.arrange.bed
done
ln -s ~/TE_TSS/p2_evolution/f_synteny/01_pos/enhancer.dominant.merge.bed  ./

ln -s AB.arrange.bed BA.arrange.bed
ln -s AD.arrange.bed DA.arrange.bed
ln -s BD.arrange.bed DB.arrange.bed


for F in A B D
do
for R in A B D
do
if [ $F != $R ];then
###tag the last gene in a block
awk -v aa=$F '$1~aa'  ${F}${R}.arrange.bed|sed 's/-/\t/'|sort -k7,7 -k8,8nr |sort -k7,7 -u |sed 's/\t/-/7'|awk 'BEGIN{OFS="\t"}{print $0"\t""last"}'  |sort -k1,1 -k2,2n > ${F}${R}_${F}_last.bed
cut -f 7  ${F}${R}_${F}_last.bed > ${F}${R}_${F}_last.id
awk -v aa=$F '$1~aa'  ${F}${R}.arrange.bed|grep -F -w -v -f  ${F}${R}_${F}_last.id|awk 'BEGIN{OFS="\t"}{print $0"\t""top"}'|cat - ${F}${R}_${F}_last.bed |sort -k1,1 -k2,2n > ${F}${R}_${F}_type.bed

###find the eRNA downdream closest gene corresponding block tag(if tag is the top,remove it )
awk -v aa=$F '$1~aa'  enhancer.dominant.merge.bed |sort -k1,1 -k2,2n |awk 'BEGIN{OFS="\t"}{$6="+";print $0}' |bedtools closest -a  - -b  ${F}${R}.arrange.bed  -D  "a" -iu |awk 'BEGIN{OFS="\t"}{print $0,"2"}'  > ${F}${R}_${F}_down.bed
cut -f 4,13 ${F}${R}_${F}_down.bed |awk '$2!="."'|sed 's/-/\t/'|awk '$3!="0"' |awk 'BEGIN{OFS="\t"}{print $1,$2"-"$3-1,$2"-"$3}'|sort -k1,1  > ${F}${R}_${F}_with_block.bed
cut -f1   ${F}${R}_${F}_with_block.bed |sort -u >  ${F}${R}_${F}_with_block.id

###if the filtered  eRNA  have another block which,find it by searching the eRNA uptream tag which is not the last one
cut -f 4,13 ${F}${R}_${F}_down.bed |sed 's/-/\t/'|awk '$3=="0"' |grep -F -w -v -f  ${F}${R}_${F}_with_block.id - |cut -f 1 |grep -F -w -f - enhancer.dominant.merge.bed |sort -k1,1 -k2,2n |awk 'BEGIN{OFS="\t"}{$6="+";print $0}'|bedtools closest -a  - -b ${F}${R}_${F}_type.bed -D  "a" -id |grep  top|cut -f 4,13 |sed 's/-/\t/'|awk 'BEGIN{OFS="\t"}{print $1,$2"-"$3,$2"-"$3+1}'|sort -k1,1  > ${F}${R}_${F}_with_block.sup.bed
cat ${F}${R}_${F}_with_block.bed  ${F}${R}_${F}_with_block.sup.bed  > ${F}${R}_${F}_all_block.bed
Rscript f_pair.r $F $R
Rscript f_pair_info.r $F $R
awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1}'  ${F}_vs_${R}.synteny.bed |sort -k1,1 -k2,2n > ${F}_vs_${R}.synteny.tmp.bed
bash ~/yilin/genome/iwgsc_cs/f_split_v2.sh  ${F}_vs_${R}.synteny.tmp.bed
awk -v aa=$F -v  bb=$R 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,aa"_vs_"bb}' ${F}_vs_${R}.synteny.tmp.bed.split.bed  |sort -k1,1 -k2,2n > ${F}_vs_${R}.synteny.split.bed
fi
done
done


###03 syntenic_block_blast
cat all_peak_DHS.bed  cs_cage_4tissue_merge_unique_best.site_frq.bed  |sort -k1,1 -k2,2n|cut -f 1-3   > cs_all_peak_and_cage.bed

for F in A B D
do
for R in A B D
do
if [ $F != $R ];then
{
name=${F}_vs_${R}
for id in `cat  ${F}_vs_${R}.synteny.2.split.bed|cut -f 1 |sort -u `
do
grep "${id}" cs_cage_4tissue_merge.dominant.bed  |bedtools slop -l 10 -r 200 -g 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai  -i  - |bedtools getfasta -fi 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -bed -  -s|sed 's/(+)//'|sed 's/(-)//'  > ${name}/cs_${F}.1.fa
grep "${id}" ${name}.synteny.2.split.bed|cut -f 5-7|bedtools getfasta -fi 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta  -bed - > ${name}/${R}.1.fa
blastn -query  ${name}/cs_${F}.1.fa -subject ${name}/${R}.1.fa -task  blastn  -outfmt 6  |awk '$11<=0.01'|sed 's/:/\t/g'|sed 's/-/\t/g'|awk -v aa=$id  'BEGIN{OFS="\t"}{if($13<$14)print $1,$2+$11,$2+$12,$4,$5+$13,$5+$14,"+",aa;else print $1,$2+$11,$2+$12,$4,$5+$14,$5+$13,"-",aa}'|cut -f 4-8 >>  ${F}_vs_${R}_homology.bed
done

###arrange blast result
awk '$6=="+"'  cs_cage_4tissue_merge.dominant.bed |cut -f 4 |grep -F -w -f - ${name}_homology.bed |sort -u |sort -k1,1 -k2,2n |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,".",$4}'|bedtools merge -i - -s -d 50 -c 4,6 -o distinct |tr ',' '\t'|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,".",$NF}' > ${name}_plus.tmp
awk '$6=="-"'  cs_cage_4tissue_merge.dominant.bed |cut -f 4 |grep -F -w -f -  ${name}_homology.bed |sort -u |sort -k1,1 -k2,2n |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,".",$4}'|bedtools merge -i - -s -d 50 -c 4,6 -o distinct |tr ',' '\t'|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$(NF-1),".",$NF}' > ${name}_minus.tmp
cat ${name}_plus.tmp ${name}_minus.tmp |sort -k1,1 -k2,2n >  ${name}_homology.merge.bed
bedtools intersect -a ${name}_homology.merge.bed  -b cs_all_peak_and_cage.bed -wa -sorted|sort -u |sort -k1,1 -k2,2n > ${name}_homology.merge.with_peak_or_cage.bed
} &
fi
done
done


###stat
cut -f 1  ${F}_vs_${R}.synteny.2.split.bed|sort -u  | awk '{print $1"\t""with_synteny"}' > ${F}_vs_${R}_with_synteny.txt
cut -f 4  ${name}_homology.merge.bed|sort -u |awk '{print $1"\t""with_homology"}' > ${F}_vs_${R}_with_homology.txt
cut -f 4 ${name}_homology.merge.with_peak_or_cage.bed|sort -u |awk '{print $1"\t""with_enhancer"}' > ${F}_vs_${R}_with_enhancer.txt




###CAGE correlation
gene=IWGSC_v1.1_HC_20170706.split.gene.bed
bedtools slop -r 2000000 -l 2000000  -g 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai  -i cs_enhancer.bed |bedtools intersect -a - -b  $gene  -wao -sorted|cut -f 4,10|sort -u > enhancer_gene_pair_2M.txt 
###cage RPM

cat cs_enhancer.bed  cs_promoter.bed |sort -k1,1 -k2,2n > en_pr.bed 

for i in `cat id.txt`
do
n=`grep -w "${i}" stat_count.txt |cut -f 2 `
bedtools intersect -a en_pr.bed -b ~/data_process/03_arrange/${i}_unique_best.site_frq.bed  -wao -sorted |awk 'BEGIN{OFS="\t"}{if($10==".")$10=0;print $4,$10}'|bedtools groupby -i - -g 1  -c 2 -o sum |awk -v aa=$n -v bb=$i 'BEGIN{print "id""\t"bb }{print $1"\t"$2*1000000/aa}' > en_pr_$i.rpm.txt &
done
wait
paste `ls en_pr_*rpm.txt`|awk '{printf $1"\t";for(i=2;i<=NF;i+=2)printf $i"\t"}{printf "\n"}' > enhancer_promoter_CTC_rep_rpm.txt


split -n l/20  enhancer_gene_TSS_pair_by_2M.txt  -d -a 3 pair
for i in `ls pair*`
do
Rscript f_cor_enhancer_gene.r $i &
done
wait
echo "done"
cat *_cage_cor_test.txt > enhancer_gene_TSS_pair_2M.txt
sort -k2,2 -k1,1 -k5,5nr  enhancer_gene_TSS_pair_2M.txt  |sort -k2,2 -k1,1 -u > enhancer_gene_TSS_pair_2M_uniq.txt
awk '$5>0.5'  enhancer_gene_TSS_pair_2M_uniq.txt  > enhancer_gene_TSS_pair_2M_filter0.5.txt



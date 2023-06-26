###CAGE correlation
gene=IWGSC_v1.1_HC_20170706.split.gene.bed
bedtools slop -r 2000000 -l 2000000  -g 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai  -i cs_enhancer.bed |bedtools intersect -a - -b  $gene  -wao -sorted|cut -f 4,10|sort -u > enhancer_gene_pair_2M.txt 

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


###RP(distance)

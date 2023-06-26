for i in `ls *ctss`
do
id=`basename $i .ctss`
Rscript 01_wheat_cager.r $id &
done
wait
echo "done"

for i in `ls *ctss`
do
id=`basename $i .ctss`
 paste  $id.bed ${id}_TC_width.txt |cut -f 2,3,4,5,8,6,21|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6,$7;printf "%.2f\n",$5}'|paste - -|sed '1d'|sort -k7,7nr|awk 'BEGIN{OFS="\t"}{if($6<=10)print $1,$2,$3,"sharp_"NR"_"$7,$5,$4;else print $1,$2,$3,"broad_"NR"_"$7,$5,$4}'|awk '{printf ("%s\t%d\t%d\t%s\t%9d\t%s\n", $1,$2,$3,$4,$5,$6)}'|sort -k1,1 -k2,2n > ${id}_all.bed
done

cat *all.bed |awk '$NF=="+"'|sort -k1,1 -k2,2n |bedtools merge -i - -s -d 20|awk 'BEGIN{OFS="\t"}{print $0,".",".","+"}'    > plus.bed
cat *all.bed |awk '$NF=="-"'|sort -k1,1 -k2,2n |bedtools merge -i - -s -d 20 |awk 'BEGIN{OFS="\t"}{print $0,".",".","-"}'  > minus.bed
cat  plus.bed minus.bed  |sort -k1,1 -k2,2n > cs_cage_4sample_all_merge.bed
for i in seedling root embryo spikelet_I ; do bedtools intersect -a cs_cage_4sample_all_merge.bed  -b ../../04_rep_merge/cs_cage_${i}_merge_unique_best.site_frq.bed   -wao  -sorted -s |awk 'BEGIN{OFS="\t"}{if($10==".")$10=0;print $1,$2,$3,$4,$5,$6,$10}'|bedtools groupby -i - -g 1,2,3,4,5,6 -c  7  -o sum  >  cs_cage_4sample_all_merge_${i}.bed  & done

wait
 paste `ls  cs_cage_4sample_all_merge_*.bed`|cut -f 1-6,7,14,21,28|awk 'BEGIN{OFS="\t"}{print "merge_tmp_"NR,$0}' >  merge.1.txt

 awk '$8>=10&&$9<=3&&$10<=3&&$11<=3' merge.1.txt > embryo_unique.bed &
 awk '$8<=3&&$9>=10&&$10<=3&&$11<=3' merge.1.txt > root_unique.bed &
 awk '$8<=3&&$9<=3&&$10>=10&&$11<=3' merge.1.txt > seedling_unique.bed &
 awk '$8<=3&&$9<=3&&$10<=3&&$11>=10' merge.1.txt > spikelet_I_unique.bed &
wait
cat  *unique.bed|cut -f 1 |grep -F -w -v -f -  merge.1.txt|cut -f 2-7 > merge.not_singleton.txt

for i in `seq 1 3 `; do cut -f 2-7  embryo_unique.bed |bedtools intersect -a - -b  cs_cage_embryo_rep${i}_unique_best.site_frq.bed  -wao -sorted -s |awk 'BEGIN{OFS="\t"}{if($10==".")$10=0;print $1,$2,$3,$4,$5,$6,$10}'|bedtools groupby -i - -g 1,2,3,4,5,6 -c  7  -o sum >  embryo_unique_rep${i}.bed  & done
for i in `seq 1 2 `;do for sam in root spikelet_I ;  do cut -f 2-7  ${sam}_unique.bed |bedtools intersect -a - -b  cs_cage_${sam}_rep${i}_unique_best.site_frq.bed  -wao -sorted -s |awk 'BEGIN{OFS="\t"}{if($10==".")$10=0;print $1,$2,$3,$4,$5,$6,$10}'|bedtools groupby -i - -g 1,2,3,4,5,6 -c  7  -o sum >  ${sam}_unique_rep${i}.bed  & done  ;  done
for i in `seq 1 4 `; do cut -f 2-7  seedling_unique.bed |bedtools intersect -a - -b  cs_cage_seedling_rep${i}_unique_best.site_frq.bed  -wao -sorted -s |awk 'BEGIN{OFS="\t"}{if($10==".")$10=0;print $1,$2,$3,$4,$5,$6,$10}'|bedtools groupby -i - -g 1,2,3,4,5,6 -c  7  -o sum >  seedling_unique_rep${i}.bed  & done
wait

paste `ls seedling_unique_rep*`|cut -f 1-6,7,14,21,28|awk 'BEGIN{OFS="\t"}{if($7>=5)$7=1;else $7=0; if($8>=5)$8=1;else $8=0;if($9>=5)$9=1;else $9=0;if($10>=5)$10=1;else $10=0 ;print $0}'|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$8+$9+$10}'|awk '$7>=3' >  seedling_unique_filter.bed
paste embryo_unique_rep1.bed  embryo_unique_rep2.bed  embryo_unique_rep3.bed |cut -f 1-6,7,14,21|awk 'BEGIN{OFS="\t"}{if($7>=5)$7=1;else $7=0; if($8>=5)$8=1;else $8=0;if($9>=5)$9=1;else $9=0; print $0}'|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7+$8+$9}'|awk '$7>=2' > embryo_unique_filter.bed
paste root_unique_rep1.bed root_unique_rep2.bed |cut -f 1-6,7,14|awk '$7>=5&&$8>=5' > root_unique_filter.bed
paste spikelet_I_unique_rep1.bed spikelet_I_unique_rep2.bed  |cut -f 1-6,7,14|awk '$7>=5&&$8>=5' > spikelet_I_unique_filter.bed

cat   *unique_filter.bed merge.not_singleton.txt|cut -f 1-6 |sort -k1,1 -k2,2n > cs_cage_4sample_all_merge_filtter.bed

awk '{print $0"\t""id_"NR}'
bedtools intersect -a cs_cage_4sample_all_merge_filtter.bed  -b ../../04_rep_merge/cs_cage_4tissue_merge_unique_best.site_frq.bed  -wao -sorted -s   |bedtools groupby -i -  -g 1-6 -c 10 -o sum|awk 'BEGIN{OFS="\t"}{if(($3-$2)>10)print $1,$2,$3,"broad_"NR"_"$7*1000000/123747891,$5,$6;else print $1,$2,$3,"sharp_"NR"_"$7*1000000/123747891,$5,$6}'  > cs_cage_4tissue_merge_cluster.bed
 bedtools intersect -a cs_cage_4tissue_merge_cluster.bed   -b ../../04_rep_merge/cs_cage_4tissue_merge_unique_best.site_frq.bed  -wao -sorted -s   |sort -k4,4 -k10,10nr|sort -k4,4 -u |awk 'BEGIN{OFS="\t"}{print $1,$8,$9,$4,$2"-"$3,$6}'|sort -k1,1 -k2,2n  > cs_cage_4tissue_merge.dominant.bed &
bedtools intersect -a cs_cage_4tissue_merge_cluster.bed   -b ../../04_rep_merge/cs_cage_4tissue_merge_unique_best.site_frq.bed  -wao -sorted -s   |sort -k4,4 -k10,10nr|sort -k4,4 -u |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9,$6}'|sort -k1,1 -k2,2n  > cs_cage_4tissue_merge_all.bed  &

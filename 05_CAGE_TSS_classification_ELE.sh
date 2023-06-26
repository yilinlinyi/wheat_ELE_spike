##statistic intergenic noncoding with or without histone modification peaks, generate a 0/1 matrix
awk '$5=="noncoding"'  all_ctc_info.1.txt |grep intergenic |cut -f 1 > intergenic_noncoding.txt
grep -F -w -f intergenic_noncoding.txt  ~/data_process/05_cager/split/anno/cs_cage_4tissue_merge_all.bed  > intergenic_noncoding_all.bed

for i in `cat peaks/id.txt`
do
bedtools intersect -a intergenic_noncoding_all.bed  -b peaks/macs_CS_${i}_rep1_peaks.bed  -wao -sorted |awk 'BEGIN{OFS="\t"}{if($7==".")$7=0;else $7=1;print $0}'|cut -f 1-7 |sort  -u -k1,1 -k2,2n  >  intergenic_noncoding_${i}_01.bed
done
paste `ls intergenic_noncoding_*01.bed`|awk '{printf $4"\t"$7"\t";for(i=14;i<=NF;i+=7)printf $i"\t"}{printf "\n"}' |cat head.txt  - |sed 's/\t$//' > intergenic_noncoding.01.txt

cut -f 1,3,4,6,7,9,10,11,12| intergenic_noncoding.01.txt |awk '{print $1"\t"$2+$3+$4+$5+$6+$7+$8+$9}'|sed '1d'|awk 'BEGIN{OFS="\t"}{if($2>0)print $1,"active";else print $1,"non"}' >  intergenic_noncoding_active_histone_01.txt

grep active intergenic_noncoding_active_histone_01|cut -f 1 > cs_ELE_TSS.id

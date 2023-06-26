jaspar=JASPAR2020_CORE_non-redundant_pfms_meme.txt
###scan +-3kb region around TSS 
enhancer_TSS=cs_ELE.bed

##enhancer
bedtools slop -r 3000 -l 3000 -g 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai  -i $enhancer_TSS |bedtools getfasta -fi  161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -bed - -s -name >  cs_ELE_TSS_3k.fasta
fimo --o  all_enhancer_fimo --max-stored-scores 1000000  $jaspar    cs_ELE_TSS_3k.fasta
##promoter
cut -f 1  cs_cage_coding_classification.gene.txt |grep -F -w -f - cs_cage_4tissue_merge.dominant.bed  |bedtools slop -r 3000 -l 3000 -g 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai  -i - |bedtools getfasta -fi  161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -bed - -s -name >  cs_promoter_3k.fasta
fimo --o  all_promoter_fimo --max-stored-scores 1000000  $jaspar  cs_promoter_3k.fasta

##background
bedtools shuffle -i cs_promoter_3k.fasta  -g ~/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai  | head -n 2000 |sort -k1,1 -k2,2n > background.bed
#bedtools getfasta -fi 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -bed bak.bed    -name+ -s |tr ':' '\t'|cut -f 1  >   background.fa
fimo --o  background_fimo  --max-stored-scores 1000000 $jaspar   background.fa


##full length RLG_famc7.3_LTR
bedtools intersect -a  RLG_famc7.3_fl.LTR.bed  -b RLG_famc7.3_enhancer.bed  -wa -sorted |sort -u > fl_7.3_with_enhancer.bed
bedtools subtract -a RLG_famc7.3_fl.LTR.bed  -b  cs_cage_4tissue_merge_cluster.bed  -A -sorted |shuf |head -n 115 >  fl_7.3_wo_enhancer.bed
awk 'BEGIN{OFS="\t"}{$4="with_"$4;print $0}'  fl_7.3_with_enhancer.bed|bedtools getfasta -fi ~/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -bed - -name   >  fl_7.3_with_enhancer.fa
awk 'BEGIN{OFS="\t"}{$4="wo_"$4;print $0}'  fl_7.3_wo_enhancer.bed |bedtools getfasta -fi ~/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -bed - -name   >  fl_7.3_wo_enhancer.fa
cat fl_7.3_with_enhancer.fa fl_7.3_wo_enhancer.fa >  fl_7.3_w_wo_enhancer.fa
fimo --o fimo_of_7.3  --max-stored-scores 1000000  $jaspar    fl_7.3_w_wo_enhancer.fa










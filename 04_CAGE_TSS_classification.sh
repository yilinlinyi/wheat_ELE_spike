##goal:assign each coding CTC a CDS, a transcripts ,a protein,a gene model
#01 used the RNA-seq of 10 wheat tissues to assemble transcripts 
###tissue:CS_Callus;CS_Embryo;CS_Flag;CS_Leaf;CS_Root;CS_Seedling;CS_Sheath;CS_Spikelet_I;CS_Spikelet_II;CS_Stem
ref=161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
##1.1
for i in `ls CS_*.Hisat2.q20.sorted.bam`
do
stringtie  $i   -p 60   --rf  -o  `basename $i .bam`.gff3
done
stringtie  --merge -o cs_10tissue_merged.gtf  *gff3

##1.2
gtf_genome_to_cdna_fasta.pl  cs_10tissue_merged.gtf  $ref > transcripts.fasta  
gtf_to_alignment_gff3.pl cs_10tissue_merged.gtf  > transcripts.gff3  
#wait
TransDecoder.LongOrfs -S -t transcripts.fasta
dir=transcripts.fasta.transdecoder_dir
##blastp:
blastp -query  $dir/longest_orfs.pep  \
-db uniprot-reviewed_yes+taxonomy_58023.fasta -max_target_seqs 1 \
-outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6 
###pfam search
hmmsearch  --cpu 10 --domtblout pfam.domtblout ~/software/interproscan-5.28-67.0/data/pfam/31.0/pfam_a.hmm $dir/longest_orfs.pep  
###predict
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3  transcripts.gff3  transcripts.fasta >  transcripts.fasta.transdecoder.genome.gff3

seqtk seq transcripts.fasta.transdecoder.pep |paste - - |grep complete |tr ' ' '\t'|cut -f 1,9|tr '\t' '\n'  > transcripts.fasta.transdecoder_compelete.pep
grep '>' transcripts.fasta.transdecoder_compelete.pep|sed 's/>//' >  complete_CDS.id
sed 's/\./\t/2' complete_CDS.id |cut -f 1 |awk '{print "ID="$1";"}' |cat - complete_CDS.id  > complete_CDS_gene.id
grep -F -f complete_CDS_gene.id transcripts.fasta.transdecoder.genome.gff3  > complete_transcript.gff3 

###unique CDS by 5UTR
awk '$3=="five_prime_UTR"' complete_transcript.gff3 | sed 's/;/\t/' | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$10,"five_prime_UTR",$7}' | sed 's/Parent=//' | sort -k1,1 -k2,2n > complete_5UTR.bed
awk '$3=="mRNA"' complete_transcript.gff3|sed 's/;/\t/'|awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$9,"mRNA",$7}'|sed 's/ID=//' |sort -k1,1 -k2,2n > complete_mRNA.bed
awk '$3=="CDS"' complete_transcript.gff3 | sed 's/;/\t/' | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$10,"CDS",$7}' | sed 's/Parent=//' | sort -k1,1 -k2,2n > complete_CDS.bed
awk '$3=="gene"' complete_transcript.gff3  | sed 's/;/\t/' |awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$9,"gene",$7}'|sed 's/ID=//'|sort -k1,1 -k2,2n  > complete_transcript_gene.bed

sort -k4,4 -k1,1 -k2,2n complete_5UTR.bed | bedtools groupby -i - -g 4 -c 1,2,3,6 -o distinct,min,max,distinct|awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,".",$5}' |sort -k1,1 -k2,2n > complete_5UTR_merge.bed
cat  transcripts.fasta.transdecoder_compelete.pep|paste - - |awk '{print $1"\t"length($2)}'|sed 's/>//'|sort -k1b,1 >  transcripts.fasta.transdecoder_compelete_pep_length.txt


#02.CTC classification

###HC_5UTR
bedtools intersect -a cs_cage_4tissue_merge.dominant.bed  -b HC_5UTR_merge_pep_length.bed  -wao -s -sorted |awk '$7!="."'|sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""5_UTR"}' > 1.list
###LC_5UTR
cut -f 1 1.list  |grep -F -w -v -f -   cs_cage_4tissue_merge.dominant.bed |bedtools intersect -a - -b LC_5UTR_merge_pep_length.bed  -wao -s -sorted |awk '$7!="."'|sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""5_UTR"}' > 2.list

###assembly_with_gene
 cat 1.list 2.list|cut -f 1 |grep -F -w -v -f -   cs_cage_4tissue_merge.dominant.bed |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene.bed -wao -s -sorted |awk '$7!="."' |grep -v NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""5_UTR_lnRNA"}' > 3.list
###assembly_without_gene
cat 1.list 2.list 3.list|cut -f 1 |grep -F -w -v -f -   cs_cage_4tissue_merge.dominant.bed |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene.bed -wao -s -sorted |awk '$7!="."' |grep  NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""5_UTR_lnRNA"}' > 4.list
##tssup500
cat 1.list 2.list 3.list 4.list  |cut -f 1 |grep -F -w -v  -f  - cs_cage_4tissue_merge.dominant.bed  |bedtools intersect -a - -b cs_HC_tssup500bp.pep.bed  -wao -sorted -s |awk '$7!="."'|sort -k4,4 -k11,11nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""tssup500"}' >  5.list
cat 1.list 2.list 3.list 4.list  5.list|cut -f 1 |grep -F -w -v  -f  - cs_cage_4tissue_merge.dominant.bed  |bedtools intersect -a - -b cs_LC_tssup500bp.pep.bed  -wao -sorted -s |awk '$7!="."'|sort -k4,4 -k11,11nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""tssup500"}' >  6.list
cat 1.list 2.list 3.list 4.list  5.list 6.list |cut -f 1 |grep -F -w -v  -f  - cs_cage_4tissue_merge.dominant.bed  |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene_tssup500.bed  -wao -sorted -s |awk '$7!="."'|grep -v NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""tssup500_lnRNA"}' > 7.list
cat 1.list 2.list 3.list 4.list  5.list 6.list 7.list  |cut -f 1 |grep -F -w -v  -f  - cs_cage_4tissue_merge.dominant.bed  |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene_tssup500.bed  -wao -sorted -s |awk '$7!="."'|grep  NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""tssup500_lnRNA"}' > 8.list
###rescue
cat 1.list 2.list 3.list  4.list 5.list 6.list 7.list 8.list  | cut -f 1 |grep --color=auto -F -w -v -f -  cs_cage_4tissue_merge_all.bed |sort -k1,1 -k2,2n > cs_cage_torescue.bed
bedtools intersect -a cs_cage_torescue.bed -b HC_5UTR_merge_pep_length.bed  -wao -s -sorted |awk '$7!="."'|sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""5_UTR"}' > 1.list.1
cut -f 1 1.list.1 |grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b LC_5UTR_merge_pep_length.bed  -wao -s -sorted |awk '$7!="."'|sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""5_UTR"}' > 2.list.1
cat 1.list.1 2.list.1 |cut -f 1|grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene.bed -wao -s -sorted |awk '$7!="."' |grep -v NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""5_UTR_lnRNA"}' > 3.list.1

cat 1.list.1 2.list.1 3.list.1 |cut -f 1|grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene.bed -wao -s -sorted |awk '$7!="."' |grep  NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""5_UTR_lnRNA"}' > 4.list.1

cat 1.list.1 2.list.1 3.list.1  4.list.1 |cut -f 1|grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b cs_HC_tssup500bp.pep.bed  -wao -sorted -s |awk '$7!="."'|sort -k4,4 -k1,11nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""tssup500"}' >  5.list.1
cat 1.list.1 2.list.1 3.list.1  4.list.1 5.list.1 |cut -f 1|grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b cs_LC_tssup500bp.pep.bed  -wao -sorted -s |awk '$7!="."'|sort -k4,4 -k1,11nr |sort -k4,4 -u |cut -f 4,10|awk '{print $1"\t"$2"\t""tssup500"}' >  6.list.1
cat 1.list.1 2.list.1 3.list.1  4.list.1 5.list.1 6.list.1 |cut  -f 1|grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b complete_5UTR_merge_pep_length.gene_tssup500.bed  -wao -sorted -s |awk '$7!="."'|grep -v NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""tssup500_lnRNA"}' > 7.list.1
cat 1.list.1 2.list.1 3.list.1  4.list.1 5.list.1 6.list.1 7.list.1 |cut -f 1|grep -F -w -v -f -   cs_cage_torescue.bed |bedtools intersect -a - -b  complete_5UTR_merge_pep_length.gene_tssup500.bed  -wao -sorted -s |awk '$7!="."'|grep  NA |sort -k4,4 -k13,13nr |sort -k4,4 -u |cut -f 4,10,14|awk '{print $1"\t"$2":"$3"\t""tssup500_lnRNA"}' > 8.list.1


####stat
cat [1-8].list*  > cs_cage_coding_classification.txt
grep -v MS cs_cage_coding_classification.txt |sed 's/\./\t/2'|cut -f 1,2,4 > cs_cage_coding_classification.gene.1.txt
grep MS cs_cage_coding_classification.txt |awk '$2!~/:NA/'|sed 's/:/\t/'|cut -f 1,3,4  > cs_cage_coding_classification.gene.2.txt
grep MS cs_cage_coding_classification.txt |awk '$2~/:NA/'|sed 's/:/\t/'|cut -f 1,2,4 > cs_cage_coding_classification.gene.3.txt
cat cs_cage_coding_classification.gene.[1-3].txt > cs_cage_coding_classification.gene.txt

grep -v MS cs_cage_coding_classification.txt  >  cs_cage_coding_classification.cds.1.txt
grep  MS cs_cage_coding_classification.txt|sed 's/:/\t/'|cut -f 1,2,4 >  cs_cage_coding_classification.cds.2.txt
cat cs_cage_coding_classification.cds.[1-2].txt > cs_cage_coding_classification.cds.txt

for i in `seq 1 8 `; do n=`cat $i.list $i.list.1 |wc -l` ;echo -e  $i"\t"$n >> stat.txt ;  done
cut -f 1  cs_cage_coding_classification.txt|grep -F -w -v -f  -  cs_cage_4tissue_merge_all.bed > cs_cage_nocoding.bed
cut  -f 1 cs_cage_coding_classification.txt|grep -F -w -f  -  cs_cage_4tissue_merge_all.bed >  cs_cage_coding.bed
Rscript f_bar.r
 bedtools slop -r 100 -l 100 -i cs_cage_4tissue_merge.dominant.bed -g ~/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.fai > cs_cage_4tissue_merge.dominant.slop200bp.bed

###rpm
###sample rpm
for i in embryo seedling root spikelet_I
do
n=`awk '{sum+=$4}END{print sum}' ~/yilin/cage/data_process/04_rep_merge/cs_cage_${i}_merge_unique_best.site_frq.bed`
bedtools intersect -a cs_cage_4tissue_merge_all.bed  -b  ~/yilin/cage/data_process/04_rep_merge/cs_cage_${i}_merge_unique_best.site_frq.bed  -wao -sorted -s |awk 'BEGIN{OFS="\t"}{if($10==".")$10=0;print $0}' |bedtools groupby -i - -g 1,2,3,4,5,6 -c 10 -o sum  |awk -v aa=$n 'BEGIN{OFS="\t"}{$7=$7*1000000/aa;print $0}' > cs_cage_4tissue_merge_${i}_rpm.txt &
done
paste `ls cs_cage_*rpm.txt`|cut -f 1-7,14,21,28|cat ~/yilin/cage/data_process/05_cager/permissive/anno/head.txt  -  > cs_cage_sample_rpm.bed


##rpm 01
 sed '1d '  cs_cage_sample_rpm.bed  |awk 'BEGIN{OFS="\t"}{if($7>=0.5)$7="1";else $7=0;if($8>=0.5)$8="1";else $8=0;if($9>=0.5)$9="1";else $9=0;if($10>=0.5)$10="1";else $10=0;print $0}'|cat head.txt  - > cs_cage_sample_rpm01.bed



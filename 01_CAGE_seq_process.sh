ref=161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta
cd $dir/01_cleandata
id=$1 
cutadapt   -a "AGATCGGAAGAGC" -a  "AAAAAAAAAAAA"  -m 20 -o  $id.rmadapt.fq  $dir/00_rawdata/$id.raw.R1.fq.gz  
cd $dir/02_mapping/
bowtie2 -p 20 -x $ref -U  $dir/01_cleandata/$id.rmadapt.fq | samtools view -@ 20 -Sb -h - | samtools sort -@ 20  -o  $id.bowtie.sort.raw.bam -
cd $dir/01_cleandata
cutadapt -u 9 -a  "AAAAAAAAAAAA"  -q 20 -m 20 -o $id.clean.fq $id.rmadapt.fq  
fastqc $id.rmadapt.fq 
fastqc $id.clean.fq 
cd $dir/02_mapping/
bowtie2 -p 20 -x  $ref  -U $dir/01_cleandata/$id.clean.fq  | samtools view -@ 20 -Sb -h - | samtools sort -@ 20    -o  $id.bowtie.sort.bam -
wait
samtools view -@ 20 -Sbh -F 4 -q 1   $id.bowtie.sort.bam   > $id.bowtie.filter1.bam
#unique
samtools view -h  $id.bowtie.filter1.bam  |grep -v "XS:i" |samtools view -Sb -h >  $id.bowtie.unique.bam  &
#best
samtools view -h  $id.bowtie.filter1.bam|grep -E "XS:i|@" |sed 's/AS:i:/AS:i:\t/'|sed 's/XS:i:/XS:i:\t/'|awk '$13!=$15||$1~/^@/'  |sed s'/AS:i:\t/AS:i:/'|sed 's/XS:i:\t/XS:i:/' |samtools view -Sb -h > $id.bowtie.best.bam &
#mutiple best
samtools view -h  $id.bowtie.filter1.bam|grep -E "XS:i|@" |sed 's/AS:i:/AS:i:\t/'|sed 's/XS:i:/XS:i:\t/'|awk '$13==$15||$1~/^@/'  |sed s'/AS:i:\t/AS:i:/'|sed 's/XS:i:\t/XS:i:/' |samtools view -Sb -h > $id.bowtie.multiple_best.bam &
##3mis tmp
samtools view -h $id.bowtie.sort.raw.bam |grep -E  "NM:i:[0-3]\b|@" |samtools view -@ 30 -Sbh >  $id.bowtie.sort.less_3mis.bam &
wait
##3mis
samtools view -@ 30 $id.bowtie.sort.less_3mis.bam |cut -f 1 |sort -k1b,1 >  $id.bias.read.id

rm_rRNA(){
i=$1
bamToBed -i $id.bowtie.$i.bam >$id.bowtie.$i.bed
cut -f 4 $id.bowtie.$i.bed|sort -u|  sed 's/\/1//' > $id.bowtie.$i.id
/public/home/xieyilin/miniconda2/bbtools/lib/filterbyname.sh  in=$dir/01_cleandata/$id.clean.fq   out=${id}.$i.mapped.fq   names=$id.bowtie.$i.id    include=t
sort_path=/public/home/xieyilin/yilin/genome/iwgsc_cs/non_nucl_seq/sortme_database
/public/home/xieyilin/miniconda2/bin/sortmerna --ref \
${sort_path}/rRNA_databases/silva-bac-16s-id90.fasta,${sort_path}/index/silva-bac-16s-id90-db:\
${sort_path}/rRNA_databases/cs_rRNA_chlo_mitro.fasta,${sort_path}/index/cs_rRNA_chlo_mitro-db:\
${sort_path}/rRNA_databases/silva-bac-23s-id98.fasta,${sort_path}/index/silva-bac-23s-id98-db:\
${sort_path}/rRNA_databases/silva-arc-16s-id95.fasta,${sort_path}/index/silva-arc-16s-id95-db:\
${sort_path}/rRNA_databases/silva-arc-23s-id98.fasta,${sort_path}/index/silva-arc-23s-id98-db:\
${sort_path}/rRNA_databases/silva-euk-18s-id95.fasta,${sort_path}/index/silva-euk-18s-id95-db:\
${sort_path}/rRNA_databases/silva-euk-28s-id98.fasta,${sort_path}/index/silva-euk-28s-id98-db:\
${sort_path}/rRNA_databases/rfam-5s-database-id98.fasta,${sort_path}/index/rfam-5s-database-id98-db:\
${sort_path}/rRNA_databases/rfam-5.8s-database-id98.fasta,${sort_path}/index/rfam-5.8s-database-id98-db \
--reads  ${id}.$i.mapped.fq  --num_alignments 1 \
--fastx --aligned ${id}.$i.mapped_rRNA --other ${id}.$i.mapped_non_rRNA --log -a 8 -m 20000  -v

cat ${id}.$i.mapped_non_rRNA.fq|paste - - - -|cut -d ' ' -f 1 |sort -k1,1  > ${id}.$i.mapped_non_rRNA.id
sed 's/^@//' ${id}.$i.mapped_non_rRNA.id|sort --parallel=30  -k1b,1 |grep -F -v -f  $id.bias.read.id -   |grep  -F  -w -f -  $id.bowtie.$i.bed  > ${id}_$i.sortme.bed
awk 'BEGIN{OFS="\t"}{if($6=="+")print $1,$2,$2+1,$6;else print $1,$3-1,$3,$6}'  ${id}_$i.sortme.bed  |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",".",$4}'|sort -k1,1 -k2,2n  -k6,6 > ${id}_$i.sortme.site.bed
bedtools groupby -i ${id}_$i.sortme.site.bed  -g 1,2,3,6  -c 6 -o count |awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,".",$4}'> ${id}_$i.sortme.site.frq.bed
bedtools subtract -a ${id}_$i.sortme.site.frq.bed  -b /public/home/xieyilin/yilin/genome/iwgsc_cs/non_nucl_seq/cs_rRNA.split.slop1k.bed   -A > ${id}_$i.site_frq.bed
}
for i in  unique  best multiple_best
do
rm_rRNA $i &
done
wait
cat ${id}_best.site_frq.bed  ${id}_unique.site_frq.bed |sort -k1,1 -k2,2n |bedtools groupby -g 1,2,3,5,6 -c 4 -o sum|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6,$4,$5}' >  ${id}_unique_best.site_frq.bed
awk 'BEGIN{OFS="\t"}{print $1,$3,$6,$4}'  ${id}_unique_best.site_frq.bed > ${id}.ctss
genomesize=~/yilin/genome/iwgsc_cs161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta.size
id1=${id}_unique_best
awk '$6=="+"' ${id1}.site_frq.bed  > ${id1}_site_frq_plus.bed
awk '$6=="-"' ${id1}.site_frq.bed  > ${id1}_site_frq_minus.bed
cut -f 1-4  ${id1}_site_frq_plus.bed >  ${id1}_site_frq_plus.bedgraph
cut -f 1-4  ${id1}_site_frq_minus.bed > ${id1}_site_frq_minus.bedgraph
bedGraphToBigWig ${id1}_site_frq_plus.bedgraph $genomesize ${id1}_site_frq_plus.bw
bedGraphToBigWig ${id1}_site_frq_minus.bedgraph  $genomesize ${id1}_site_frq_minus.bw


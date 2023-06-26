
awk '$3=="long_terminal_repeat"' cs_all_LTRharvest.arrange.gff3 |cut -f 1,4,5,9 -|sort -k4,4 |paste - - |awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3,$5,$6,$7}' >  ltr_pair.txt

cat ltr_pair.txt|while read line
do
	n=`expr $n + 1`
	echo $line|tr ' ' '\t'|cut -f 2-4 > ${i}.pair.1.bed
	echo $line|tr ' ' '\t'|cut -f 5-7  > ${i}.pair.2.bed
	bedtools getfasta -fi ~/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta  -bed ${i}.pair.1.bed  > ${i}.pair.1.fa  &
	bedtools getfasta -fi ~/yilin/genome/iwgsc_cs/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta  -bed ${i}.pair.2.bed  > ${i}.pair.2.fa  &
	wait
	cat ${i}.pair.1.fa ${i}.pair.2.fa >  ${i}.pair.fa
	prank -d=${i}.pair.fa  -o=${i}_prank  -f=nexus -DNA
	~/miniconda2/pkgs/emboss-6.6.0-h6debe1e_0/bin/distmat   ${i}_prank.best.nex  -nucmethod 2 -outfile ${i}_prank.dist
	dist=`tail -n 2 ${i}_prank.dist |head -n 1 |cut -f 3`
	echo ${i} $n $dist   >> ${i}_dist.txt
done


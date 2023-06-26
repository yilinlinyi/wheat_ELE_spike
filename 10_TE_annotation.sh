##CLARITE
RepeatMasker -pa 10 wheat4_v2.fa -gff -nolow -xm -lib /picb/rsgeno/jingyu/wheat/TE_polymorphism/00genome/clariTeRep.fna -e ncbi

###get required file format for clari_te
genome=$sp.fa
genomesize=$sp.fa.fai
repeatmasker_gff=$sp.fa.out.gff
xm_file=$sp.fa.out.xm
for i in `cat Chr.txt`
do
{
 ###get gene annotation in embl format
    echo "get embl"
    echo ${i}${j} |seqtk subseq ${genome} - >${i}${j}.fa
    seqret -sequence ${i}${j}.fa -outseq ${i}${j}.embl -osformat embl
    embl_file=${i}${j}_anno.embl
    head -n 2 ${i}${j}.embl >${i}${j}_anno.embl
    echo -e "AC   unknown"";""\n""XX""\n""XX""\n""FH   Key            Location/Qualifiers""\n""FH" >> $embl_file
    count=0
    grep ${i}${j} $repeatmasker_gff |while read line
    do
        count=$(( $count + 1 ))
        location=`echo $line |awk 'BEGIN{OFS=".."}{print$4,$5}'`
        locustag=`echo $line |awk -v count=$count 'BEGIN{OFS="_"}{print$1,"REPEATMASKER",count}'`
        id=`echo $line |awk -v count=$count 'BEGIN{OFS="_"}{print$1,"REPEATMASKER",$4,$5,"Repeat_region",count}'`
        echo -e "FT   repeat_region   $location""\n""FT                   /locus_tag=\"$locustag\"""\n""FT                   /id=\"$id\"" >> $embl_file
    done
    echo "XX" >> $embl_file
    tail -n +3 ${i}${j}.embl >> $embl_file
#    rm -f ${i}${j}.embl
    ###repeatmasker
    echo "get xm_file"
    grep ${i}${j} $xm_file | grep -v "Simple_repeat" |sed 's/Unspecified/Unknown/g' >${i}${j}.fas.out.xm

    ###CLARI_TE
    echo "clari_te"
    perl /picb/rsgeno/jingyu/software/CLARI-TE/clariTE.pl -LTR clariTeRep_LTRposition.txt -gene $embl_file -classi clariTeRep_classification.txt -fasta ${i}${j}.fa -dir ./ ${i}${j}.fas.out.xm
} 
done


###de novo annotate full length LTR

gt suffixerator -db $sp.fa  -indexname  CS_1.0.fa  -tis -suf -lcp -des -dna
gt ltrharvest -index CS_1.0.fa -overlaps best -seed 30 -minlenltr 100 -maxlenltr 2000 -mindistltr 3000 -maxdistltr 25000 -similar 85 -mintsd 4 -maxtsd 20 -motif tgca -motifmis 1 -vic 60 -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3   -longoutput -out cs_1.0_chr_all.fsa -gff3  cs_all_LTRharvest.gff3
gt gff3 -sort cs_all_LTRharvest.gff3  > cs_all_LTRharvest.sorted.gff3
gt -j 8  ltrdigest  -trnas /public/home/xieyilin/yilin/TE_annotation/plantRNA-result.fasta  -hmms Pfam-A.hmm  -outfileprefix  cs_ltr  cs_all_LTRharvest.sorted.gff3 CS_1.0.fa > cs_all_LTRdigest.gff3

#gff arrangement
n=`less cs_all_LTRharvest.gff3|head -n 400 |grep "##seq"|wc -l`
less cs_all_LTRharvest.gff3 |grep "#[a-z]"|sed '1d'|head -n $n |cut -d ' '  -f 4 > id.1
less cs_all_LTRharvest.gff3 |grep "#[a-z]"|sed '1d'|tail -n $n > id.2
paste id.1  id.2 |sort -k1b,1  > id
grep -v "#" cs_all_LTRharvest.gff3 |sort -k1b,1 |join id  -  |tr ' ' '\t'|sed 's/#//'|cut -f 2,3,4,5,6,7,8,9,10 > cs_all_LTRharvest.arrange.gff3  
sort -k1,1 -k4,4n  cs_all_LTRharvest.arrange.gff3 > cs_all_LTRharvest.sorted.gff3  
less cs_all_LTRharvest.sorted.gff3 |awk '$3=="repeat_region" '  > cs_all_LTR_flt.gff3 

###full length LTR classification
awk '$3=="repeat_region"'| cs_fl_LTR.gff3 |awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$9,$6,$7}'  > cs_fl_LTR_region.bed
bedtools getfasta -s -fi $sp.fa  -bed cs_fl_LTR_region.bed  > cs_fl_LTR_region.fa
seqtk seq clariTeRep.fna |paste - - |grep "classif=D"|tr '\t' '\n' > clariTeRep_tir.fna
makeblastdb -in clariTeRep.fna  -dbtype nucl
blastn -query cs_fl_LTR_region.fa    -db  clariTeRep.fna   -out cs_fl_ltr_vs_clari_db_best.tab -outfmt 7 -max_target_seqs 1 -num_threads 16

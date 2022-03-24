windowMaker -g ~/tools/dm6_chrom.sizes -w 25 > dm6_25bp_window.bed
windowMaker -g ~/tools/dm6_chrom.sizes -w 1000 > dm6_1kb_window.bed

while read line;
do
awk '{print $4"\t"$1}' $line.multi.bed | uniq -c > $line.multi2.txt
awk '{print $2}' $line.multi2.txt | uniq -u > $line.multi_ID2.txt
awk 'NR==FNR{Arr[$0]++;next} ($2 in Arr){print $2"\t"$1}' $line.multi_ID2.txt $line.multi2.txt > $line.multi_ID_count2.txt

#get forward and reverse strand bed files

awk '{if($6=="+") print}' $line.multi.bed > $line.multi_for.bed
awk '{if($6=="-") print}' $line.multi.bed > $line.multi_rev.bed

awk 'NR==FNR{Arr[$1]=$2;next} ($4 in Arr){print $1"\t"$2"\t"$3"\t"$4"\t"int(1000/Arr[$4])/1000"\t"$6}' $line.multi_ID_count2.txt $line.multi_for.bed > $line.multi_for_frac2.bed
awk 'NR==FNR{Arr[$1]=$2;next} ($4 in Arr){print $1"\t"$2"\t"$3"\t"$4"\t"int(1000/Arr[$4])/1000"\t"$6}' $line.multi_ID_count2.txt $line.multi_rev.bed > $line.multi_rev_frac2.bed
done < sample_ID.txt

while read line;
do
zcat $line.fastq.gz | awk '{if(NR%4==1) printf "%s\t",$1; else if(NR%4==2) printf "%s\n",$1;}' > $line.fasta_ID
sed -i 's/@//g' $line.fasta_ID
bash ./run_dist $line.multi_for_frac2.bed > $line.multi_for_frac2.dist
bash ./run_dist $line.multi_rev_frac2.bed > $line.multi_rev_frac2.dist
awk 'NR==FNR{Arr[$1]++;next} ($1 in Arr){print $0}' $line.multi_for_frac2.dist $line.fasta_ID > $line.multi_for_frac2_ID.txt
awk 'NR==FNR{Arr[$1]++;next} ($1 in Arr){print $0}' $line.multi_rev_frac2.dist $line.fasta_ID > $line.multi_rev_frac2_ID.txt
done < sample_ID.txt

#Generate a final statistics table for intra-chromosomal multi-mapped reads
while read line;
do
echo -e "READ_ID\tCHR\tPOS_MIN\tPOS_MAX\tDISTANCE\tSEQUENCE" > $line.multi_for_frac2_final.txt
awk 'NR==FNR{Arr[$1]=$2;next} ($1 in Arr){print $0"\t"Arr[$1]}' $line.multi_for_frac2_ID.txt $line.multi_for_frac2.dist >> $line.multi_for_frac2_final.txt
echo -e "READ_ID\tCHR\tPOS_MIN\tPOS_MAX\tDISTANCE\tSEQUENCE" > $line.multi_rev_frac2_final.txt
awk 'NR==FNR{Arr[$1]=$2;next} ($1 in Arr){print $0"\t"Arr[$1]}' $line.multi_rev_frac2_ID.txt $line.multi_rev_frac2.dist >> $line.multi_rev_frac2_final.txt
echo -e "READ_ID,SEQUENCE,MAX_DISTANCE,CHR,POS_FIRST,POS_LAST" > $line.multi_for_frac2_final.csv
echo -e "READ_ID,SEQUENCE,MAX_DISTANCE,CHR,POS_FIRST,POS_LAST" > $line.multi_rev_frac2_final.csv
awk '{print $1","$6","$5","$2","$3","$4}' $line.multi_for_frac2_final.txt >> $line.multi_for_frac2_final.csv
awk '{print $1","$6","$5","$2","$3","$4}' $line.multi_rev_frac2_final.txt >> $line.multi_rev_frac2_final.csv
done < sample_ID.txt

#Select multi-mapped intachromosomal reads with max distance <=2Mb, and merge with uniquely mapped reads
while read line;
do
awk '{if($5<=2000000) print $1}' $line.multi_for_frac2_final.txt > $line.multi_for_frac2_ID_2Mb.txt
awk 'NR==FNR{Arr[$1]++;next} ($4 in Arr){print $0}' $line.multi_for_frac2_ID_2Mb.txt $line.multi_for_frac2.bed > $line.multi_for_frac2_2Mb.bed
awk '{if($5<=2000000) print $1}' $line.multi_rev_frac2_final.txt > $line.multi_rev_frac2_ID_2Mb.txt
awk 'NR==FNR{Arr[$1]++;next} ($4 in Arr){print $0}' $line.multi_rev_frac2_ID_2Mb.txt $line.multi_rev_frac2.bed > $line.multi_rev_frac2_2Mb.bed

bamToBed -i $line.uniq.bam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t1\t"$6}'> $line.uniq_frac.bed
awk '{if($6=="+") print}' $line.uniq_frac.bed > $line.uniq_for_frac.bed
awk '{if($6=="-") print}' $line.uniq_frac.bed > $line.uniq_rev_frac.bed

cat $line.uniq_for_frac.bed $line.multi_for_frac2_2Mb.bed > $line.for_frac.bed
cat $line.uniq_rev_frac.bed $line.multi_rev_frac2_2Mb.bed > $line.rev_frac.bed

sortBed -i $line.for_frac.bed > $line.for_frac_sort.bed
sortBed -i $line.rev_frac.bed > $line.rev_frac_sort.bed

mapBed -null 0 -a dm6_25bp_window.bed -b $line.for_frac_sort.bed | awk '{if($4 >0 ) print}' > $line.for_frac_25bp_bins.bdg
mapBed -null 0 -a dm6_1kb_window.bed -b $line.for_frac_sort.bed | awk '{if($4 >0 ) print}' > $line.for_frac_1kb_bins.bdg

mapBed -null 0 -a dm6_25bp_window.bed -b $line.rev_frac_sort.bed | awk '{if($4 >0 ) print}' > $line.rev_frac_25bp_bins.bdg
mapBed -null 0 -a dm6_1kb_window.bed -b $line.rev_frac_sort.bed | awk '{if($4 >0 ) print}' > $line.rev_frac_1kb_bins.bdg
done < sample_ID.txt

while read line;
do
for=$(awk '{s+=$5} END {print s}' $line.for_frac_sort.bed)
awk -v tot="$for" '{print $1"\t"$2"\t"$3"\t"$4/tot*10000000}' $line.for_frac_25bp_bins.bdg > $line.for_frac_25bp_bins_10M.bdg
awk -v tot="$for" '{print $1"\t"$2"\t"$3"\t"$4/tot*10000000}' $line.for_frac_1kb_bins.bdg > $line.for_frac_1kb_bins_10M.bdg

rev=$(awk '{s+=$5} END {print s}' $line.rev_frac_sort.bed)
awk -v tot="$rev" '{print $1"\t"$2"\t"$3"\t"$4/tot*10000000}' $line.rev_frac_25bp_bins.bdg > $line.rev_frac_25bp_bins_10M.bdg
awk -v tot="$rev" '{print $1"\t"$2"\t"$3"\t"$4/tot*10000000}' $line.rev_frac_1kb_bins.bdg > $line.rev_frac_1kb_bins_10M.bdg
done < sample_ID.txt


while read line;
do
# cutadapt -j 32 -m 20 -M 22 -o ${line}_20_22nt.fastq.gz $line.fastq.gz
# cutadapt -j 32 -m 23 -M 29 -o ${line}_23_29nt.fastq.gz $line.fastq.gz
zcat ${line}_20_22nt.fastq.gz | awk '{if(NR%4==1) print $1}' > ${line}_20_22nt_SEQID.txt
zcat ${line}_23_29nt.fastq.gz | awk '{if(NR%4==1) print $1}' > ${line}_23_29nt_SEQID.txt
sed -i 's/@//g' ${line}_20_22nt_SEQID.txt
sed -i 's/@//g' ${line}_23_29nt_SEQID.txt
awk 'NR==FNR{Arr[$1]++;next} ($4 in Arr){print $0}' ${line}_20_22nt_SEQID.txt $line.for_frac_sort.bed > $line.for_frac_sort_20_22nt.bed
awk 'NR==FNR{Arr[$1]++;next} ($4 in Arr){print $0}' ${line}_20_22nt_SEQID.txt $line.rev_frac_sort.bed > $line.rev_frac_sort_20_22nt.bed
awk 'NR==FNR{Arr[$1]++;next} ($4 in Arr){print $0}' ${line}_23_29nt_SEQID.txt $line.for_frac_sort.bed > $line.for_frac_sort_23_29nt.bed
awk 'NR==FNR{Arr[$1]++;next} ($4 in Arr){print $0}' ${line}_23_29nt_SEQID.txt $line.rev_frac_sort.bed > $line.rev_frac_sort_23_29nt.bed
mapBed -null 0 -a dm6_25bp_window.bed -b $line.for_frac_sort_20_22nt.bed | awk '{if($4 >0 ) print}' > $line.for_frac_25bp_bins_20_22nt.bdg
mapBed -null 0 -a dm6_1kb_window.bed -b $line.for_frac_sort_20_22nt.bed | awk '{if($4 >0 ) print}' > $line.for_frac_1kb_bins_20_22nt.bdg
mapBed -null 0 -a dm6_25bp_window.bed -b $line.rev_frac_sort_20_22nt.bed | awk '{if($4 >0 ) print}' > $line.rev_frac_25bp_bins_20_22nt.bdg
mapBed -null 0 -a dm6_1kb_window.bed -b $line.rev_frac_sort_20_22nt.bed | awk '{if($4 >0 ) print}' > $line.rev_frac_1kb_bins_20_22nt.bdg
mapBed -null 0 -a dm6_25bp_window.bed -b $line.for_frac_sort_23_29nt.bed | awk '{if($4 >0 ) print}' > $line.for_frac_25bp_bins_23_29nt.bdg
mapBed -null 0 -a dm6_1kb_window.bed -b $line.for_frac_sort_23_29nt.bed | awk '{if($4 >0 ) print}' > $line.for_frac_1kb_bins_23_29nt.bdg
mapBed -null 0 -a dm6_25bp_window.bed -b $line.rev_frac_sort_23_29nt.bed | awk '{if($4 >0 ) print}' > $line.rev_frac_25bp_bins_23_29nt.bdg
mapBed -null 0 -a dm6_1kb_window.bed -b $line.rev_frac_sort_23_29nt.bed | awk '{if($4 >0 ) print}' > $line.rev_frac_1kb_bins_23_29nt.bdg
done < sample_ID.txt

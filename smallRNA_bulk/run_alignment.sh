while read line;
do
bowtie2 -x /home/fgao/reference_genome/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome -U ${line}.fastq.gz -S $line.sam -p 32 -a
awk '{if($_~"^@" || $_~"XS:i:") print }' $line.sam > $line.multi.sam
awk '{if($_!~"XS:i:") print }' $line.sam > $line.uniq.sam

samtools view -bS $line.multi.sam > $line.multi.bam
samtools sort -O BAM -o $line.multi.sort.bam $line.multi.bam
samtools index $line.multi.sort.bam
samtools view -bS $line.uniq.sam > $line.uniq.bam
samtools sort -O BAM -o $line.uniq.sort.bam $line.uniq.bam
samtools index $line.uniq.sort.bam

bamToBed -i ${line}.multi.bam > ${line}.multi.bed

done < sample_ID_sel.txt

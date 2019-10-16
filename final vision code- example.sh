#!/usr/bin/bash

cp -r /localdisk/data/BPSM/Assignment1 ~

# fastqc quality control
fastqc /localdisk/data/BPSM/Assignment1/fastq/*fq.gz ./

# copy the genome and bowtie2 build_index
cp /localdisk/data/BPSM/Assignment1/Tbb_genome/Tb927_genome.fasta.gz ./
# unzip the fasta.gz to make sure this can be indexed
gunzip Tb927_genome.fasta.gz
bowtie2-build Tb927_genome.fasta index

#bowtie2 alignment
bowtie2 -x index -1 216_L8_1.fq.gz -2 216_L8_2.fq.gz -S 216_L8.sam -p 6
bowtie2 -x index -1 218_L8_1.fq.gz -2 218_L8_2.fq.gz -S 218_L8.sam -p 6
bowtie2 -x index -1 219_L8_1.fq.gz -2 219_L8_2.fq.gz -S 219_L8.sam -p 6
bowtie2 -x index -1 221_L8_1.fq.gz -2 221_L8_2.fq.gz -S 221_L8.sam -p 6
bowtie2 -x index -1 222_L8_1.fq.gz -2 222_L8_2.fq.gz -S 222_L8.sam -p 6
bowtie2 -x index -1 220_L8_1.fq.gz -2 220_L8_2.fq.gz -S 220_L8.sam -p 6

# samtools - converting the sam to indexed bam
samtools view -bS 216_L8.sam > 216_L8.bam
samtools sort -o 216_L8_sorted.bam 216_L8.bam
samtools index 216_L8_sorted.bam

samtools view -bS 218_L8.sam > 218_L8.bam
samtools sort -o 218_L8_sorted.bam 218_L8.bam
samtools index 218_L8_sorted.bam

samtools view -bS 219_L8.sam > 219_L8.bam
samtools sort -o 219_L8_sorted.bam 219_L8.bam
samtools index 219_L8_sorted.bam

samtools view -bS 220_L8.sam > 220_L8.bam
samtools sort -o 220_L8_sorted.bam 220_L8.bam
samtools index 220_L8_sorted.bam

samtools view -bS 221_L8.sam > 221_L8.bam
samtools sort -o 221_L8_sorted.bam 221_L8.bam
samtools index 221_L8_sorted.bam

samtools view -bS 222_L8.sam > 222_L8.bam
samtools sort -o 222_L8_sorted.bam 222_L8.bam
samtools index 222_L8_sorted.bam

# bedtools - generate the counts data
cp /localdisk/data/BPSM/Assignment1/Tbbgenes.bed ./
bedtools multicov -bams 216_L8_sorted.bam -bed Tbbgenes.bed > Slender_216_counts.txt
bedtools multicov -bams 218_L8_sorted.bam -bed Tbbgenes.bed > Slender_218_counts.txt
bedtools multicov -bams 219_L8_sorted.bam -bed Tbbgenes.bed > Slender_219_counts.txt
bedtools multicov -bams 220_L8_sorted.bam -bed Tbbgenes.bed > Stumpy_220_counts.txt
bedtools multicov -bams 221_L8_sorted.bam -bed Tbbgenes.bed > Stumpy_221_counts.txt
bedtools multicov -bams 222_L8_sorted.bam -bed Tbbgenes.bed > Stumpy_222_counts.txt

# statistical analysis
cut -f 7 Slender_216_counts.txt >> Slender_216_counts_only.txt
cut -f 7 Slender_218_counts.txt >> Slender_218_counts_only.txt
cut -f 7 Slender_219_counts.txt >> Slender_219_counts_only.txt
cut -f 7 Stumpy_220_counts.txt >> Stumpy_220_counts_only.txt
cut -f 7 Stumpy_221_counts.txt >> Stumpy_221_counts_only.txt
cut -f 7 Stumpy_222_counts.txt >> Stumpy_222_counts_only.txt

paste Slender_216_counts_only.txt Slender_218_counts_only.txt Slender_219_counts_only.txt >> Slender_counts_1.txt

paste Stumpy_220_counts_only.txt Stumpy_221_counts_only.txt Stumpy_222_counts_only.txt >> Stumpy_counts_1.txt

# Calculate the statistical average for each life cycle
while read v1 v2 v3
do
    sum=$(($v1+$v2+$v3))
    average=$((${sum}/3))
    echo -e "Slender_average\t${average}" >> Slender_average.txt
done < Slender_counts_1.txt

while read v1 v2 v3
do
    sum=$(($v1+$v2+$v3))
    average=$((${sum}/3))
    echo -e "Stumpy_average\t${average}" >> Stumpy_average.txt
done < Stumpy_counts_1.txt

# Slender_average.txt Stumpy_average.txt

cut -f 4 Slender_218_counts.txt >> gene.txt
paste gene.txt Slender_average.txt Stumpy_average.txt >> final.txt
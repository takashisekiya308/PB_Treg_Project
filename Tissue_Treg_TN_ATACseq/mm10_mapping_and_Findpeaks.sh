#!/bin/sh
fastq=*_R1.fastq
for i in $fastq
do
yolo=$(echo "$i" | rev | cut -c 10- | rev)
fastq_quality_trimmer -Q33 -t 20 -l 30 -i ${yolo}_R1.fastq -o ${yolo}_R1_trimmed.fastq
fastq_quality_trimmer -Q33 -t 20 -l 30 -i ${yolo}_R2.fastq -o ${yolo}_R2_trimmed.fastq
bowtie2 -p 8 -q -N 0 -x /Users/sekiya/Desktop/reference_genome/mm10_ebwt/mm10 -1 ${yolo}_R1_trimmed.fastq -2 ${yolo}_R2_trimmed.fastq -S ${yolo}.sam
grep -v "XS:" ${yolo}.sam > ${yolo}_unique.sam
samtools view -S -b ${yolo}_unique.sam >${yolo}.bam
bedtools intersect -a ${yolo}.bam -b /Users/sekiya/Desktop/reference_genome/mm10-blacklist.v2.bed -v > ${yolo}_blacklist.bam
samtools sort ${yolo}_blacklist.bam > ${yolo}_blacklist_sorted.bam
samtools index ${yolo}_blacklist_sorted.bam
makeTagDirectory ${yolo}_homer -single ${yolo}_blacklist_sorted.bam
findPeaks ${yolo}_homer -style factor -o auto
cd ${yolo}_homer
sed '/^#/d' peaks.txt | awk -v 'OFS=\t' '{print $2, $3, $4}' > peaks.bed
cd ..
genomeCoverageBed -ibam ${yolo}_blacklist_sorted.bam -bg -trackline -trackopts 'name="${yolo}_blacklist_sorted" color=250,0,0' > ${yolo}_blacklist_sorted.bedGraph
done

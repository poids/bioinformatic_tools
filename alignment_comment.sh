#!/bin/bash

#Description:
#Script aligns fastq files in zipped format with HISAT2
#Converts Sam output into BAM
#Then creates BAM Index Files for visualization (.bai)

#Written & Maintained by Vasco Morais (April 2019)

echo $RNA_HOME

export RNA_PRACTICE_DATA_DIR=$RNA_HOME/practice
echo $RNA_PRACTICE_DATA_DIR

RNA_PATH=$RNA_HOME/practice
ADAPTER=${RNA_REFS_DIR}'/illumina_multiplex.fa'
MIX='_r'
OUTPUT='_Build37-ErccTranscripts-chr22'

echo $RNA_PATH
echo $ADAPTER


echo '-------------------------------'
echo 'BEGIN SCRIPT'

for m in 1 2
do
if [ $m = 1 ]; then
	REP='hcc1395_normal_rep'
	ID='HCC1395_normal_rep'
	OUTPUT='HCC1395_normal.bam'
elif [ $m = 2 ]; then
	REP='hcc1395_tumor_rep'
	ID='HCC1395_tumor_rep'
	OUTPUT='HCC1395_tumor.bam'
fi
for r in 1 2 3;
do echo 'Mix: '${m};
echo 'Read: '${r};i
echo 'Alignment:'
hisat2 -p 8 --rg-id=${ID}${r} --rg SM:${ID}${r} --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/trimmed/${REP}${r}_r${m}_1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/trimmed/${REP}${r}_r${m}_2.fastq.gz -S ./${ID}${r}.sam
echo '---------';
echo 'Sam2Bam'
samtools sort -@ 8 -o ${ID}${r}.bam ${ID}${r}.sam
echo '---------';
done;
echo 'Merge HISAT2 BAM Files'
java -Xmx2g -jar $RNA_HOME/student_tools/picard.jar MergeSamFiles OUTPUT=$OUTPUT INPUT=${ID}1.bam INPUT=${ID}2.bam INPUT=${ID}3.bam

echo '---------';
echo 'Indexing BAM Files'
echo '---------';

find $OUTPUT -exec echo samtools index {} \; | sh
done;

echo finished!

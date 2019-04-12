#!/bin/bash

#DESCRIPTION:
#Tool for analyzing readcount/pileup over a genomic region of interest

#REQUIREMENTS:
#samtools
#reference genome (.fasta)
#alignment mapping data (.bam)
#genomic region of interest (chr:start-stop)

#Written and Maintained by Vasco Morais (April 2019)

while [ -n "$1" ]; do # while loop starts
    case "$1" in
    -ref)
        REF="$1"
        echo "Reference Genome: ${REF}"
        shift
        ;;
    -chr)
        CHR="$2"
        echo "Chromosome: ${REF}"
        shift
        ;;
    -start)
        START="$3"
        echo "Start: ${START}"
        shift
        ;;
    -stop)
        STOP="$4"
        echo "Stop: ${STOP}"
        shift
        ;;
    *) echo "Option $1 not recognized"
    esac
    shift
done


# find *.bam -exec echo samtools index {} \; | sh

# samtools faidx $RNA_REF_FASTA


# samtools mpileup -f $RNA_REF_FASTA -r 22:18918457-18918467 $RNA_ALIGN_DIR/UHR.bam $RNA_ALIGN_DIR/HBR.bam

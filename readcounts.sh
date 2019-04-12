#!/bin/bash

#DESCRIPTION:
#Tool for analyzing readcount/pileup over a genomic region of interest

#REQUIREMENTS:
#samtools
#reference genome (.fasta)
#alignment mapping data (.bam)
#genomic region of interest (chr:start-stop)

#Written and Maintained by Vasco Morais (April 2019)

while (( "$#" )); do # while loop starts
    case "$1" in
        -ref)
            REF="$2"
            echo "Reference Genome: ${REF}"
            shift 2
            ;;
        -chr)
            CHR="$2"
            echo "Chromosome: ${CHR}"
            shift 2
            ;;
        -start)
            START="$2"
            echo "Start: ${START}"
            shift 2
            ;;
        -stop)
            STOP="$2"
            echo "Stop: ${STOP}"
            shift 2
            ;;
        -read)

        --) # end argument parsing
            shift
            break
            ;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1
            ;;
        *) # preserve positional arguments
            READS="$1"
            echo "Option $1 not recognized"
            shift
            ;;
    esac
    # shift
done


find .bam -exec echo samtools index {} \; | sh

samtools faidx ${REF}


# samtools mpileup -f $RNA_REF_FASTA -r 22:18918457-18918467 $RNA_ALIGN_DIR/UHR.bam $RNA_ALIGN_DIR/HBR.bam

perl -ne '@path=split("/", $_); @output=split(".bam", $path[-1]); printf "$output[0]\n";'
#!/bin/bash

#DESCRIPTION:
#Tool for analyzing readcount/pileup over a genomic region of interest and producing per-base readcount

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
        -pos)
            POS="$2"
            CHR=$(echo ${POS} | perl -ne '@chr=split(":", $_); printf "$chr[0]";')
            START=$(echo ${POS} | perl -ne '@pos=split(":", $_); @start=split("-", $pos[1]); printf "$start[0]";')
            STOP=$(echo ${POS} | perl -ne '@pos=split(":", $_); @stop=split("-", $pos[1]); printf "$stop[0]";')
            echo "Chromosome: ${CHR}"
            echo "Start: ${START}"
            echo "Stop: ${STOP}"
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
        -i)
            INPUT="$2"
            OUTPUT=$(echo ${INPUT} | perl -ne '@path=split("/", $_); @output=split(".bam", $path[-1]); printf "$output[0]\n";')
            echo "Input: ${INPUT}"
            echo "Output: ${OUTPUT}"
            shift 2
            ;;
        --) # end argument parsing
            shift
            break
            ;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1
            ;;
        *) # preserve positional arguments
            #READS="$1"
            echo "Option $1 not recognized"
            shift
            ;;
    esac
    # shift
done

### Index Bam ###
echo '---'
echo 'INDEX BAM FILE'
echo '---'

samtools index ${INPUT} ${INPUT}.bai

#find .bam -exec echo samtools index {} \; | sh

### Build Reference Genome Index ###
echo '---'
echo 'BUILD REF GENOME INDEX'
echo '---'

samtools faidx ${REF}

### Count Pile-Up on ROI ###
echo '---'
echo 'COUNT PILE-UP ON REGION OF INTEREST'
echo '---'

samtools mpileup -f $RNA_REF_FASTA -r ${CHR}:${START}-${STOP} ${INPUT}

#DOES IT REQUIRE MORE THAN ONE BAM FILE???

### Perform BAM read count on ROI ###
echo '---'
echo 'PERFORM BAM READCOUNT ON REGION OF INTEREST'
echo '---'

#Create bedfile with genomic coordinates in header
echo "${CHR} ${START} ${STOP}" > ${OUTPUT}_snvs.bed
#Run bam-readcout
bam-readcount -l ${OUTPUT}_snvs.bed -f ${REF} ${INPUT} 2>/dev/null 1>${OUTPUT}_bam-readcounts.txt

#Calculate per-base readcount
cat UHR_bam-readcounts.txt | perl -ne '@data=split("\t", $_); @Adata=split(":", $data[5]); @Cdata=split(":", $data[6]); @Gdata=split(":", $data[7]); @Tdata=split(":", $data[8]); print "UHR Counts\t$data[0]\t$data[1]\tA: $Adata[1]\tC: $Cdata[1]\tT: $Tdata[1]\tG: $Gdata[1]\n";' > ${OUTPUT}_per-base-readcount.txt

#Show per-base results
cat ${OUTPUT}_per-base-readcount.txt



#Remove temporary files?
rm ${INPUT}.bai
rm ${REF}.fai
rm ${OUTPUT}_snvs.bed
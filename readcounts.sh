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
        -h|--help)
            echo 'Usage: [OPTIONS] [PARAMS]... [BAM FILE]'
            echo 'Description: Tool for analyzing readcount/pileup over a genomic region of interest and producing per-base readcount'
            echo ''
            echo '  readcounts.sh [options]* -ref <ref-genome.fa> {-pos | [-chr -start -stop]} -i <input1.bam> <input2.bam>...'
            echo ''
            echo '  Input:'
            echo '      -ref            Path/to/reference/genome (.fasta or .fa)'
            echo '      -pos            Genomic coordinates <chr:start-stop> (e.g., 22:15867-15997)'
            echo '      -chr            Chromsome (e.g., 22)'
            echo '      -start          Starting Genomic Coordinate'
            echo '      -stop           Ending Genomic Coordinate'
            echo '      -i              Path/to/BAMfiles (.bam) - enter one or multiple files'
            echo ''
            echo '**Note user can either provide genomic coordinate using or -pos flag'
            echo '  with standard coordinate format or by providing each -chr, -start,'
            echo '  and -stop argument individually'
            echo ''
            echo '  Options:'
            echo '      -r              Removes fasta (.fai) and bam (.bai) indexes generated during script'
            echo '      -h, --help      Display help'
            echo ''
            echo 'Written and Maintained by Vasco Morais (April 2019)'
            exit 1
            ;;
        -r) #Removes termporary files if included
            RM=true
            echo "Will remove temporary files"
            shift 1
            ;;
        -ref)
            REF="$2"
            shift 2
            ;;
        -pos)
            POS="$2"
            CHR=$(echo ${POS} | perl -ne '@chr=split(":", $_); printf "$chr[0]";')
            START=$(echo ${POS} | perl -ne '@pos=split(":", $_); @start=split("-", $pos[1]); printf "$start[0]";')
            STOP=$(echo ${POS} | perl -ne '@pos=split(":", $_); @stop=split("-", $pos[1]); printf "$stop[1]";')
            shift 2
            ;;
        -chr)
            CHR="$2"
            shift 2
            ;;
        -start)
            START="$2"
            shift 2
            ;;
        -stop)
            STOP="$2"
            shift 2
            ;;
        -i)
            shift 1
            if [ $# -gt 1 ] ; then
                INPUT=$@
            else
                INPUT="$1"
            fi
            shift $#
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
done

for i in $INPUT; do
    OUTPUT=$(echo ${i} | perl -ne '@path=split("/", $_); @output=split(".bam", $path[-1]); printf "$output[0]\n";')

    echo "Reference Genome: ${REF}"
    echo "Chromosome: ${CHR}"
    echo "Start: ${START}"
    echo "Stop: ${STOP}"
    echo "Input: ${i}"
    echo "Output: ${OUTPUT}"


    ### Index Bam ###
    echo '---'
    echo 'INDEX BAM FILE'
    echo '---'

    samtools index ${i} ${i}.bai

    ### Build Reference Genome Index ###
    echo '---'
    echo 'BUILD REF GENOME INDEX'
    echo '---'

    samtools faidx ${REF}

    ### Count Pile-Up on ROI ###
    echo '---'
    echo 'COUNT PILE-UP ON REGION OF INTEREST'
    echo '---'

    samtools mpileup -f $RNA_REF_FASTA -r ${CHR}:${START}-${STOP} ${i}

    ### Perform BAM read count on ROI ###
    echo '---'
    echo 'PERFORM BAM READCOUNT ON REGION OF INTEREST'
    echo '---'

    #Create bedfile with genomic coordinates in header
    echo "${CHR} ${START} ${STOP}" > ${OUTPUT}_snvs.bed
    #Run bam-readcout
    bam-readcount -l ${OUTPUT}_snvs.bed -f ${REF} ${i} 2>/dev/null 1>${OUTPUT}_bam-readcounts.txt

    #Calculate per-base readcount
    cat ${OUTPUT}_bam-readcounts.txt | perl -ne '@data=split("\t", $_); @Adata=split(":", $data[5]); @Cdata=split(":", $data[6]); @Gdata=split(":", $data[7]); @Tdata=split(":", $data[8]); print "Counts\t$data[0]\t$data[1]\tA: $Adata[1]\tC: $Cdata[1]\tT: $Tdata[1]\tG: $Gdata[1]\n";' > ${OUTPUT}_per-base-readcount.txt

    #Show per-base results
    cat ${OUTPUT}_per-base-readcount.txt

    #Remove indexes?
    if [ "$RM" = true ] ; then
        rm ${i}.bai
        rm ${REF}.fai
    fi

    #Remove bedfile
    rm ${OUTPUT}_snvs.bed
done
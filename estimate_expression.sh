#!/bin/bash -e

#DESCRIPTION:
#Used to generate expression estimates from SAM/BAM alignment files

#REQUIREMENTS:
#Stringtie
#HTSeq-Count
#gene annotations (.gtf/.gff)
#alignment mapping data (.bam/.sam)

#Use Stringtie package to generate TPM/FPKM/Coverage matrices
#Use HTSeq-Count to generate raw counts matrix

#Written and Maintained by Vasco Morais (June 2019)

while (( "$#" )); do # while loop starts
	case "$1" in
		-h|--help
			echo 'help'
			exit 1
			;;
		-s) #Use Stringtie
			STRING=true
			S_PATH="./expression/stringtie/ref_only"
			mkdir -p ${S_PATH} #create output folder
			shift 1
			;;
		-h) #Use HTSeq-Count
			HTSEQ=true
			H_PATH="./expression/htseq_counts"
			mkdir -p ${H_PATH} #create output folder
			shift 1
			;;
		-r) #GTF/GFF annotations
			REF="$2"
			shift 2
			;;
		-a) #Alignment files (.bam/.sam)
			ALIGN_DIR="$2"
			shift 2
			;;
		--) # end argument parsing
			shift
			break
			;;
		-*|--*=) # unsupported flags
			echo "Error: unsupported flag $1" >&2
			exit 1
			;;
		*) # preserve positional arguments
			#Maybe not needed because one of the flags will be default
			echo "Option $1 not recognized"
			shift
			;;
	esac
done	

#Required conditons to run script
if [[ $STRING != true && $HTSEQ != true ]]; then
	echo 'No tool selected'
	exit 1
fi

INPUT=`ls ${ALIGN_DIR} | grep '[0-9].bam$' | cut -d '.' -f 1`

#Check to see if alignment files are present
if [ 0 -lt $(echo ${INPUT} | wc -w) ]; then
	echo 'Alignment Files Found'
else
	echo 'Alignment Files Not Found!'
	exit 1
fi

#Checks to see if annotation file is present
if [ -f ${REF} ]; then
	echo 'Annotation file found'
else
	echo "Annotation file not found"
	exit 1
fi

echo '====================='
echo 'Estimating Expression'
echo '====================='

echo ${INPUT}
echo '----------'

for i in ${INPUT}
do
	echo ${i}
	stringtie -p 8 -G ${REF} -e -B -o ${S_PATH}/${i}/transcripts.gtf -A ${S_PATH}/${i}/gene_abundances.tsv ${ALIGN_DIR}/${i}.bam
done

echo '============================'
echo 'Building Expression Matrices'
echo '============================'

#Switching to output directory
cd ${S_PATH}

echo '---'
#Downloads script required to build matrices
echo 'Downloading expression matrix script'
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/stringtie_expression_matrix.pl
chmod +x stringtie_expression_matrix.pl
echo '---'

#Makes comma separated list of all input files
ALIGN_LIST=`echo ${INPUT} | tr [:space:] , | sed 's/\(.*\),/\1 /'`

#Builds TPM Expression Matrix
echo 'Build TPM Expression Matrix'
./stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs=${ALIGN_LIST} --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv
echo '---'

#Builds FPKM Expression Matrix
echo 'Build FPKM Expression Matrix'
./stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs=${ALIGN_LIST} --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv
echo '---'

#Builds Coverage Matrix
echo 'Build Coverage Matrix'
./stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs=${ALIGN_LIST} --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv
echo '---'

#Removes stringtie_expression_matrix script
rm stringtie_expression_matrix.pl

#Switching back to main directory
cd ../../..

echo '========================='
echo 'Building Raw Count Matrix'
echo '========================='

#Switching to output directory
cd ${H_PATH}



#Switching back to main directory
cd ../..

echo 'finished!'
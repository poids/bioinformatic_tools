#!/bin/bash

#Requires RSeQC
#For pair-end reads

INPUT=$1
REF=$2
OUTPUT=$3

bam_stat.py -i ${INPUT}
clipping_profile.py -i ${INPUT} -o ${OUTPUT} -s "PE"
geneBody_coverage.py -r ${REF} -i ${INPUT} -o ${OUTPUT}
infer_experiment.py -r ${REF} -i ${INPUT}
inner_distance.py -r ${REF} -i ${INPUT} -o ${OUTPUT}
junction_annotation.py -r ${REF} -i ${INPUT} -o ${OUTPUT}
junction_saturation.py -r ${REF} -i ${INPUT} -o ${OUTPUT}
read_distribution.py -r ${REF} -i ${INPUT}
read_duplication.py -i ${INPUT} -o ${OUTPUT}
read_GC.py -i ${INPUT} -o ${OUTPUT}
read_NVC.py -i ${INPUT} -o ${OUTPUT}
read_quality.py -i ${INPUT} -o ${OUTPUT}
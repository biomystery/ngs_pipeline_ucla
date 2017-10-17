#!/bin/bash
if [ -z "$1" ]; then
    echo -e "\nUsage: filterfun.sh 1.bamfile\n"
    exit 1
elif [ -z "$2" ]; then
    echo -e "\nUsage: filterfun.sh 2.dir\n"
    exit 1
elif [ -z "$3" ]; then
    echo -e "\nUsage: filterfun.sh 3.processor num\n"
    exit 1
fi



prefix=".filtered"
filetype='.bam'
data=$1
WORKING_DIR=$2; NPROC_PER_SAMPLE=$3
nameStr=$(echo "$data"| cut -f1 -d".")
#echo $nameStr
#echo $data
#echo $WORKING_DIR'/'$nameStr$prefix$filetype
#echo $NPROC_PER_SAMPLE
eval "samtools view -b -F 2820 -q 30 -@ $NPROC_PER_SAMPLE $data > $WORKING_DIR/$nameStr$prefix$filetype"



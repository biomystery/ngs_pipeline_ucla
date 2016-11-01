#!/bin/bash
echo -e "############################################################"
echo -e "(`date`) Welcome to the RNAseq Auto analysis ver 0.0"
echo -e "############################################################"


echo -e "(`date`) Initating input parameters ....."
echo -e " "

#input 1: consolidated FASTQ Folder 
PARENT_DIR=$PWD #  project dir

# log files
LOG_FILE=$PARENT_DIR"/run.log.txt"; LOG_ERR_FILE=$PARENT_DIR"/run.err.txt"
echo -e "(`date`) Setting folder \n" | tee -a $LOG_FILE
echo -e "(`date`) Project folder is $PARENT_DIR \n" | tee -a $LOG_FILE
echo -e "(`date`) Log files are $LOG_ERR_FILE and $LOG_FILE \n" | tee -a $LOG_FILE
echo -e "------------------------------" | tee -a $LOG_FILE
echo -e "(`date`) NOW please upload the password.txt file \n" | tee -a $LOG_FILE
echo -e "(`date`) It should only contain 1 line with : separator \n" | tee -a $LOG_FILE
echo -e "(`date`) example: SxaQSEQsWA146L6:wK4rq3yX4Em7 \n" | tee -a $LOG_FILE

#read -p"Press [Enter] key to when you finished construct the file"

# barcode files
echo -e "------------------------------" | tee -a $LOG_FILE
echo -e "(`date`) NOW please upload the barcode.txt file \n" | tee -a $LOG_FILE
echo -e "(`date`) line format: samplename<tab>barcode \n" | tee -a $LOG_FILE
echo -e "(`date`) example: \n" | tee -a $LOG_FILE
echo -e "(`date`) example: \n" | tee -a $LOG_FILE
echo -e "(`date`) Wt-0	ATTCCTT
Wt-TNF-05	ACTGATA
Wt-TNF-1	GAGTGGA
Wt-TNF-3	CGTACGT
\n" | tee -a $LOG_FILE
#read -p"Press [Enter] key to when you finished construct the file"


BARCODE_FILE=$PARENT_DIR/barcode.txt
echo -e "(`date`) barcode file is $BARCODE_FILE \n" | tee -a $LOG_FILE

# number of samples 
SAMPLE_NO=`wc -l < $BARCODE_FILE` 
echo -e "There are $SAMPLE_NO samples in this experiment: \n" | tee -a $LOG_FILE

while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    barcode=$(echo -e "$line" | cut -f2 -d$'\t')
    echo "$prefix $barcode"
    done<$BARCODE_FILE | tee -a $LOG_FILE


#number  of processors per sample for fastqc 
NPROC_PER_SAMPLE=2
echo -e " $NPROC_PER_SAMPLE processors per sample \n" | tee -a $LOG_FILE
TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE)) # calculate number of total processor for the user
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | tee -a $LOG_FILE

echo -e "Check the password file (password.txt)" | tee -a $LOG_FILE
PASSWD_INFO_INPUT=$(head -n 1 password.txt)
echo -e "Your password is : $PASSWD_INFO_INPUT"


echo -e "(`date`)Starting Step 6a: bamcompare to input " | tee -a $LOG_FILE


STEP='06a_bamcompare'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;


LOG_FILE_STEP=$WORKING_DIR"/log.txt"
LOG_ERR_FILE=$WORKING_DIR"/err.txt"
date | tee -a $LOG_FILE_STEP
bamComparefun (){
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
	if (`echo $nameStr | grep -q b1`); then input=b1-input.bam; else input=b2-input.bam; fi
	bamCompare --bamfile1 "$data" --bamfile2 "$input" --binSize 10 --normalizeTo1x 2150570000  --smoothLength 30 --extendReads 200 -p $NPROC_PER_SAMPLE -o $WORKING_DIR/$nameStr".bw" --ignoreForNormalization chrX &
    done
}

cd $PARENT_DIR'/04a_deDupped'
ls *.bam | bamComparefun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP




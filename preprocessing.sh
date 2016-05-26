#!/bin/bash
echo -e "############################################################"
echo -e "(`date`) initating input parameters ....."
echo -e "############################################################"

#input 1: consolidated FASTQ Folder 

PARENT_DIR=$PWD #  project dir
# log files
LOG_FILE=$PARENT_DIR"/run.log.txt"; LOG_ERR_FILE=$PARENT_DIR"/run.err.txt"
echo -e "(`date`)Setting folder \n" | tee -a $LOG_FILE
echo -e "(`date`) project folder is $PARENT_DIR \n" | tee -a $LOG_FILE
echo -e "(`date`) log files are $LOG_ERR_FILE and $LOG_FILE \n" | tee -a $LOG_FILE

# barcode files
BARCODE_FILE=$PARENT_DIR/barcode.txt
echo -e "(`date`) barcode file is $BARCODE_FILE \n"

# samples 
SAMPLE_NO=`wc -l < $BARCODE_FILE` 
echo -e "There are $SAMPLE_NO samples in this experiment \n" | tee -a $LOG_FILE

#number  of processors per sample for fastqc 
NPROC_PER_SAMPLE=2
echo -e " $NPROC_PER_SAMPLE processors per sample \n" | tee -a $LOG_FILE
TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE)) # calculate number of total processor for the user
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | tee -a $LOG_FILE

echo -e "Please input the password info (example: )\n" | tee -a $LOG_FILE
PASSWD_INFO_INPUT="SxaQSEQsWA148L3:wa4kD9Ye8hr8"


echo -e "############################################################"| tee -a $LOG_FILE
echo -e "`date` Running the pipelines " | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE


echo -e "############################################################"| tee -a $LOG_FILE
echo -e "0.1 `date` starting downloading" | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

STEP="00raw"; WORK_DIR=$PARENT_DIR"/"$STEP
mkdir $WORK_DIR; cd $WORK_DIR

grab_bscrc.sh SxaQSEQsWA148L3:wa4kD9Ye8hr8 # why freeze here
wait;echo -e "(`date`) downloaded the raw qseq data" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "0.2 `date` starting decompressing " | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

RAW_DIR=$(echo $PASSWD_INFO_INPUT|cut -f1 -d":"); RAW_DIR=$WORK_DIR/$RAW_DIR; cd $RAW_DIR
ls -1 *.gz | parallel -j $TOTAL_PROC_NO --eta gunzip 
wait;echo -e "(`date`) decompressed the raw qseq data" | tee -a $LOG_FILE


echo -e "############################################################"| tee -a $LOG_FILE
echo -e "0.3 `date` starting converting qseq to fastq " | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

ls -1 * | while read data; do  echo ${data:0:10}; done | tee -a $LOG_FILE
ls -1 * | while read data; do  qseq2fastq.pl $data ${data:0:10}'.fastq' & done | tee -a $LOG_FILE
wait;echo -e "(`date`) 0.3 convering finished" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "0.4 `date` starting demultiplex " | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

STEP="consolidate";
WORK_DIR=$WORK_DIR'/'$STEP; mkdir $WORK_DIR

demultiplexFunMuticore (){
    count=1;ncore=$TOTAL_PROC_NO;
    prefix=$WORK_DIR"/"
    while read data
    do
        suffix=${data:5:11}
        echo 'processing:'$suffix
        idxfile=$(echo "$data" |sed -r 's/_1_/_2_/')
        echo "using idx file: "$idxfile
        echo "calulating no.$count sampling ......"
        if [ `echo $count" % "$ncore | bc` -eq 0 ] #mod ncore
        then
            cat $data | fastx_barcode_splitter.pl --bcfile $BARCODE_FILE --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1
        else
            cat $data | fastx_barcode_splitter.pl --bcfile $barcodefile --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1 &
        fi
        count=$((count+1))
    done
}

ls -1 *_1_*.fastq | demultiplexFunMuticore 1>>$LOG_FILE 2>>$LOG_ERR_FILE

wait;echo -e "(`date`) 0.4 demultiplex finished" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "0.5 `date` merge fastq " | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

cd $WORK_DIR
while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    #eval "ls ${prefix}* "
    echo $prefix
    eval "cat ${prefix}* > ${prefix}.fastq && rm ${prefix}_* &"
    #eval "rm ${prefix}_* " 
done<$BARCODEFILE | tee -a $LOG_FILE
wait;echo -e "(`date`) 0.5 demultiplex finished" | tee -a $LOG_FILE


#------------------------------------------------------------
# 5. remove tmp files & mv files to the server 
#------------------------------------------------------------
# define serverDir
#serverDir="/mnt/biggie/backed_up/frank/projects/dynamic-decoding/caRNA/batch3/2000/"
#rm *.fastq # rm tmp fastq files 
#ls -1 *.txt | parallel -j 48 --eta gzip -9  # gzip txt  raw file




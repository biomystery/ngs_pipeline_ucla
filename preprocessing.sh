#!/bin/bash

#------------------------------------------------------------
# use read function to take input 
#------------------------------------------------------------
echo "Please input the password info (example: )\n"


#------------------------------------------------------------
# input parameters, check
#------------------------------------------------------------


E_WRONG_ARGS=85
script_parameters="-w -p -m -z"
# -w = working directory
# -p = password inf 


if [ $# -ne 4 ]
then
  echo "Usage: `basename $0` $script_parameters"
  # `basename $0` is the script's filename.
  exit $E_WRONG_ARGS
fi



# define paths & steps here 

WORK_DIR="/home/frank/kim/2000/"
STEP0_DIR="00raw/"
barcodefile=./barcode.txt

############################################################
# 1. Download & convert
############################################################
mkdir cd $WORK_DIR$STEP0_DIR; cd $WORK_DIR$STEP0_DIR
grab_bscrc.sh SxaQSEQsWA148L3:wa4kD9Ye8hr8 # why freeze here
cd SxaQSEQsWA148L3
ls -1 *.gz | parallel -j 60 --eta gunzip

diff <(ls -1 *_1_*.txt | while read data; do echo ${data:6:4}; done) <(ls -1 *_2_*.txt | while read data; do echo ${data:6:4}; done)


############################################################
# 2.  convert qseq txt to fastq 
############################################################

ls -1 * | while read data; do  echo ${data:0:10}; done
ls -1 * | while read data; do  qseq2fastq.pl $data ${data:0:10}'.fastq' & done

# method 2, not working currently 
qseq2fastqfunc (){
    while read data
    do
        qseq2fastq.pl $data ${data:0:10}'.fastq'
    done
}

ls -1 *.txt | parallel -j 48 --eta qseq2fastqfunc


#------------------------------------------------------------
# 3 Demultiplex FASTQ
#------------------------------------------------------------
STEP="02_consolidate"
echo $STEP
prefix=$WORK_DIR$STEP'/'
mkdir $prefix

# method 1 / not working currently 
demultiplexFun (){
    while read data
    do
        suffix=${data:5:11}
        echo 'processing:'$suffix
        idxfile=$(echo "$data" |sed -r 's/_1_/_2_/')
        echo "using idx file: "$idxfile
        cat $data | fastx_barcode_splitter.pl --bcfile $barcodefile --prefix ../02_consolidated/ --suffix $suffix --idxfile $idxfile --mismatches 1 
    done
}
export -f demultiplexFun
ls -1 *_1_*.fastq | parallel -j 48 --eta demultiplexFun

# method 2 
demultiplexFunMuticore (){
    count=1;ncore=60;
    while read data
    do
        suffix=${data:5:11}
        echo 'processing:'$suffix
        idxfile=$(echo "$data" |sed -r 's/_1_/_2_/')
        echo "using idx file: "$idxfile
        echo "calulating no.$count sampling ......"
        if [ `echo $count" % "$ncore | bc` -eq 0 ] #mod ncore
        then
            cat $data | fastx_barcode_splitter.pl --bcfile $barcodefile --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1
        else
            cat $data | fastx_barcode_splitter.pl --bcfile $barcodefile --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1 &
        fi
        count=$((count+1))
    done
}

ls -1 *_1_*.fastq | demultiplexFunMuticore > demultiplex.log.txt

# check if all files demultiplexed
ls -1 ../../02_consolidate/Wt-TNF-8*.fastq | wc -l ;  ls -1 *_2_*.fastq | wc -l ; 
diff <(ls -1 ../../02_consolidate/Wt-TNF-8*.fastq |  while read data; do echo ${data: -11}; done) <(ls -1 *_2_*.fastq | while read data; do echo ${data:5:11}; done)

# head -c 2 # check the first two characters 

# 4. Merge FASTQ
barcodefile=../00raw/SxaQSEQsWA148L3/barcode.txt

cd ../../02_consolidate/
while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    #eval "ls ${prefix}* "
    echo $prefix
    eval "cat ${prefix}* > ${prefix}.fastq && rm ${prefix}_* &"
    #eval "rm ${prefix}_* " 
done<$barcodefile

# cat alpha-0_*.fastq >alpha-0.fastq
#------------------------------------------------------------
# 5. remove tmp files & mv files to the server 
#------------------------------------------------------------
# define serverDir
serverDir="/mnt/biggie/backed_up/frank/projects/dynamic-decoding/caRNA/batch3/2000/"
rm *.fastq # rm tmp fastq files 
ls -1 *.txt | parallel -j 48 --eta gzip -9  # gzip txt  raw file




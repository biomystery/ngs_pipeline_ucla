

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
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | te\Wt-TNF-1	GAGTGGA
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
NPROC_PER_SAMPLE=4
echo -e " $NPROC_PER_SAMPLE processors per sample \n" | tee -a $LOG_FILE
TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE)) # calculate number of total processor for the user
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | tee -a $LOG_FILE

echo -e "Check the password file (password.txt)" | tee -a $LOG_FILE
PASSWD_INFO_INPUT=$(head -n 1 password.txt)
echo -e "Your password is : $PASSWD_INFO_INPUT"

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) start running the pipelines " | tee -a $LOG_FILE
echo -e "(`date`) 1.1  starting downloading" | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

STEP="00raw"; WORKING_DIR=$PARENT_DIR"/"$STEP
LOG_FILE_STEP=$WORKING_DIR"/log.txt"
mkdir -p $WORKING_DIR; cd $WORKING_DIR
date| tee -a $LOG_FILE_STEP
grab_bscrc.sh $PASSWD_INFO_INPUT | tee - a $LOG_FILE_STEP

wait;echo -e "(`date`) downloaded the raw qseq data" | tee -a $LOG_FILE_STEP
RAW_DIR=$(echo $PASSWD_INFO_INPUT|cut -f1 -d":"); RAW_DIR=$WORKING_DIR/$RAW_DIR; cd $RAW_DIR
echo -e "(`date`) there are total: `ls -1 | wc -l`  raw qseq data" | tee -a $LOG_FILE_STEP

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.2 starting decompressing " | tee -a $LOG_FILE
date| tee -a $LOG_FILE_STEP
ls -1 *.gz | parallel -j $TOTAL_PROC_NO gunzip |tee -a $LOG_FILE
wait;echo -e "(`date`) decompressed the raw qseq data" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.3 starting converting qseq to fastq " | tee -a $LOG_FILE

qseq2fastqPar (){
    count=1;ncore=$TOTAL_PROC_NO;
    prefix=$WORKING_DIR"/"
    while read data
    do

        if [ `echo $count" % "$ncore | bc` -eq 0 ] 
        then
	    qseq2fastq.pl $data ${data:0:10}'.fastq' 
        else
	    qseq2fastq.pl $data ${data:0:10}'.fastq' &
        fi
        count=$((count+1))
    done
}

date| tee -a $LOG_FILE_STEP
ls -1 * | qseq2fastqPar | tee -a $LOG_FILE_STEP
wait;echo -e "(`date`) 1.3 convering finished" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e " (`date`) 1.4 starting demultiplex " | tee -a $LOG_FILE


STEP="consolidate";
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

demultiplexFunMuticore (){
    count=1
    ncore=$TOTAL_PROC_NO;
    prefix=$WORKING_DIR"/"
    while read data
    do
        suffix=${data:5:11}
        echo 'processing:'$suffix
        idxfile=$(echo "$data" |sed -r 's/_1_/_2_/')
        echo "using idx file: "$idxfile
        echo "calulating no.$count sampling ......"
	echo "suffix is $suffix"
	echo "save into folder $prefix"
	echo "(`date`) $count"
        if [ `echo $count" % "$ncore | bc` -eq 0 ] #mod ncore
        then
	    cat $data | fastx_barcode_splitter.pl --bcfile $BARCODE_FILE --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1
       else
           cat $data | fastx_barcode_splitter.pl --bcfile $BARCODE_FILE --prefix $prefix --suffix $suffix --idxfile $idxfile --mismatches 1 &
        fi
        count=$((count+1))
    done
}

date | tee -a $LOG_FILE_STEP
ls -1 *_1_*.fastq | demultiplexFunMuticore 1|tee -a $LOG_FILE_STEP 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP

wait;echo -e "(`date`) 1.4 demultiplex finished" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.5 start merge fastq " | tee -a $LOG_FILE


cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    eval "cat ${prefix}_* > ${prefix}.fastq"  
    echo "(`date`)$prefix"
done<$BARCODE_FILE | tee -a $LOG_FILE_STEP 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.6  start remove  original fastq " | tee -a $LOG_FILE
mkdir -p $WORKING_DIR'/unassigned'
mv unmatched* $WORKING_DIR'/unassigned'
rm *_*.fastq | tee -a $LOG_FILE_STEP 
ls -lght >> $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Part 1.7: qc of consolidated FASTQ files" | tee -a $LOG_FILE

WORKING_DIR=$WORKING_DIR'/fastqc';mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
		      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {}     
                      1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
date | tee -a $LOG_ERR_FILE
wait;echo -e "(`date`) Step 1.7 QC consolidated fastq Finshed!" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2: trimming" | tee -a $LOG_FILE
WORKING_DIR=$PARENT_DIR'/02trim'
mkdir -p $WORKING_DIR
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
ls -1 *.fastq | xargs -n1 -P $TOTAL_PROC_NO -i \
                      cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  {} \
                      -o $WORKING_DIR"/"{}".trim.fastq" \
                      1>>$LOG_FILE_STEP 2>> $LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP
wait;echo -e "(`date`) Step 2: trimming  Finshed!"| tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2.1: QC of trimed QC files" | tee -a $LOG_FILE

cd $WORKING_DIR
WORKING_DIR=$WORKING_DIR'/fastqc';mkdir -p $WORKING_DIR

date | tee -a $LOG_FILE_STEP
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i\
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {}\
                      1>>$LOG_FILE_STEP 2>>$LOG_ERR_FILE
date | tee -a $LOG_ERR_FILE
wait;echo -e "(`date`)Step 2.1 QC trimed fastq Finshed!" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.1: alignment" | tee -a $LOG_FILE

WORKING_DIR=$PARENT_DIR'/03alignment/a_all'; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP

alignfun (){
    while read data;do
    samplename=${data%.*}
    echo -e "(`date`) align $samplename" 
    bowtie2 -t --non-deterministic --mm --phred33 -p $NPROC_PER_SAMPLE --very-sensitive -x mm10 -U $data | samtools sort -@ $NPROC_PER_SAMPLE -T $samplename -o $WORKING_DIR"/"$samplename".bam" &
done
}
ls -1 *.fastq | alignfun |tee -a $LOG_FILE_STEP 2>> $LOG_FILE_STEP


wait;echo -e "`date`: Step 3.1 alignment Finshed!" | tee -a $LOG_FILE


#------------------------------------------------------------
cd $WORKING_DIR
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.2: index bam file" | tee -a $LOG_FILE

date | tee -a $LOG_FILE_STEP
ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i\
                    samtools index {}\
                    1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP
wait;echo -e "(`date`) Step 3.2 index bam file Finshed!" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.1: filtering the aligned bam files" | tee -a $LOG_FILE
STEP='04filter'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i \
                     filterfun.sh {} $WORKING_DIR $NPROC_PER_SAMPLE \
    | tee -a $LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  Step 4.1 fiter bam finished" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE

cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --progress -j $SAMPLE_NO samtools index {} 
                    1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP

prefix=".txt"		
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done | tee -a $LOG_FILE_STEP

date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  finished Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.3a: deduplicate" | tee -a $LOG_FILE

STEP='04a_deDupped'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"

date | tee -a $LOG_FILE_STEP
mkdir -p $WORKING_DIR"/tmp"
#JAVA_OPTS=-Xmx64g
dedupfun (){
    while read data; do
    namestr=$(echo "$data"|cut -f1 -d".")
    java -Xmx64g -jar /opt/picard/picard.jar MarkDuplicates I=$data O=$WORKING_DIR"/"$namestr".bam" REMOVE_DUPLICATES=true ASSUME_SORTED=true M=$WORKING_DIR"/"$namestr".txt" TMP_DIR=$WORKING_DIR"/tmp"
done
}
ls -1 *.bam | dedupfun 1>>$LOG_FILE_STEP 2>>$LOG_FILE_STEP
rm -r $WORKING_DIR"/tmp"
date | tee -a $LOG_FILE_STEP

cd $WORKING_DIR
date | tee -a $LOG_FILE_STEP
ls *.bam | parallel --progress -j $SAMPLE_NO samtools index {}
1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP

prefix=".txt"
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done | tee -a $LOG_FILE_STEP

date | tee -a $LOG_FILE_STEP

echo -e "(`date`)  Step 4.3 deduplicate finished" | tee -a $LOG_FILE


echo -e "############################################################" | tee -a $LOG_FILE
echo -e " (`date`) Starting Step 4.3:  QC the bam" | tee -a $LOG_FILE

#java -Xmx16g -jar /opt/picard/picard.jar BamIndexStats I=./AlphaKO-0.bam
#bam_stat.py -i ./AlphaKO-0.bam -q 30


wait;echo -e "(`date`) Step 4.3 Finshed!" | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 4.4: spp cross-corelation " | tee -a $LOG_FILE
STEP='04b_spp'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
LOG_FILE_STEP=$WORKING_DIR"/log.txt"
LOG_ERR_FILE=$WORKING_DIR"/err.txt"

sppfun (){
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	Rscript /opt/phantompeakqualtools/run_spp.R -c="$data" -savp=$WORKING_DIR"/"$nameStr".pdf" -p=$NPROC_PER_SAMPLE &
    done
}

ls -1 *.bam | sppfun 1>> $LOG_FILE_STEP 2>> $LOG_ERR_FILE



STEP='qulimap'
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR;

qualmapfun (){
    int1=1;int2=1;
    while read data; do
	nameStr=$(echo "$data"| cut -f1 -d".")
	mkdir -p $WORKING_DIR'/'$nameStr

	if [ `echo $int1" % 4" | bc` -eq 0 ]
	then
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G
	else
	    echo -e "caculating $int1/$SAMPLE_NO samples \n"
	    int1=$((int1+int2))
	    JAVA_OPTS="-Djava.awt.headless=true" qualimap bamqc -nt 2 -c  -bam "$data"  -outdir $WORKING_DIR/$nameStr --java-mem-size=4G & 
	fi
    done
}



ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 5: peak calling " | tee -a $LOG_FILE

echo -e "############################################################" | tee -a $LOG_FILE







echo -e "(`date`)Starting Step 6: make tracks " | tee -a $LOG_FILE

STEP='06tracks'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;

LOG_FILE_STEP=$WORKING_DIR"/log.txt"
LOG_ERR_FILE=$WORKING_DIR"/err.txt"
date | tee -a $LOG_FILE_STEP
trackfun (){
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
	bamCoverage --bam "$data" --binSize 10 --normalizeTo1x 2150570000  --smoothLength 30 --extendReads 200 -p $NPROC_PER_SAMPLE -o $WORKING_DIR/$nameStr".bw" --ignoreForNormalization chrX &
    done
}
ls *.bam | trackfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE_STEP
date | tee -a $LOG_FILE_STEP

STEP='07multiqc'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;
multiqc . -o ./07multiqc &


# END







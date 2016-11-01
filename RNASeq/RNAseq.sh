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
echo -e "There are $SAMPLE_NO samples in this experiment \n" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."

#number  of processors per sample for fastqc 
NPROC_PER_SAMPLE=2
echo -e " $NPROC_PER_SAMPLE processors per sample \n" | tee -a $LOG_FILE
TOTAL_PROC_NO=$((SAMPLE_NO*NPROC_PER_SAMPLE)) # calculate number of total processor for the user
echo -e " total: $TOTAL_PROC_NO processors will be using for this analysis \n" | tee -a $LOG_FILE

#echo -e "Please input the password info (example: )\n" | tee -a $LOG_FILE
echo -e "Check the password file (password.txt)" | tee -a $LOG_FILE
PASSWD_INFO_INPUT=$(head -n 1 password.txt)
echo -e "Your password is : $PASSWD_INFO_INPUT"
#read -p"Press [Enter] key to continue..."

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) start running the pipelines " | tee -a $LOG_FILE
echo -e "(`date`) 1.1  starting downloading" | tee -a $LOG_FILE
echo -e "############################################################"| tee -a $LOG_FILE

STEP="00raw"; WORKING_DIR=$PARENT_DIR"/"$STEP
mkdir -p $WORKING_DIR; cd $WORKING_DIR
#read -p"Press [Enter] key to continue..."
grab_bscrc.sh $PASSWD_INFO_INPUT

wait;echo -e "(`date`) downloaded the raw qseq data" | tee -a $LOG_FILE
RAW_DIR=$(echo $PASSWD_INFO_INPUT|cut -f1 -d":"); RAW_DIR=$WORKING_DIR/$RAW_DIR; cd $RAW_DIR
echo -e "(`date`) there are total: `ls -1 | wc -l`  raw qseq data" | tee -a $LOG_FILE

echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.2 starting decompressing " | tee -a $LOG_FILE
ls -1 *.gz | parallel -j $TOTAL_PROC_NO  gunzip #--eta
wait;echo -e "(`date`) decompressed the raw qseq data" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."

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

ls -1 * | qseq2fastqPar | tee -a $LOG_FILE
wait;echo -e "(`date`) 1.3 convering finished" | tee -a $LOG_FILE

#read -p"Press [Enter] key to continue..."
echo -e "############################################################"| tee -a $LOG_FILE
echo -e " (`date`) 1.4 starting demultiplex " | tee -a $LOG_FILE


STEP="consolidate";
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR

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

ls -1 *_1_*.fastq | demultiplexFunMuticore 1|tee -a $LOG_FILE 2>>$LOG_ERR_FILE

wait;echo -e "(`date`) 1.4 demultiplex finished" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.5 start merge fastq " | tee -a $LOG_FILE


cd $WORKING_DIR
while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    eval "cat ${prefix}_* > ${prefix}.fastq  && rm ${prefix}_* "
    echo "(`date`)$prefix"
done<$BARCODE_FILE | tee -a $LOG_FILE


#read -p"Press [Enter] key to continue..."
echo -e "############################################################"| tee -a $LOG_FILE
echo -e "(`date`) 1.6  start remove unmatched fastq " | tee -a $LOG_FILE

ls -lght >> $LOG_FILE
mkdir -p $WORKING_DIR'/unmatched'
mv unmatched_*  $WORKING_DIR'/unmatched'

wait;echo -e "(`date`) 1.6 unmatched removed " | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Part 1.7: qc of consolidated FASTQ files" | tee -a $LOG_FILE
WORKING_DIR=$WORKING_DIR'/fastqc'
mkdir -p $WORKING_DIR
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait;echo -e "(`date`) Step 1.7 QC consolidated fastq Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
#------------------------------------------------------------
# 5. remove tmp files & mv files to the server 
#------------------------------------------------------------
# define serverDir
#serverDir="/mnt/biggie/backed_up/frank/projects/dynamic-decoding/caRNA/batch3/2000/"
#rm *.fastq # rm tmp fastq files 
#ls -1 *.txt | parallel -j 48 --eta gzip -9  # gzip txt  raw file

############################################################
# mapping
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2: trimming" | tee -a $LOG_FILE
WORKING_DIR=$PARENT_DIR'/02trim'
mkdir -p $WORKING_DIR


ls -1 *.fastq | xargs -n1 -P $TOTAL_PROC_NO -i \
                      cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  {} \
                      -o $WORKING_DIR"/"{}".trim.fastq" \
                      1>>$LOG_ERR_FILE 2>> $LOG_FILE        
wait;echo -e "(`date`) Step 2: trimming  Finshed!"| tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2.1: QC of trimed QC files" | tee -a $LOG_FILE

cd $WORKING_DIR
WORKING_DIR=$WORKING_DIR'/fastqc';mkdir -p $WORKING_DIR
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait;echo -e "(`date`)Step 2.1 QC trimed fastq Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."



echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.1: alignment" | tee -a $LOG_FILE

WORKING_DIR=$PARENT_DIR'/03alignment'; mkdir -p $WORKING_DIR;


# --readFilesCommand gunzip -c \
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      STAR --genomeDir /opt/ngs_indexes/star/mm10.primary_assembly.gencode.vM6_refchrom.50bp \
                      --runThreadN $NPROC_PER_SAMPLE --readFilesIn {} \
                      --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
                      --limitBAMsortRAM 17179869184 --outFilterType BySJout\
                      --outFilterMultimapNmax 20 --alignSJoverhangMin 8\
                      --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
                      --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20\
                      --alignIntronMax 1000000 --seedSearchStartLmax 30\
                      --outFileNamePrefix $WORKING_DIR'/'{}\
                      --genomeLoad LoadAndKeep \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
#| tee -a $LOG_FILE
wait;echo -e "`date`: Step 3.1 alignment Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."

#------------------------------------------------------------
cd $WORKING_DIR
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.2: index bam file" | tee -a $LOG_FILE

ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i \
                    samtools index {} \
                    1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "(`date`) Step 3.2 index bam file Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
#------------------------------------------------------------


echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.1: filtering the aligned bam files" | tee -a $LOG_FILE
STEP='04filter'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;


ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i \
                     filterfun.sh {} $WORKING_DIR $NPROC_PER_SAMPLE \
                     | tee -a $LOG_FILE

echo -e "(`date`)  Step 4.1 fiter bam finished" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE


cd $WORKING_DIR
ls *.bam | parallel -j $SAMPLE_NO samtools index {} 
                    1>>$LOG_ERR_FILE 2>>$LOG_FILE

prefix=".txt"		
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done
echo -e "(`date`)  finished Step 4.2: indexing the filtered bam" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
echo -e "############################################################" | tee -a $LOG_FILE
echo -e " (`date`) Starting Step 4.3: qualimap QC the mapping" | tee -a $LOG_FILE

STEP='a_qulimap'
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir -p $WORKING_DIR;

qualmapfun (){
    int1=1;int2=1;
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        mkdir -p $WORKING_DIR'/'$nameStr
        
        if [ `echo $int1" % 5" | bc` -eq 0 ]
        then
            echo -e "caculating $int1/$SAMPLE_NO samples \n"
            int1=$((int1+int2))
            JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $WORKING_DIR/$nameStr --java-mem-size=4G            
        else
            echo -e "caculating $int1/$SAMPLE_NO samples \n"
            int1=$((int1+int2))
            JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $WORKING_DIR/$nameStr --java-mem-size=4G &
        fi
    done
}



ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "(`date`) Step 4.3 Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 5: generate the counts file " | tee -a $LOG_FILE


STEP='05counts'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;


featureCounts -T $TOTAL_PROC_NO -s 2 -t exon -g gene_id -a /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -o $WORKING_DIR/counts-gene.txt *.bam 1>>$LOG_ERR_FILE 2>>$LOG_FILE

wait;echo -e "(`date`) Step 5 Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."
echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 6: make tracks " | tee -a $LOG_FILE

STEP='06tacks'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir -p $WORKING_DIR;

trackfun (){
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        bam2wig.py -i "$data" -s /opt/ngs_indexes/genomes/mm/mm10.chrom.sizes -t 500000000 -d '+-,-+' -o $WORKING_DIR/$nameStr".bw" & 
    done
}
ls *.bam | trackfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE

echo -e "(`date`)finished tracks & now deleting wig files " | tee -a $LOG_FILE

cd $WORKING_DIR; rm *.wig 
wait;echo -e "(`date`) Step 6 Finshed!" | tee -a $LOG_FILE
#read -p"Press [Enter] key to continue..."

#echo -e "############################################################" | tee -a $LOG_FILE
#echo -e "(`date`)Starting Step 7: compress all fastq files" | tee -a $LOG_FILE


#cd $PARENT_DIR
#find . -type f -name "*.fastq" | parallel -j $TOTAL_PROC_NO --eta gzip -9 | tee -a $LOG_FILE
#wait;echo -e "(`date`) Step 7 Finshed!" | tee -a $LOG_FILE
#read -p "Press [Enter] key to continue..."

echo -e "############################################################" | tee -a $LOG_FILE
echo -e "(`date`) Congratuations! U finished the pipeline in total " | tee -a $LOG_FILE



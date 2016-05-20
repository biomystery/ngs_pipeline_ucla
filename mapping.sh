#!/bin/bash
# The Script for mapping
############################################################
# input parameters
############################################################
#input 1: consolidated FASTQ Folder 
echo -n "1/2. Enter fastq folder \n" 
read FASTQ_DIR # /home/frank/caRNA/2000/01_consolidated
echo -e "input is $FASTQ_DIR \n"

cd $FASTQ_DIR;cd ../../; PARENT_DIR=$PWD; cd $FASTQ_DIR;
# log files
LOG_FILE=$PARENT_DIR"/run.log.txt"; LOG_ERR_FILE=$PARENT_DIR"/run.err.txt"

SAMPLE_NO=`ls -1 *.fastq|wc -l` 
echo -e "There are $SAMPLE_NO samples in the folder ($FASTQ_DIR) \n" | tee -a $LOG_FILE

#input 2: number of processors per sample for fastqc 
echo -e "2/2. Enter the number of processor per sample \n" 
read NPROC_PER_SAMPLE
echo -e "input is $NPROC_PER_SAMPLE\n"
let TOTAL_PROC_NO=$SAMPLE_NO*$NPROC_PER_SAMPLE # calculate number of total processor for the user



#------------------------------------------------------------
#  1.1 QC of sequencing reads
#------------------------------------------------------------
WORKING_DIR=$FASTQ_DIR'/fastqc';mkdir $WORKING_DIR

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 1.1: QC of consolidated FASTQ files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
echo -e "Step 1.1 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
# 2. FASTQ trimming
#------------------------------------------------------------
WORKING_DIR=$PARENT_DIR'/02trim'; mkdir $WORKING_DIR;

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 2: trimming\n" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
ls -1 *.fastq | xargs -n1 -P $PROCESSORS_NO -i \
                      cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  {} \
                      -o $WORKING_DIR{}".trim.fastq" \
                      1>>$LOG_ERR_FILE 2>> $LOG_FILE        
wait;echo -e "Step 2 Finshed!"| tee -a $LOG_FILE

# 1.1 qc of trimed fastqc
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 2.1: QC of trimed QC files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

cd $WORKING_DIR;WORKING_DIR=$WORKING_DIR'/fastqc';mkdir $WORKING_DIR
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait;echo -e "Step 2.1 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
#3. Mapping/alignment. 
#------------------------------------------------------------
# 3.1 compress the trimed fastqc
WORKING_DIR=$PARENT_DIR'/03alignment'; mkdir $WORKING_DIR;
	--readFilesCommand gunzip -c & 

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 3, total xx steps:\n" | tee -a $LOG_FILE
echo -e "Starting Step 3.1: alignment\n" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

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
                       1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "Step 3.1 Finshed!" | tee -a $LOG_FILE

# 3.2 compress the trimed fastqc
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 3.2: compress the trimed QC files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

ls -1 *.fastq | parallel -j $TOTAL_PROC_NO --eta gzip -9 | tee -a $LOG_FILE
wait;echo -e "Step 3.2 Finshed!" | tee -a $LOG_FILE


#------------------------------------------------------------
# 4. index the bam file 
#------------------------------------------------------------
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 4: index sam file" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i \
                    samtools index {} \
                    1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "Step 4 Finshed!" | tee -a $LOG_FILE

#5.  Filtering (multiple maping reads)
mkdir $homeDir$step5
cd $homeDir$step4
prefix=".filtered"
filetype='.bam'
filterfun (){
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        #echo $nameStr
        samtools view -F 2820 -q 30 -@ 2 -b $data > $homeDir$step5$nameStr$prefix$filetype & 
    done
}
ls *.bam | filterfun


cd $homeDir$step5
ls *.bam |while read data; do samtools index "$data"  & done

prefix=".txt"		
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done

#6.  Mapped read QC
mkdir $homeDir$step6
cd $homeDir$step5
qualmapfun (){
    int1=1;int2=3;
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        mkdir $homeDir$step6$nameStr
        
        if [ `echo $int1" % 2" | bc` -eq 0 ]
        then
            int1=$((int1+int2))
            JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $homeDir$step6$nameStr --java-mem-size=4G            
        else
            int1=$((int1+int2))
            JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $homeDir$step6$nameStr --java-mem-size=4G &
        fi
    done
}

#test if


    

ls *.bam | qualmapfun
#${prefix}"${data:9:10}"




#7.  Counting
mkdir $homeDir$step7
cd $homeDir$step5 #filtered 
featureCounts -T 44 -s 2 -t exon -g gene_id -a /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -o ../07count/counts.txt *.bam &

# 8. Tracks
cd $homeDir$step5
mkdir $homeDir$step8
trackfun (){
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        bam2wig.py -i "$data" -s /opt/ngs_indexes/genomes/mm/mm10.chrom.sizes -o track.bw -t 500000000 -d '+-,-+' -o $homeDir$step8$nameStr & 
    done
}
ls *.bam | trackfun





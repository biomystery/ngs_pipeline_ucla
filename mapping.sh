#!/bin/bash
# The Script for mapping
echo -e "############################################################"
echo -e "(`date`) initating input parameters ....."
echo -e "############################################################"

#input 1: consolidated FASTQ Folder 
echo -e "(`date`)Setting folder \n" 
#read FASTQ_DIR # /home/frank/caRNA/2000/01_consolidated
FASTQ_DIR=$PWD # /home/frank/caRNA/2000/01_consolidated
echo -e "(`date`) raw fastq folder is $FASTQ_DIR \n"
cd $FASTQ_DIR;cd ../; PARENT_DIR=$PWD; cd $FASTQ_DIR;
echo -e "(`date`) project folder is $PARENT_DIR \n"

# log files
LOG_FILE=$PARENT_DIR"/run.log.txt"; LOG_ERR_FILE=$PARENT_DIR"/run.err.txt"

SAMPLE_NO=`ls -1 *.fastq|wc -l` 
echo -e "There are $SAMPLE_NO samples in the folder ($FASTQ_DIR) \n" | tee -a $LOG_FILE

#input 2: number of processors per sample for fastqc 
#echo -e " Enter the number of processor per sample \n" 
#read NPROC_PER_SAMPLE
NPROC_PER_SAMPLE=2
echo -e " $NPROC_PER_SAMPLE process per sample \n" | tee -a $LOG_FILE
let TOTAL_PROC_NO=$SAMPLE_NO*$NPROC_PER_SAMPLE # calculate number of total processor for the user


echo -e "############################################################"
echo -e "`date` Running the pipelines \n" 
echo -e "############################################################"

#------------------------------------------------------------
#  1 QC of the sequencing reads
#------------------------------------------------------------
WORKING_DIR=$FASTQ_DIR'/fastqc';mkdir $WORKING_DIR

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "Starting Step 1.1: QC of consolidated FASTQ files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
echo -e "(`date`) Step 1.1 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
# 2. FASTQ trimming
#------------------------------------------------------------
WORKING_DIR=$PARENT_DIR'/02trim'; mkdir $WORKING_DIR;

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2: trimming\n" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
ls -1 *.fastq | xargs -n1 -P $PROCESSORS_NO -i \
                      cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  {} \
                      -o $WORKING_DIR{}".trim.fastq" \
                      1>>$LOG_ERR_FILE 2>> $LOG_FILE        
wait;echo -e "(`date`) Step 2 Finshed!"| tee -a $LOG_FILE

# 1.1 qc of trimed fastqc
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 2.1: QC of trimed QC files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

cd $WORKING_DIR;WORKING_DIR=$WORKING_DIR'/fastqc';mkdir $WORKING_DIR
ls -1 *.fastq | xargs -n1 -P $SAMPLE_NO -i \
                      fastqc -t $NPROC_PER_SAMPLE -outdir $WORKING_DIR {} \
                      1>>$LOG_FILE 2>>$LOG_ERR_FILE
wait;echo -e "(`date`)Step 2.1 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
#3. Mapping/alignment. 
#------------------------------------------------------------
# 3.1 compress the trimed fastqc
WORKING_DIR=$PARENT_DIR'/03alignment'; mkdir $WORKING_DIR;
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3, total xx steps:\n" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.1: alignment\n" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

# --readFilesCommand gunzip -c \
ls -1 *.gz | xargs -n1 -P $SAMPLE_NO -i \
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
wait;echo -e "`date`: Step 3.1 Finshed!" | tee -a $LOG_FILE

# 3.2 compress the trimed fastqc
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 3.2: compress the trimed QC files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

ls -1 *.fastq | parallel -j $TOTAL_PROC_NO --eta gzip -9 | tee -a $LOG_FILE
wait;echo -e "(`date`) Step 3.2 Finshed!" | tee -a $LOG_FILE


#------------------------------------------------------------
# 3.1. index the bam file 
#------------------------------------------------------------
cd $WORKING_DIR
echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 3.3: index sam file" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
ls -1 *.bam | xargs -n1 -P $SAMPLE_NO -i \
                    samtools index {} \
                    1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "(`date`) Step 3.3 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
#4.  Filtering (multiple maping reads)
#------------------------------------------------------------
STEP='04filter'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir $WORKING_DIR;


filterfun (){
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
}

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`) Starting Step 4: filtering the aligned bam files" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

export -f filterfun
ls *.bam |parallel --progress -j $SAMPLE_NO  filterfun {} $WORKING_DIR $NPROC_PER_SAMPLE | tee -a $LOG_ERR_FILE 


cd $WORKING_DIR
ls *.bam | parallel --progress -j $SAMPLE_NO samtools index {} 
                    1>>$LOG_ERR_FILE 2>>$LOG_FILE

prefix=".txt"		
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done

#------------------------------------------------------------
#4.1. QC the mapping 
#------------------------------------------------------------
STEP='a_qulimap'
WORKING_DIR=$WORKING_DIR'/'$STEP; mkdir $WORKING_DIR;

qualmapfun (){
    int1=1;int2=3;
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        mkdir $WORKING_DIR'/'$nameStr
        
        if [ `echo $int1" % 2" | bc` -eq 0 ]
        then
            echo -e "caculating $int1/$SAMPLE_NO samples \n"
            int1=$((int1+int2))
            JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $WORKING_DIR$nameStr --java-mem-size=4G            
        else
            echo -e "caculating $int1/$SAMPLE_NO samples \n"
            int1=$((int1+int2))
            JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $WORKING_DIR$nameStr --java-mem-size=4G &
        fi
    done
}



qualmapfun2 (){
    data=$1
    nameStr=$(echo "$data"| cut -f1 -d".")
    homeDir=WORKING_DIR
    mkdir $homeDir$nameStr
    JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir $homeDir$nameStr --java-mem-size=4G            
}

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e " (`date`) Starting Step 4.1: qualimap QC the mapping" | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

ls *.bam | qualmapfun 1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "(`date`) Step 4.1 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
#5.  generate counts 
#------------------------------------------------------------
STEP='05counts'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir $WORKING_DIR;

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 5: generate the counts file " | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE

featureCounts -T $TOTAL_PROC_NO -s 2 -t exon -g gene_id -a /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -o $WORKING_DIR/counts-gene.txt *.bam | tee -a $LOG_FILE
#1>>$LOG_ERR_FILE 2>>$LOG_FILE
wait;echo -e "(`date`) Step 5 Finshed!" | tee -a $LOG_FILE

#------------------------------------------------------------
# 6. Tracks
#------------------------------------------------------------
STEP='06tacks'
WORKING_DIR=$PARENT_DIR'/'$STEP; mkdir $WORKING_DIR;

echo -e "--------------------\n" | tee -a $LOG_FILE
echo -e "(`date`)Starting Step 6: make tracks " | tee -a $LOG_FILE
echo -e "--------------------\n" | tee -a $LOG_FILE
trackfun (){
    while read data; do
        nameStr=$(echo "$data"| cut -f1 -d".")
        bam2wig.py -i "$data" -s /opt/ngs_indexes/genomes/mm/mm10.chrom.sizes -o track.bw -t 500000000 -d '+-,-+' -o $homeDir$step8$nameStr & 
    done
}
ls *.bam | trackfun





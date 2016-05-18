# the full pipeline for caRNA seq analysis of the 2000 
rsync -a zhangcheng@164.67.9.46:~/xx_barcode ./

# define paths & steps here 
#proj_path="/mnt/biggie/no_backup/frank/kim/hiseq4000/"
fastq_path="/home/frank/caRNA/2000/"
step0="00raw/"
step1="01fastqc/"
step2="02trim/"
step_back_fastq="00fasq"
step3="03fastqc/"
step4="04alignment/"
step5="05filtered/"
step6="06qualimap/"
step7="07count/"
step8="08tracks/"
homeDir=$HOME'/'
barcodefile=$fastq_path"xx_barcode/barcode.txt"

############################################################
# Download & convert
############################################################
mkdir cd $fastq_path$step0; cd $fastq_path$step0
grab_bscrc  SxaQSEQsVA107L7:KN6ga9yV3cG6

#ls -1 s_7_1* | while read data; do echo ${data:6:4}; done  > name1.txt
#ls -1 s_7_2* | while read data; do echo ${data:6:4}; done  > name2.txt
#diff name1.txt name2.txt
# with anonymous pipe 
#diff <(ls -1 s_7_1* | while read data; do echo ${data:6:4}; done) <(ls -1 s_7_2* | while read data; do echo ${data:6:4}; done)
#diff <(ls -1 *.txt | while read data; do echo ${data:6:10}; done) <(ls -1 *.fastq | while read data; do echo ${data:6:4}; done)
#diff <(ls -1 *.f | while read data; do echo ${data:6:10}; done) <(ls -1 *.fastq | while read data; do echo ${data:6:4}; done)

diff <(ls -1 *_1_*.fastq | while read data; do echo ${data:6:4}; done) <(ls -1 *_2_*.fastq | while read data; do echo ${data:6:4}; done)

ls -1 * | while read data; do  qseq2fastq.pl $data ${data:0:10}'.fastq' & done

qseq2fastqfunc (){
    while read data
    do
        qseq2fastq.pl $data ${data:0:10}'.fastq'
    done
}

ls -1 *.txt | parallel -j 48 --eta qseq2fastqfunc

# Demultiplex FASTQ

mkdir ../02_consolidated

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
#ls -1 *_1_*.fastq | parallel -j 48 --eta demultiplexFun 

demultiplexFunMuticore (){
    count=1;ncore=48;
    while read data
    do
        suffix=${data:5:11}
        echo 'processing:'$suffix
        idxfile=$(echo "$data" |sed -r 's/_1_/_2_/')
        echo "using idx file: "$idxfile
        echo "calulating no.$count sampling ......"
        if [ `echo $count" % "$ncore | bc` -eq 0 ] #mod ncore
        then
            cat $data | fastx_barcode_splitter.pl --bcfile $barcodefile --prefix ../02_consolidated/ --suffix $suffix --idxfile $idxfile --mismatches 1
        else
            cat $data | fastx_barcode_splitter.pl --bcfile $barcodefile --prefix ../02_consolidated/ --suffix $suffix --idxfile $idxfile --mismatches 1 &
        fi
        count=$((count+1))
    done
}

ls -1 *_1_*.fastq | demultiplexFunMuticore

# check if all files demultiplexed
diff <(ls -1 ../02_consolidated/ifnar-0*.fastq |  while read data; do echo ${data: -11}; done) <(ls -1 *_2_*.fastq | while read data; do echo ${data:5:11}; done)

# head -c 2 # check the first two characters 

# 4. Merge FASTQ
cd ../02_consolidated/
while read line; do
    prefix=$(echo -e "$line" | cut -f1 -d$'\t')
    #eval "ls ${prefix}* "
    #eval "cat ${prefix}* > ${prefix}.fastq && rm ${prefix}_* &"
    eval "rm ${prefix}_* " 
done<$barcodefile

# cat alpha-0_*.fastq >alpha-0.fastq

# 5. remove tmp files 

#------------------------------------------------------------
#  1. QC of sequencing reads
#------------------------------------------------------------

echo  ${proj_path}${step1_fastqc}
ls  ${proj_path}${step1_fastqc}
ls  ${fastq_path}
mkdir -p $proj_path

mkdir $proj_path'01fastqc'

qcfun (){
    while read data; do
        fastqc -t 2 -outdir $proj_path'01fastqc' $data &
    done
}

cd $fastq_path
ls *.fastq | qcfun


#find /home/shared/dropbox/supriya/caRNA/batch3/ -name '*.fastq'|  xargs -n 1 fastqc -t 4 -outdir /mnt/biggie/no_backup/frank/caRNA/batch3/01fastqc &

# 2.  trimming

mkdir $proj_path$step2
trimfun () {
    while read data; do 
        cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  $data -o $proj_path'02trim/'$data".trim.fastq" &
    done
}

ls | trimfun &

mkdir $proj_path$step_back_fastq
cd $fastq_path
ls | while read data; do sudo gzip -9 $data &  done
sudo mv ${fastq_path}'*' $proj_path$step_back_fastq & 

# 3. qc of trimed fastqc
mkdir $proj_path$step3
mkdir $homeDir$step3

#ls | fastqc  -t 4 -outdir /mnt/biggie/no_backup/frank/caRNA/batch3/03fastqc 
qcfun (){
    while read data; do
        fastqc  -t 2 -outdir $homeDir$step3 $data &
    done
}
cd $proj_path$step2
ls $proj_path$step2 | qcfun
ls *.fastq | qcfun

		
#	4. Mapping/alignment. 
mkdir $proj_path$step4

mkdir $HOME'/'$step4

starfun (){
    while read data; do
	STAR --genomeDir /opt/ngs_indexes/star/mm10.primary_assembly.gencode.vM6_refchrom.50bp --runThreadN 2 --readFilesIn $data --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 17179869184 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --seedSearchStartLmax 30 --outFileNamePrefix $HOME'/'$step4$data --genomeLoad LoadAndKeep & #--readFilesCommand gunzip -c &
    done
}

cd $proj_path$step2
sudo ls *.fastq | starfun 

cd $proj_path$step2
sudo ls *.fastq | while read data; do sudo gzip -9 $data &  done

#find -name '*.fastq'|  xargs -n 1 gzip -9 &
#find -name '*.fastq'|  while read data;do  gzip -9 "$data" &  done 
#gzip -9  sample.fastq &
#gzip -9  trimmed.fastq &
		
# samtools index Aligned.sortedByCoord.out.bam
ls *.bam |while read data; do samtools index "$data"  & done

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





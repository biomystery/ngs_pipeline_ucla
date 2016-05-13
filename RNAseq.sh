#This applies to 
	
#  1. QC of sequencing reads
mkdir /mnt/biggie/no_backup/frank/caRNA/batch3/01fastqc

qcfun (){
    while read data; do
        fastqc -t 4 -outdir /mnt/biggie/no_backup/frank/caRNA/batch3/01fastqc &
    done
}

ls /home/shared/dropbox/supriya/caRNA/batch3/*.fastq | qcfun


#find /home/shared/dropbox/supriya/caRNA/batch3/ -name '*.fastq'|  xargs -n 1 fastqc -t 4 -outdir /mnt/biggie/no_backup/frank/caRNA/batch3/01fastqc &

# 2.  Read trimming

mkdir /mnt/biggie/no_backup/frank/caRNA/batch3/02trim

trimfun () {
    local STR=".trim.fastq"
    local DIR="/mnt/biggie/no_backup/frank/caRNA/batch3/02trim/"
    cutadapt -f fastq -e 0.1 -O 6 -q 20 -m 35 -a AGATCGGAAGAGC  $1 -o $DIR$1$STR 
}

wrapfun (){
    while read data; do 
	trimfun "$data" &
    done
}

cd /home/shared/dropbox/supriya/caRNA/batch3/
ls  | wrapfun 


# 3. qc of trimed fastqc
mkdir /mnt/biggie/no_backup/frank/caRNA/batch3/03fastqc
ls | fastqc  -t 4 -outdir /mnt/biggie/no_backup/frank/caRNA/batch3/03fastqc 
qcfun (){
    while read data; do
        fastqc  -t 4 -outdir /mnt/biggie/no_backup/frank/caRNA/batch3/03fastqc $data &
    done
}

ls /mnt/biggie/no_backup/frank/caRNA/batch3/02trim/*.fastq | qcfun

		
#	4. Mapping
mkdir /mnt/biggie/no_backup/frank/caRNA/batch3/04align


s1="../04align/"
starfun (){
    while read data; do
	STAR --genomeDir /opt/ngs_indexes/star/mm10.primary_assembly.gencode.vM6_refchrom.50bp --runThreadN 2 --readFilesIn $data --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate -- limitBAMsortRAM 17179869184 --outFilterType BySJout --outFilterMultimapNmax 20 -- alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 -- outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 -- seedSearchStartLmax 30 --outFileNamePrefix ${s1}"$data" --genomeLoad LoadAndKeep --readFilesCommand gunzip -c &
    done
}

        echo  ${s1}"$data"



cd /mnt/biggie/no_backup/frank/caRNA/batch3/02trim/
ls /mnt/biggie/no_backup/frank/caRNA/batch3/02trim/*.fastq | starfun &




find -name '*.fastq'|  xargs -n 1 gzip -9 &
find -name '*.fastq'|  while read data;do  gzip -9 "$data" &  done 
gzip -9  sample.fastq &
gzip -9  trimmed.fastq &
		
samtools index Aligned.sortedByCoord.out.bam
		
		
#	5.  Filtering (multiple maping reads)
mkdir ../05filtered
prefix="../05filtered/filtered."
ls *.bam |while read data; do samtools view -F 2820 -q 30 -@ 2 -b "$data" > ${prefix}"$data" & done

cd ../05filtered
ls *.bam |while read data; do samtools index "$data"  & done

		
prefix=".txt"		
#	6.  Mapped read QC
ls *.bam |while read data; do echo "$data"${prefix}  ; done
ls *.bam |while read data; do samtools flagstat "$data" > "$data"${prefix}  & done


mkdir ../06qualimap

prefix="../06qualimap"		
ls *.bam | while read data; do   mkdir ${prefix}"${data:9:10}"; done

ls *.bam | while read data; do JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam "$data" -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir ${prefix}"${data:9:10}" --java-mem-size=4G & done

JAVA_OPTS="-Djava.awt.headless=true" qualimap rnaseq -bam .bam -gtf /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -p strand-specific-reverse -outdir qualimap/report --java-mem-size=4G 

#7.  Counting
prefix="../07count/"
s=".count.txt"
mkdir ${prefix}

ls *.bam | while read data; do featureCounts -T 2 -s 2 -t exon -g gene_id -a /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -o ${prefix}"${data:9:10}"${s} "$data" & done

featureCounts -T 44 -s 2 -t exon -g gene_id -a /opt/ngs_indexes/models/mm/mm10/gencode.vM6.refchrom.annotation.gtf -o ../07count/counts.txt *.bam &

# 8. Tracks
mkdir ${prefix}
prefix="../08tracks/"
ls *.bam | while read data; do bam2wig.py -i "$data" -s /opt/ngs_indexes/genomes/mm/mm10.chrom.sizes -o track.bw -t 500000000 -d '+-,-+' -o ${prefix}"${data:9:10}" & done





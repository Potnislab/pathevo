########################################################################################
##############          making bam formated sample files              ###################
###########################################################################################

# Data preparation

# Raw reads cleaning

#!/bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bbmap/37.36
#input fastq files

for fq1 in /path_blklisted_genes_removed_samples/*_blk_filter.1.fq
do
echo "working with file $fq1"

base=$(basename $fq1 _1.fastq.gz)
echo "base name is $base"

fq1=/path_blklisted_genes_removed_samples/${base}_blk_filter.1.fq
fq2=/path_blklisted_genes_removed_samples/${base}_blk_filter.2.fq


bbduk.sh -Xmx1g in1=$fq1 in2=$fq2 out1=${base}.clean1.fq out2=${base}.clean2.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

bbduk.sh -Xmx1g in1=${base}.clean1.fq in2=${base}.clean2.fq out1=${base}.clean_trim1.fq out2=${base}.clean_trim2.fq qtrim=rl trimq=10


done

## 1.  Blacklisted_reads_removed Samples were mapping to reference genome using bwa-mem 

### On terminal window

module load bwa/0.7.12  ##load module

bwa index  Xp_AL65_AL22.fna  #Xp_AL65_AL22.fna is non-redudant pangenome made with Superpang


## indexing of reference file used AL65 and AL22 non-redundant pangenome
 
##in queue

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12

#input fastq files
for fq1 in /path/to/clean_reads/*.blk_filter.1.fq   ## made by the script mentioned in the https://github.com/Potnislab/AtDep_2021_metagenome/blob/main/Raw%20read%20processing/quality-trimming_script
do
echo "working with file $fq1"
base=$(basename $fq1 .blk_filter.1.fq)
echo "base name is $base"
fq1=/path/to/1/${base}.blk_filter.1.fq
fq2=/path/to/2/${base}.blk_filter.2.fq
bwa mem /path/to/indexed_reference/Xp_AL65_AL22.fna  $fq1 $fq2 > ${base}.sam    
done

## 2. Picard Sortsam  (https://gatk.broadinstitute.org/hc/en-us/articles/360046788792-SortSam-Picard-#--VALIDATION_STRINGENCY )
  
#!/bin/bash 
 
module load picard/1.79 
 
 
for file in /path/to/*.sam;    ## using the output files from previous step as input
    do tag=${file%.sam}; 
 
java -Xmx2g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar I="$file" O=$tag.sort.sam SORT_ORDER=coordinate 
 
done 

## O= output_file_name.sort.sam 



## 3. remove low quality alignments 
 
#!/bin/bash
module load samtools/1.11 
 
for file in /path/to/*.sort.sam;   ##using output of sort sam as input in this step
    do tag=${file%.sort.sam}; 
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -Sb "$file" > $tag.q20.bam 
done 

## 4.Picard remove duplicates

#!/bin/bash

module load picard/2.24.0


for file in *.LH3.q20.bam;
    do tag=${file%.LH3.q20.bam};

picard MarkDuplicates -I "$file" -O $tag.LH3.rmd.sort.bam -M $tag.dupstat.txt --REMOVE_DUPLICATES true --REMOVE_SEQUENCING_DUPLICATES true

done


##################################
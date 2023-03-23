# METAPOP:  https://github.com/metaGmetapop/metapop                 

## A. making bam formated sample files #####

## 1.  Read mapping to reference genome using bwa-mem 

On terminal window

        module load bwa/0.7.12  ##load module
        bwa index Xanthomonas_perforans_LH3.fna  ## indexing of reference file
 
run in the queue
Script for quality trimmed reads are available at https://github.com/Potnislab/AtDep_2021_metagenome/blob/main/Raw%20read%20processing/quality-trimming_script

        #!/bin/bash
        source /opt/asn/etc/asn-bash-profiles-special/modules.sh
        module load bwa/0.7.12
        
        #input fastq files
        for fq1 in /path/to/clean_reads/*.paired.1.fq   
        do
        echo "working with file $fq1"
        base=$(basename $fq1 .paired.1.fq)
        echo "base name is $base"
        fq1=/path/to/1/${base}.paired.1.fq
        fq2=/path/to/2/${base}.paired.2.fq
        bwa mem /path/to/indexed_reference/Xanthomonas_perforans_LH3.fna $fq1 $fq2 > ${base}.sam    ##.LH3.sam is output file
        done

## 2. Picard Sortsam 
(https://gatk.broadinstitute.org/hc/en-us/articles/360046788792-SortSam-Picard-#--VALIDATION_STRINGENCY )
using the output files from previous step as input
  
       #!/bin/bash
       module load picard/1.79
       for file in /path/to/*.sam;
       do tag=${file%.sam};
       java -Xmx2g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar I="$file" O=$tag.sort.sam SORT_ORDER=coordinate
       done 


O= output_file_name.sort.sam 

## 3. remove low quality alignments
using output of sort sam as input in this step
 
       #!/bin/bash
       module load samtools/1.11
       for file in /path/to/*.sort.sam;   
       do tag=${file%.sort.sam};
       samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -Sb "$file" > $tag.q20.bam
       done 

## B. Install METAPOP as 

        conda create --name metapop2 python=3.7

        #To activate this environment, use
        $ conda activate metapop2
        #To deactivate an active environment, use
        $ conda deactivate

## install dependencies 

        conda update --all
        conda install mamba=0.22.1
        mamba install -y bcftools samtools prodigal numpy pysam r-ggrepel r-base r-data.table r-ggplot2 r-rcolorbrewer r-doparallel r-cowplot r-bit64 r-gggenes r-stringr r-vegan r-compositions r-pheatmap -c bioconda -c conda-forge -c r
        
# pip installation of metapop        
        pip install metapop
        
# Input required: read count tsv file
make a tsv file containing all the samples with number of reads counts using
        
        samtools view -c SAMPLE.bam

or counting only mapped (primary aligned) reads
         
        samtools view -c -F 260 SAMPLE.bam

options 
-c  count reads and print the total number
-F 260  output primary aligned mapped reads read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion bit 3 + bit 9 = 4 + 256 = 260

tsv file format will be specified as columns
'BAM File Names'     'Number Reads or Bps in Read Library'

move  all the *.q20.bam in one folder 

         mv /path/to/q20.bam /path/to/Sample_.q20.bam


### Command

         #!/bin/bash
         metapop --input_samples /path/to/bam/directory/Sample_.q20.bam --reference /path/to/reference_directory \
         --norm /path/to/read_counts.tsv  --output /path/to/output_directory \
         --plot_all --snp_scale both

options
--plot_all FLAG : Metapop will print all contigs from microdiversity results. This will likely take a long time. Default prints top 3 genomes by highest % of genes under positive selection in each sample.
--snp_scale [local, global, both] : Metapop will print microdiversity results using SNPs assessed at the local (per sample) or global (across all samples) levels, or both. Defaults to local.

output files contain the Nucleotid diversity and values for selection pressure for different samples
it would be located in /path/to/output_directory/10.Microdiversity


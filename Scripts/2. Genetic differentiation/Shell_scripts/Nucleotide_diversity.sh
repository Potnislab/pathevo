#########################################################################################
##############                               MetaPop                   ###################
###########################################################################################

## METAPOP:  https://github.com/metaGmetapop/metapop 

##1. Install METAPOP as 

conda create --name metapop2 python=3.7

###To activate this environment, use
conda activate metapop2

# To deactivate an active environment, use
#
#     $ conda deactivate

##install dependencies 

conda update --all
conda install mamba=0.22.1

mamba install -y bcftools samtools prodigal numpy pysam r-ggrepel r-base r-data.table r-ggplot2 r-rcolorbrewer r-doparallel r-cowplot r-bit64 r-gggenes r-stringr r-vegan r-compositions r-pheatmap -c bioconda -c conda-forge -c r

pip install metapop


####################
## make a tsv file containing all the samples with number of reads counts using
samtools view -c SAMPLE.bam

## or counting only mapped (primary aligned) reads
samtools view -c -F 260 SAMPLE.bam

## options
## -c  count reads and print the total number
## -F 260  output primary aligned mapped reads read unmapped & not primary alignment criteria 3 & 9 are selected for exclusion bit 3 + bit 9 = 4 + 256 = 260

## tsv file format will be specified as columns
## 'BAM File Names'     'Number Reads or Bps in Read Library'


## move  all the *.q20.bam in one folder 

mv /path/to/q20.bam /path/to/Sample_.q20.bam


#we used the samples after removing blacklisted genes here.

###Command
#!/bin/bash
metapop --input_samples /path/to/bam/directory_including_samples_with_removed_blk_genes/*.q20.bam --reference /path/to/reference_directory/Xp_AL65_AL22.fna \
--norm /path/to/read_counts.tsv  --output /path/to/output_directory \
--plot_all --snp_scale both


##options
## --plot_all FLAG : Metapop will print all contigs from microdiversity results. This will likely take a long time. Default prints top 3 genomes by highest % of genes under positive selection in each sample.
## --snp_scale [local, global, both] : Metapop will print microdiversity results using SNPs assessed at the local (per sample) or global (across all samples) levels, or both. Defaults to local.



## output files contain the Nucleotid diversity and values for selection pressure for different samples

## it is located in /path/to/output_directory/Microdiversity

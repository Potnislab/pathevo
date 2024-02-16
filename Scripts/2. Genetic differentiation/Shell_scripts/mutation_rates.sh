##############################################
#############   MUTATION RATES         #######
##########                             #######
##############################################

##1. download rhometa repository using git clone command

git clone https://github.com/sid-krish/rhometa.git

## dependancies installation

conda create -name Nextflow

conda activate Nextflow

##1.installation using

conda install -c bioconda nextflow

## OR

conda install -c "bioconda/label/cf201901" nextflow 


##2 installed freebayes

conda install -c "bioconda/label/cf201901" freebayes


###### THETA ESTIMATES: MUTATION RATES #################

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/3-2021.11
module load samtools
module load bcftools/1.13

for file in /path/to/*.q20.bam;    ## samples were without blacklisted genes
do tag=${file%.q20.bam};

nextflow run theta_est.nf --bam "$file" --fa /path/to/Xp_AL22_AL65.fna --output_dir Theta_mutation_rates_output;
done

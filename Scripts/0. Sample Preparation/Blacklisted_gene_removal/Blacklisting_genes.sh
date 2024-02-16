# first step download of MIDAS repository using gitclone command 

git clone https://github.com/snayfach/MIDAS

##go inside the repository and download reference database from MIDAS

cd MIDAS
wget http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz 

##unzip the tar.gz file

tar -zxvf midas_db_v1.2.tar.gz

##install the dependancies

sudo python MIDAS/setup.py install

##update environment variables using this
export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
export PATH=$PATH:/path/to/MIDAS/scripts
export MIDAS_DB=/path/to/midas_database

##install numpy (required dependancy): https://numpy.org/install/

conda create -n my-env #create an environment for numpy 
conda activate my-env  #activate the environment

# The actual install command

conda install numpy
## script for gene gain/loss: 

#!/bin/bash
module load anaconda/3-2020.07
#SBATCH --partition=general
#SBATCH --job-name=my_conda_job
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=20000

export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
export PATH=$PATH:/path/to/MIDAS/scripts
export PYTHONPATH=$PYTHONPATH:/mnt/beegfs/apps/dmc/apps/anaconda_3-2020.07/lib/python3.8/site-packages

for fq1 in path_to_input_sample_files/*.paired.1.fq
   
do

echo "working with file $fq1"

base=$(basename $fq1.paired.1.fq)

echo "base name is $base"

fq1=/path_to_input_sample_files/${base}.paired.1.fq
fq2=/path_to_input_sample_files/${base}.paired.2.fq

run_midas.py species midas_${base} -1 $fq1 -2 $fq2 -d /path_to_MIDAS_reference_genome/midas_db_v1.2
done


######################### 
# Then we looked into the output, we took the species prevalence with than zero mean abundance from previous step output file.
#those species tags were {Pseudomonas_oleovorans_57108, Methylobacterium_sp_59341, Pseudomonas_fulva_57974, Brevundimonas_nasdae_60942, Methylobacterium_radiotolerans_54853 \
#	Alpha_proteobacterium_59626, Methylobacterium_sp_58573, Sphingomonas_phyllosphaerae_58541, Sphingomonas_sp_60678, Sphingomonas_taxi_60871, Sphingomonas_sp_61953
#Pantoea_agglomerans_54643, Pantoea_vagans_57743, Methylobacterium_oryzae_55240, Pseudomonas_fulva_58092, Methylobacterium_mesophilicum_62490
# Pantoea_sp_60701, Microbacterium_paraoxydans_56209, Pseudomonas_sp_59673, Xanthomonas_arboricola_57436, Sphingomonas_sp_58049, Xanthomonas_axonopodis_56719, Methylobacterium_extorquens_57587
#Xanthomonas_axonopodis_61257, Aureimonas_ureilytica_58716, Xanthomonas_axonopodis_57683, Pseudomonas_rhizosphaerae_61010, Pantoea_agglomerans_56951, Pseudomonas_fluorescens_61150, Leclercia_adecarboxylata_62497
# Enterobacter_cloacae_58148, Cronobacter_zurichensis_60329, Stenotrophomonas_maltophilia_62375, Streptomyces_sp_60263, Streptomyces_sp_58511, Rhodanobacter_sp_59306


##we took all the listed reference genome in midas database and concatenated all the listed genes files in the genome directories:

zcat ./Pan_genome/*/centroids.ffn.gz  >> all_species_pan_genomes_genes.fasta

##build database for blast
source /apps/profiles/modules_dmc.sh.dyn
module load blast+
makeblastdb -in all_species_pan_genomes_genes.fasta -out all_species_pan_genomes_genes_db -dbtype nucl


#blasted the database with the Xanthomonas perforans database from Midas
blastn -db all_species_pan_genomes_genes.fasta -query Xanthomonas_perforans_55843/centroids.ffn -outfmt 6 -out blast_output.txt 


#Then, we took all the gene with more than 97% identity from the blast result and extracted the list of those genes out of our databases, we called these genes as blacklisted genes
while read gene; do samtools faidx Xanthomonas_perforans_55843/centroids.ffn $gene >> Xp_within_host_gene_changes.fasta; done < Xp_within_host_gene_changes.txt

# "Xp_within_host_gene_changes.txt" file has list of genes
# "Xp_within_host_gene_changes.fasta" output file containing all the gene sequences with more than 97% identity
# "centroids.ffn" is input file or database of Xp

#next step was to remove all the blacklisted genes sequences from samples using bbduk

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bbmap

for fq1 in /path_to_input_sample_files/*.paired.1.fq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 .paired.1.fq)
    echo "base name is $base"

    fq1=/path_to_input_sample_files/${base}.paired.1.fq
    fq2=/path_to_input_sample_files/${base}.paired.2.fq


bbduk.sh -Xmx30g in1=$fq1 in2=$fq2 out1=${base}_blk_filter.1.fq out2=${base}_blk_filter.2.fq outm1=${base}_matched.1.fq outm2=${base}_matched.2.fq \
ref=pangenome_97_blk_cen.fasta k=31 hdist=1 stats=${base}.stats.txt

done













## 1. Installation: download of MIDAS repository using gitclone command 

git clone https://github.com/snayfach/MIDAS

##go inside the repository and download reference database from MIDAS

cd MIDAS
wget http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz 



##install the dependancies

sudo python MIDAS/setup.py install

##update environment variables using this
export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
export PATH=$PATH:/path/to/MIDAS/scripts
export MIDAS_DB=/path/to/midas_database

##install numpy (required dependancy): https://numpy.org/install/
# Best practice, use an environment rather than install in the base env

conda create -n my-env

conda activate my-env

# If you want to install from conda-forge

conda config --env --add channels conda-forge

# The actual install command

conda install numpy



#2.	Build custom database

#using Pangenome of AL65 and AL22 (made by superpang) for building database 

#According to github page, we need following files: indir: Path to directory of genomes. Each subdirectory should be named with a genome identifier. 
#Each subdirectory should contain the following files; 
#1. <genome_id>.fna: Genomic DNA sequence in FASTA format; 
# 2.<genome_id>.faa: Protein sequences in FASTA format; 
# 3. <genome_id>.ffn: Gene sequences in FASTA format; 
# 4. <genome_id>.genes: Tab delimited file with genomic coordinates of genes. 
#The file should be tab-delimited file with a header and the following fields:	gene_id (CHAR); scaffold_id (CHAR);	start (INT); end (INT); strand (+ or -); 	gene_type (CDS or RNA); mapfile: Path to mapping file that specifies which genomes belonging to the same species.
#5. The file should be tab-delimited file with a header and 3 fields:	genome_id (CHAR): corresponds to subdirectory within INDIR; species_id (CHAR): : species identifier for genome_id; rep_genome (0 or 1): indicator if genome_id should be used for SNP calling


#The first three files was already made before using prokka before running SuperPang.

#And the fourth one was made using gff  file from prokka using the following software.

#It was installed using conda. https://github.com/shenwei356/csvtk#installation 


conda create -n csvtk


#
# To activate this environment, use
#
#     $ conda activate csvtk
#
# To deactivate an active environment, use
#
#     $ conda deactivate


conda install -c bioconda csvtk


# we were able to convert a .gff file obtained from Prodigal to a .genes file with the correct columns by using the program csvtk and doing the steps below:
#To count number of columns in a .gff file:

csvtk dim --cols -t nrpan_AL22_AL65.gff

#Result 9
#Notes: The CSV parser requires all the lines have same number of fields/columns. Even lines with spaces will cause error. Use '-I/--ignore-illegal-row' to skip these lines if necessary. By default, csvtk handles CSV files, use flag -t for tab-delimited files.
#To select fields/columns:

csvtk cut -f 1,3-5,7,9 --ignore-illegal-row -t nrpan_AL22_AL65.gff > nrpan_AL22_AL65.genes


#Check the number of columns on the resulting .genes file:

csvtk dim --cols -t nrpan_AL22_AL65.genes

#Result: 6

#To rename fields/columns in genes file:

csvtk rename -f 1-6 -t -n scaffold_id,gene_type,start,end,strand,gene_id nrpan_AL22_AL65.genes > nrpan_AL22_AL65_renamed_columns.genes






### keep all the files in one folder required for building database.
#run the following command:


#!/bin/bash
source /apps/profiles/modules_dmc.sh.dyn
module load vsearch/2.14.1

module load anaconda/3-2020.11
#SBATCH --partition=general
#SBATCH --job-name=my_conda_job
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=20000
export PATH=$PATH:/path-to/MIDAS/hmmer-3.3.2
export PYTHONPATH=$PYTHONPATH:/mnt/beegfs/apps/dmc/apps/anaconda_3-2020.07/lib/python3.8/site-packages
build_midas_db.py /path-to/MIDAS/genomes  genomes.mapfile  midas_custom_db  --compress




#STEP 3. Next steps were to counts the SNVs and Gene gain and loss



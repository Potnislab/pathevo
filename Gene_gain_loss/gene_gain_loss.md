#        Gene gain and loss : MIDAS                       

## to identify the genes gain and loss, we used MIDAS (https://github.com/snayfach/MIDAS)
## first step download of MIDAS repository using gitclone command 

    git clone https://github.com/snayfach/MIDAS

##go inside the repository and download reference database from MIDAS

    cd MIDAS
    wget http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz 

## unzip the tar.gz file

    tar -zxvf midas_db_v1.2.tar.gz

## install the dependancies

    sudo python MIDAS/setup.py install

## update environment variables using this

    export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS
    export PATH=$PATH:/path/to/MIDAS/scripts
    export MIDAS_DB=/path/to/midas_database

## install numpy (required dependancy): https://numpy.org/install/
# Best practice, use an environment rather than install in the base env

    conda create -n my-env

    conda activate my-env

# If you want to install from conda-forge

    conda config --env --add channels conda-forge

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
    for fq1 in /path/to/paired_end/reads/*.paired.1.fq
    do
    echo "working with file $fq1"
    base=$(basename $fq1 .paired.1.fq)
    echo "basename is $base"
    fq1=/path/to/${base}.paired.1.fq
    fq2=/path/to/${base}.paired.2.fq
    run_midas.py snps midas_${base}_Snp --species_id Xanthomonas_perforans_55843 -1 $fq1 -2 $fq2 -d /path/to/database/midas_db_v1.2
    done
 
# midas_${base}_Snp: output files
# species_id: reference id from MIDAS database
# -d: path to reference database


# merging of the result output files

      module load anaconda/3-2020.07
      #SBATCH --partition=general
      #SBATCH --job-name=my_conda_job
      #SBATCH --cpus-per-task 16
      #SBATCH --mem-per-cpu=20000
  
      export PYTHONPATH=$PYTHONPATH:/scratch/aubaxk001/MIDAS/MIDAS/MIDAS
      export PATH=$PATH:/scratch/aubaxk001/MIDAS/MIDAS/MIDAS/scripts
      export PYTHONPATH=$PYTHONPATH:/mnt/beegfs/apps/dmc/apps/anaconda_3-2020.07/lib/python3.8/site-packages
      merge_midas.py genes ./merge_genes --species_id Xanthomonas_perforans_55843 -i ./midas_AALE_S2,./midas_BSCE_S4,./................................./midas_BALM_S5,./midas_CGAM_S13,./midas_HGAE_S16,./midas_SALM_S11 -d /path/to/MIDAS/midas_db_v1.2 -t list


output file = ./merge_genes
-i : with path details of input files, list names of output files from last step as an input for this step
species_id: reference id from MIDAS database
-d: path to reference database






################## GENE CHANGES: GAIN OR LOSS  ################# 

#The command for genes 

#!/bin/bash
module load anaconda/3-2020.07
#SBATCH --partition=general
#SBATCH --job-name=my_conda_job
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=20000
export PATH=$PATH:/path-to/Midas/midas_custom_db
export PYTHONPATH=$PYTHONPATH:/mnt/beegfs/apps/dmc/apps/anaconda_3-2020.07/lib/python3.8/site-packages

# here we have used normal samples (also includes the blacklisted genes reads)
for fq1 in /path-to-input-samples/*.paired.1.fq
do
echo "working with file $fq1"
base=$(basename $fq1 .paired.1.fq)
echo "base name is $base"
fq1=/path-to-input-samples/${base}.paired.1.fq
fq2=/path-to-input-samples/${base}.paired.2.fq
run_midas.py genes midas_${base}_Xp --species_id Xanthomonas_perforans -1 $fq1 -2 $fq2 -d /path-to/Midas/midas_custom_db
done


## merging of the result output files

module load anaconda/3-2020.07
#SBATCH --partition=general
#SBATCH --job-name=my_conda_job
#SBATCH --cpus-per-task 16
#SBATCH --mem-per-cpu=20000

export PYTHONPATH=$PYTHONPATH:/path_to/MIDAS
export PATH=$PATH:path_to/MIDAS/scripts
export PYTHONPATH=$PYTHONPATH:/mnt/beegfs/apps/dmc/apps/anaconda_3-2020.07/lib/python3.8/site-packages
 
merge_midas.py genes ./merge_genes --species_id Xanthomonas_perforans  \
-i 1EE_S42_Xp,1EM_S18_Xp,1XE_S43_Xp,1XM_S19_Xp,2EE_S44_Xp,2EM_S20_Xp,2XE_S45_Xp,2XM_S21_Xp,3EE_S46_Xp,3EM_S22_Xp, \ 3XE_S47_Xp,3XM_S23_Xp,4EE_S48_Xp,4EM_S24_Xp,4XE_S49_Xp,4XM_S25_Xp,5EE_S50_Xp,5EM_S26_Xp,5XE_S51_Xp,5XM_S27_Xp,6EE_S52_Xp,6EM_S28_Xp,6XE_S53_Xp,6XM_S29_Xp -d /path/to/MIDAS/midas_db_v1.2 -t list


## output file = ./merge_genes
## -i : with path details of input files, list names of output files from last step as an input for this step
## species_id: reference id from MIDAS database
## -d: path to reference database




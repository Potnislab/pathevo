#
########################################################################################
##############          Pangenome    with SuperPang           ###################
###########################################################################################

#Superpang. https://github.com/fpusan/SuperPang 
#Paper available in Biorxiv https://www.biorxiv.org/content/10.1101/2022.03.25.485477v1 
#Conda create environment in hpc home https://its.tntech.edu/display/MON/Installing+Miniconda+or+Anaconda+Environments+in+Your+HPC+Account
#Superpang install https://anaconda.org/fpusan/superpang 

#note: the version is important. So far the version 0.9.4beta1 works well. But the previous version has a bug.


#INSTALLATION

conda create -name Superpang


#To install this package with conda run: 
conda create -n SuperPang -c conda-forge -c bioconda -c fpusan superpang


#To activate this environment, use
#
#     $ conda activate SuperPang
#
# To deactivate an active environment, use
#
#     $ conda deactivate


#Step1: SuperPang required checkm output, so we need to run checkm first

#“If you are on a machine with <40 GB of memory, the --reduced_tree flag can be used which reduces the memory requirements to approximately 14 GB.”


#!/bin/bash
module load drep/3.2.2
checkm lineage_wf -t 14 --reduced_tree -x fna ./ ./checkm

# “fna” is the suffix of input files
# “./” is the input file path
# “./checkm” is the output file path and name

#Result will be shown in the last lines of the job file (checkmshSCRIPT.oxxxxxx) and also the output directory: 
 

#The checkm results used as SuperPang input should include the following headers and corresponding information, if not, it will show:“Option –-checkm: the input has to contain headers of “Bin Id', 'Completeness', and 'Contamination’”

#Bin_Id Completeness    Contamination

 

#Step 2: first remember to activate superpang before run it if installed in home.

Superpang.sh:
#!/bin/bash
SuperPang.py --fasta AL65.fna AL22.fna --checkm checkm_output.txt --output-dir output_withcheckm --force-overwrite --verbose-mOTUpan -t 10

#Note: if you have more genomes, the command might get killed. Therefore, it is good to add a file containing the names of the input genomes with paths in single line


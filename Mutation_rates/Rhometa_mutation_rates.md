#    RHOMETA- MUTATION RATES ESTIMATES    

## 1. download rhometa repository using git clone command

    git clone https://github.com/sid-krish/rhometa.git

##    dependancies installation

    conda create -name Nextflow

    conda activate Nextflow

## 1. installation using

    conda install -c bioconda nextflow
    #or 
    conda install -c "bioconda/label/cf201901" nextflow 


## 2. installed freebayes

    conda install -c "bioconda/label/cf201901" freebayes


## 2.  THETA ESTIMATION: MUTATION RATES 

        #!/bin/bash
        source /opt/asn/etc/asn-bash-profiles-special/modules.sh
        module load anaconda/3-2021.11
        module load samtools
        module load bcftools/1.13

        for file in /path/to/*.q20.bam;    ##making files is already mentioned in the Metapop.sh
        do tag=${file%.q20.bam};

        nextflow run theta_est.nf --bam "$file" --fa /path/to/reference.fasta --output_dir Theta_mutation_rates_output;
        done

# 3.  lookuptables

        #!/bin/bash
        module load anaconda/3-2019.07
        nextflow run lookup_table_gen.nf -theta 0.00006 --output_dir lookup_tables




# 4.  Recombination rates

        #!/bin/bash
        source /opt/asn/etc/asn-bash-profiles-special/modules.sh
        module load samtools
        module load freebayes/1.0.2
        module load bcftools/1.13

        for file in /path/to/*.q20.bam;
        do tag=${file%.q20.bam};

        nextflow run rho_est.nf --bam "$file" --fa /path/to/reference.fa --lookup_tables /path/to/lookup_tables --output_dir RhoEst_nrpan_q20bam_v2_theta_0.0001; 
        done

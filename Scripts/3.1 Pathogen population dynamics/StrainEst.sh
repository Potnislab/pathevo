#StrainEst
#This work was done using the High Performance Computing (HPC) system at the Alabama Supercomputer Center.

#StrainEst documentation used:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741664/


#Directory for reference genomes: {path}/ReferenceGenomes contains: XpAL22.fasta, XpAL33.fna, XpAL57.fna, XpAL65.fna, XpGEV993.fna, XpLH3.fna
        Assemblies: XpAL22.fasta (N/A), XpAL33.fna (GCA_007713965.1), XpAL57.fna (GCF_007713955.1), XpAL65.fna (GCA_007714115.1), XpGEV993.fna (GCA_001010005.1), XpLH3.fna (GCA_001908855.1)
#Directory for sequence representative: {path}/ReferenceGenomes/SR/Xp91-118.fna (GCF_000192045.2)


#!/bin/bash
module load anaconda/2-5.2.0
strainest mapgenomes ./ReferenceGenomes/*.f* ./ReferenceGenomes/SR/Xp91-118.fna MR.fna


#!/bin/bash
module load anaconda/2-5.2.0
strainest map2snp ./ReferenceGenomes/SR/Xp91-118.fna MR.fna snp.dgrp


#!/bin/bash
module load anaconda/2-5.2.0
strainest snpdist snp.dgrp snp_dist.txt hist.pdf


#!/bin/bash
module load anaconda/2-5.2.0
strainest snpclust snp.dgrp snp_dist.txt snp_clust.dgrp clusters.txt


#!/bin/bash
module load anaconda/2-5.2.0
strainest mapgenomes ./ReferenceGenomes/*.f* ./ReferenceGenomes/SR/Xp91-118.fna MA.fna


#!/bin/bash
module load anaconda/2-5.2.0
bowtie2-build -f MA.fna MA_Xp_SC



##NOW RUN ON WHATEVER SAMPLES WITH THE 6 MA_Xp_SC files. 



##AtDep_2021


#To get the .sam files:

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bowtie2/2.2.9
for fq1 in {path}/*.paired.1.fq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 .paired.1.fq)
    echo "base name is $base"

    fq1={path}/${base}.paired.1.fq
    fq2={path}/${base}.paired.2.fq

bowtie2 --very-fast --no-unal -x MA_Xp_SC -1  $fq1 -2 $fq2 -S ./${base}.sam

done



#To get the .bam files:


#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load samtools/1.11
for f in ./*.sam
 do
   echo $f
   base=$(echo ${f} | sed 's/.sam//')
   echo $base

   # Convert file from SAM to BAM format
   samtools view -b $f > ${base}.bam

   # Sort BAM file
   samtools sort ${base}.bam -o ${base}.sorted.bam

   # index the bam file
   samtools index ${base}.sorted.bam

 done



#Run StrainEst Loop

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/2-5.2.0
for file in ./*.sorted.bam
        do
        echo "working with file1 $file"

        base=$(basename $file .sorted.bam)
        echo "base name is $base"

strainest est snp_clust.dgrp $file ./${base}

done

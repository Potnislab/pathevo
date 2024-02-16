
########################################################################################
##############                  Calculation of Fst                    ###################
###########################################################################################

###      Step 1. Mpileup

# Be careful with the order of input files, which would be used in the cmh test script.

#!/bin/bash

module load samtools/1.11
samtools mpileup -B -Q 0 -f nrpan_AL65_AL22.fna \
1EM_S18_L004.nrpan.q20.bam \
2EM_S20_L004.nrpan.q20.bam \
3EM_S22_L004.nrpan.q20.bam \
4EM_S24_L004.nrpan.q20.bam \
5EM_S26_L004.nrpan.q20.bam \
6EM_S28_L004.nrpan.q20.bam \
1XM_S19_L004.nrpan.q20.bam \
2XM_S21_L004.nrpan.q20.bam \
3XM_S23_L004.nrpan.q20.bam \
4XM_S25_L004.nrpan.q20.bam \
5XM_S27_L004.nrpan.q20.bam \
6XM_S29_L004.nrpan.q20.bam \
1EE_S42_L004.nrpan.q20.bam \
2EE_S44_L004.nrpan.q20.bam \
3EE_S46_L004.nrpan.q20.bam \
4EE_S48_L004.nrpan.q20.bam \
5EE_S50_L004.nrpan.q20.bam \
6EE_S52_L004.nrpan.q20.bam \
1XE_S43_L004.nrpan.q20.bam \
2XE_S45_L004.nrpan.q20.bam \
3XE_S47_L004.nrpan.q20.bam \
4XE_S49_L004.nrpan.q20.bam \
5XE_S51_L004.nrpan.q20.bam \
6XE_S53_L004.nrpan.q20.bam \
> atDep21.nrpan.mpileup


# atDep21.nrpan.mpileup is output file name


###      Step 2. Generating a synchronized mpileup file and supsampling

#download the popoolation pacakge in terminal (because this pacakge contains specific script for subsampling)

wget https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip

#unzip the pacakge
unzip popoolation2_1201.zip


#!/bin/bash

java -ea -Xmx12g -jar /path_to_a_script_of_popoolation2/mpileup2sync.jar --input atDep21.nrpan.mpileup --output atDep21.nrpan.sync --fastq-type sanger --min-qual 20 --threads 12


##     Step 3: Subsampling to uniform coverage

git clone https://github.com/magicDGS/PoPoolation2_magicDGS

perl /Path_to_script_/PoPoolation2_magicDGS/subsample-synchronized.pl --method withoutreplace --max-coverage 2% --target-coverage 10 --input atDep21.nrpan.sync --output atDep21.nrpan.ss10.sync

###      Step 4. Identify indel regions

#!/bin/bash

perl /Path_to_script_/PoPoolation2_magicDGS/indel_filtering/identify-indel-regions.pl --input atDep21.nrpan.mpileup --output indel-regions.gtf --indel-window 5 -–min-count 2

##       Step 5.  Filtering indels

# download package for the script for filtering indels

wget https://sourceforge.net/projects/popoolation/files/popoolation_1.2.2.zip 

unzip popoolation_1.2.2.zip #unzip package


perl /Path_to_script_/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl –input atDep21.nrpan.ss10.sync --gtf indel-regions.gtf --output atDep21.nrpan.ss10.idf.sync

##   Step 6. FSt-sliding:
### Fst analysis

#!/bin/bash
perl /path_to_script_/PoPoolation2_magicDGS/fst-sliding.pl --input atDep21.nrpan.ss10.idf.sync --output atDep21.nrpan.ss10.idf.fst --min-count 2 --min-coverage 4 --max-coverage 120 \
--pool-size 400 --window-size 1000 --step-size 1000 --suppress-noninformative 


## Step 7. convert to igv

#!/bin/bash
perl /path_to_script_/PoPoolation2_magicDGS/export/pwc2igv.pl --input   atDep21.nrpan.ss10.idf.fst --output   Atdep.ss10.idf.igv



# CMH analysis script

## Data preparation

### Raw reads cleaning

    #!/bin/bash

    source /opt/asn/etc/asn-bash-profiles-special/modules.sh
    module load bbmap/37.36

    #input fastq files

    for fq1 in ./*_1.fastq.gz
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.fastq.gz)
    echo "base name is $base"

    fq1=./${base}_1.fastq.gz
    fq2=./${base}_2.fastq.gz


    bbduk.sh -Xmx1g in1=$fq1 in2=$fq2 out1=${base}.clean1.fq out2=${base}.clean2.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

    bbduk.sh -Xmx1g in1=${base}.clean1.fq in2=${base}.clean2.fq out1=${base}.clean_trim1.fq out2=${base}.clean_trim2.fq qtrim=rl trimq=10


    done

### Mapping and sorting

    #!/bin/bash
    source /opt/asn/etc/asn-bash-profiles-special/modules.sh
    module load bwa/0.7.12
    module load picard/1.79
    #input fastq files

    for fq1 in ./*.clean_trim1.fq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 .clean_trim1.fq)
    echo "base name is $base"

    fq1=./${base}.clean_trim1.fq
    fq2=./${base}.clean_trim2.fq

    bwa mem nrpan_AL65_AL22.fna $fq1 $fq2 > ${base}.nrpan.sam

    java -Xmx2g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/SortSam.jar I=${base}.nrpan.sam O=${base}.nrpan.sort.sam SORT_ORDER=coordinate

    done

### Filtering low quality reads and duplicates

    #!/bin/bash

    module load samtools/1.11
    module load picard/2.24.0

    for file in *.nrpan.sort.sam;
    do
    echo "working with file $file"

    base=$(basename $file .nrpan.sort.sam)
    echo "base name is $base"
    
    samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -Sb "$file" > ${base}.nrpan.q20.bam

    picard MarkDuplicates -I "${base}.nrpan.q20.bam" -O ${base}.nrpan.rmd.sort.bam -M ${base}.dupstat.txt --REMOVE_DUPLICATES true --REMOVE_SEQUENCING_DUPLICATES true

    done

### Mpileup

Be careful with the order of input files, which would be used in the cmh test script.

    #!/bin/bash

    module load samtools/1.11
    samtools mpileup -B -Q 0 -f nrpan_AL65_AL22.fna \
    1EM_S18_L004.nrpan.rmd.sort.bam \
    2EM_S20_L004.nrpan.rmd.sort.bam \
    3EM_S22_L004.nrpan.rmd.sort.bam \
    4EM_S24_L004.nrpan.rmd.sort.bam \
    5EM_S26_L004.nrpan.rmd.sort.bam \
    6EM_S28_L004.nrpan.rmd.sort.bam \
    1XM_S19_L004.nrpan.rmd.sort.bam \
    2XM_S21_L004.nrpan.rmd.sort.bam \
    3XM_S23_L004.nrpan.rmd.sort.bam \
    4XM_S25_L004.nrpan.rmd.sort.bam \
    5XM_S27_L004.nrpan.rmd.sort.bam \
    6XM_S29_L004.nrpan.rmd.sort.bam \
    1EE_S42_L004.nrpan.rmd.sort.bam \
    2EE_S44_L004.nrpan.rmd.sort.bam \
    3EE_S46_L004.nrpan.rmd.sort.bam \
    4EE_S48_L004.nrpan.rmd.sort.bam \
    5EE_S50_L004.nrpan.rmd.sort.bam \
    6EE_S52_L004.nrpan.rmd.sort.bam \
    1XE_S43_L004.nrpan.rmd.sort.bam \
    2XE_S45_L004.nrpan.rmd.sort.bam \
    3XE_S47_L004.nrpan.rmd.sort.bam \
    4XE_S49_L004.nrpan.rmd.sort.bam \
    5XE_S51_L004.nrpan.rmd.sort.bam \
    6XE_S53_L004.nrpan.rmd.sort.bam \
    > atDep21.nrpan.mpileup

### Generating a synchronized mpileup file and supsampling

    #!/bin/bash

    java -ea -Xmx12g -jar /scratch/aubnxp/popoolation2/mpileup2sync.jar --input atDep21.nrpan.mpileup --output atDep21.nrpan.sync --fastq-type sanger --min-qual 20 --threads 12

    perl /scratch/aubnxp/popoolation2/PoPoolation2_magicDGS/subsample-synchronized.pl --method withoutreplace --max-coverage 2% --target-coverage 10 --input atDep21.nrpan.sync --output atDep21.nrpan.ss10.sync

### Identify and remove indel regions

    #!/bin/bash

    perl /scratch/aubnxp/popoolation2/PoPoolation2_magicDGS/indel_filtering/identify-indel-regions.pl --input atDep21.nrpan.mpileup --output
    indel-regions.gtf --indel-window 5 -–min-count 2

    perl /scratch/aubnxp/popoolation2/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl –input atDep21.nrpan.ss10.sync --gtf indel-regions.gtf --output atDep21.nrpan.ss10.idf.sync

### Fst analysis

    #!/bin/bash
    perl /scratch/aubnxp/popoolation2/PoPoolation2_magicDGS/fst-sliding.pl --input atDep21.nrpan.ss10.idf.sync --output atDep21.nrpan.ss10.idf.fst --min-count 2 --min-coverage 4 --max-coverage 120 \
    --pool-size 400 --window-size 1000 --step-size 1000 --suppress-noninformative 

### CMH test

CMH test was done for different comparisons in regard to different hosts and environments. The numbers of comparisons correspond to the order when generating mpileup file.

> To detect allele frequency change in response to ozone environment, first compare it for resistant cultivar X10R in mid season,

    perl ./cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseX10R-ambelevM.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 7-8,10-9,11-12 --remove-temp

> then in end season.

    perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseX10R-ambelevE.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 19-20,22-21,23-24 --remove-temp 

> Next, we compare it for susceptible cultivar ECW in mid season,

    perl ./cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseECW-ambelevM.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 1-2,4-3,5-6 --remove-temp 

> then in end season.

    perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseECW-ambelevE.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 13-14,16-15,17-18 --remove-temp

> On the other hand, to detect allele frequency change in response to host defense, first compare it under ambient environment in mid season,

    perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseAmb-X10RecwM.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 1-7,4-10,5-11 --remove-temp

> then in end season.

    perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseAmb-X10RecwE.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 13-19,16-22,17-23 --remove-temp

> Next, we compare it under elevated ozone in mid season,

    perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseElev-X10RecwM.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 2-8,3-9,6-12 --remove-temp

> then in end season.

    perl ./cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseElev-X10RecwE.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 14-20,15-21,18-24 --remove-temp

### Convert cmh to GWAS

    for file in *.cmh; do tag=${file%.cmh}; perl ./export/cmh2gwas.pl --input $file --output $tag.gwas --min-pvalue 1.0e-20; done

### Plot gwas in R
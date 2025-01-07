
#!/bin/bash

# Load the samtools module
module load samtools/1.11

# Loop through all BAM files in the directory
for BAM_FILE in /home/aubaxk001/METAPOP_New_reads/bam/ino/*.pan.rmd.sort.bam; do
    # Extract the base name (remove file extension)
    BASE_NAME=/home/aubaxk001/METAPOP_New_reads/bam/ino/$(basename "$BAM_FILE" .pan.rmd.sort.bam)    # Construct the output file name
    OUTPUT_FILE="${BASE_NAME}.mpileup"
    # Run samtools mpileup
    samtools mpileup -B -Q 0 -f pan_al65_22.fasta  "$BAM_FILE" > "$OUTPUT_FILE"
    # Print progress
    echo "Processed: $BAM_FILE -> $OUTPUT_FILE"
done

 version popoolation_1.2.2  

#java -ea -Xmx12g -jar /scratch/aubaxk001/Pepper_metagenome/Blk-reads/Maping/Population2/popoolation2_1201/mpileup2sync.jar --input 1EE_S42.mpileup   --output 1EE_S42.sync --fastq-type sanger --min-qual 20 --threads 12

#perl ./PoPoolation2_magicDGS/subsample-synchronized.pl --method withoutreplace --max-coverage 2% --target-coverage 10 --input 1EE_S42.sync --output 1EE_S42.ss10.sync

perl Variance-sliding.pl --measure D --input    pileup/1EE_S42.mpileup --output 1EE_S42.file  --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000


#!/bin/bash

# Load the samtools module
module load samtools/1.11
samtools mpileup -B -Q 0 -f pan_al65_22.fasta  1EM_S18.pan.rmd.sort.bam  4EM_S24.pan.rmd.sort.bam  5EM_S26.pan.rmd.sort.bam  > replicate_wise/Sus_Amb_Mid.pileup
samtools mpileup -B -Q 0 -f pan_al65_22.fasta 1EE_S42.pan.rmd.sort.bam  4EE_S48.pan.rmd.sort.bam  5EE_S50.pan.rmd.sort.bam  > replicate_wise/Sus_Amb_End.pileup


samtools mpileup -B -Q 0 -f pan_al65_22.fasta 1XM_S19.pan.rmd.sort.bam  4XM_S25.pan.rmd.sort.bam  5XM_S27.pan.rmd.sort.bam  > replicate_wise/Resis_Amb_Mid.pileup
samtools mpileup -B -Q 0 -f pan_al65_22.fasta 1XE_S43.pan.rmd.sort.bam   4XE_S49.pan.rmd.sort.bam  5XE_S51.pan.rmd.sort.bam  > replicate_wise/Resis_Amb_End.pileup


samtools mpileup -B -Q 0 -f pan_al65_22.fasta 2EM_S20.pan.rmd.sort.bam  3EM_S22.pan.rmd.sort.bam 6EM_S28.pan.rmd.sort.bam > replicate_wise/Sus_Elev_Mid.pileup
samtools mpileup -B -Q 0 -f pan_al65_22.fasta 2EE_S44.pan.rmd.sort.bam  3EE_S46.pan.rmd.sort.bam   6EE_S52.pan.rmd.sort.bam > replicate_wise/Sus_Elev_End.pileup


samtools mpileup -B -Q 0 -f pan_al65_22.fasta 2XM_S21.pan.rmd.sort.bam  3XM_S23.pan.rmd.sort.bam  6XM_S29.pan.rmd.sort.bam > replicate_wise/Resis_Elev_Mid.pileup
samtools mpileup -B -Q 0 -f pan_al65_22.fasta 2XE_S45.pan.rmd.sort.bam  3XE_S47.pan.rmd.sort.bam 6XE_S53.pan.rmd.sort.bam > replicate_wise/Resis_Elev_End.pileup





#!/bin/bash

# Set parameters

# Directory containing mpileup files
OUTPUT_DIR="output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .mpileup files in the input directory
for FILE in pileup/*.mpileup; do
    # Extract the base filename (e.g., 1EE_S42 from 1EE_S42.mpileup)
    BASENAME=$(basename "$FILE" .mpileup)

    # Set output file name
    OUTPUT_FILE="$OUTPUT_DIR/$BASENAME.file"

    # Run the Variance-sliding.pl script
    perl Variance-sliding.pl \
        --measure D \
        --input "$FILE" \
        --output "$OUTPUT_FILE" \
        --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000

done

echo "All files processed!"



perl Variance-sliding.pl --measure D --input    test.pileup --output test.file --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000
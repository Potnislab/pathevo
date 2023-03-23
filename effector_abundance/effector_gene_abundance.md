# ShortBred: Relative abundaance of type three effectors      
using ShortBred (https://github.com/biobakery/shortbred). It needs usearch for running, which was installed on hpsc account using following commands

# 1. Make a bin folder at home
      cd ~
      mkdir bin

# 2. Install usearch

     wget https://drive5.com/downloads/usearch6.1.544_i86linux32.gz

unzip
     gunzip usearch6.1.544_i86linux32.gz

and make it executable
     chmod +x usearch6.1.544_i86linux32.gz
 
and make link

      ln -s usearch6.1.544_i86linux32 usearch

add export path in your script whenever you run the code

## 1. Shortbred Identify:
to make the reference markers file, for that you could do the following options
a. take all the sequences of genes of interest in one file and the genes for which marker sequences you want to make in other file
b. A file containing the genes specific to your species of interest and in the other file all the gene sequences present in other species of genus. 
By getting all other sequenes can be done by blasting the gene sequences from species of interest on NCBI website.

code for making database

      #!/bin/bash
      source /opt/asn/etc/asn-bash-profiles-special/modules.sh
      module purge
      module load shortbred/0.9.4
      module load blast+/2.6.0
      export PATH=$PATH:/path/to/usearch6.1.544_i86linux32
      shortbred_identify.py --goi /path/to/gene_of_interest_sequence.fna --ref /path/to/othersequences.fna \
      --markers /path/to/output.faa --usearch ~/bin/usearch6.1.544_i86linux32

# 2. ShortBred Quantify:
script for estimating relative abundance of effector genes
the input read files came from the script mentioned in https://github.com/Potnislab/AtDep_2021_metagenome/blob/main/Raw%20read%20processing/host-contamination-removal_script 
      
      #!/bin/bash
      source /opt/asn/etc/asn-bash-profiles-special/modules.sh
      module load shortbred/0.9.4  
      export PATH=$PATH:~/path/bin/usearch6.1.544_i86linux32
      for F in /path/to/*trimmed_kneaddata.paired.1.fastq;    
      do
      N=$(basename $F trimmed_kneaddata.paired.1.fastq) ;
      shortbred_quantify.py --markers /path/to/all_markers_sequences_geneInterest.faa --wgs $F  \
      --results $N.txt --tmp $N --threads 2 --usearch ~/path/bin/usearch6.1.544_i86linux32
      done

--markers [path to reference markers.faa file] 
--wgs [path for sample.fastq file] 
results [file name.txt] --tmp [output] 


# Visualization: Using HClust2

# installation
by git clone

      git clone https://github.com/SegataLab/hclust2

Or by conda
      #name the environment
      $ conda create –name hclust
      # activating the environment
      conda activate hclust
      #downloading the hclust2 using bioconda
      conda install -c bioconda hclust2

The output file got from shortbred quantify should be combined and the abundance values should be normalized.Usually the normalize the values for each sample by dividing the abundance of all genes with abundance of one conserved gene. 
In this case it was done using Avrbs2 abundance.

make your dataset by including first two rows with metadata and tsv format

       hclust2.py \
       -i input.txt \  ##input file
       -o Avrbs2.png \  ## name of output file
       --skip_rows 1 \   ##rows name which you want to skip from your dataset if there is no row then donot specify it 
       --f_dist_f correlation \   ##for making it distance matrix for features
       --s_dist_f braycurtis \  
       --dpi 300 \   ##specify the image quality
       --title 'Relative_abundance_of_Effectors' \  ##title name
       --flabel_size 8 \  ##font size of features
       --slabel_size 8 \  ##font size of samples
       --metadata_rows 2,3,4 \  ##rows in which your legend data is included
       --legend_file Avrbs2.legend.png \  ##legend’s file name
       --colorbar_font_size 12 \  ##colorbar font size
       --no_sclustering \  ##if you donot want dendrogram on samples
       --no_fclustering \  ## if you donot want the clustering on features
       --flinkage #if you want dendrogram with linkage

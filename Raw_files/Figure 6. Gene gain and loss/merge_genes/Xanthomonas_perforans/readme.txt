
Description of output files and file formats from 'merge_midas.py genes'

Output files
############
genes_depth.txt  
  average-read depth of each gene per sample
genes_copynum.txt
  copy-number of each gene per sample
  estimated by dividing the read-depth of a gene by the median read-depth of 15 universal single copy genes
genes_presabs.txt  
  the presence (1) or absence (0) of each gene per sample
  estimated by applying a threshold to gene copy-number values
genes_reads.txt
  number of reads mapped to each gene per sample
genes_summary.txt
  alignment summary statistics per sample

Output formats
############
genes_depth.txt, genes_copynum.txt, genes_presabs.txt, genes_reads.txt
  tab-delimited matrix files
  field names are sample ids
  row names are gene ids
genes_summary.txt
  sample_id: sample identifier
  pangenome_size: number of non-redundant genes in reference pan-genome
  covered_genes: number of genes with at least 1 mapped read
  fraction_covered: proportion of genes with at least 1 mapped read
  mean_coverage: average read-depth across genes with at least 1 mapped read
  marker_coverage: median read-depth across 15 universal single copy genes
  aligned_reads: number of reads that aligned to pangenome
  mapped_reads: number of aligned reads after applying filters for mapping quality, base quality, alignment fraction, and percent identity

Additional information for species can be found in the reference database:
 /scratch/aubaxk001/Midas/midas_blk_db/pan_genomes/Xanthomonas_perforans

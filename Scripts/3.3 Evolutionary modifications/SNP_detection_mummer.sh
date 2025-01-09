##### MUMMER####
#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load mummer/3.22

#first run SNP detection. It will tell which SNPs are present in spcific genome based on the reference genome (Pangenome of AL22 and AL65)

nucmer --prefix=ref_qry ref.fasta qry.fasta

show-snps -Clr ref_qry.delta > ref_qry.snps


#ref_qry pangenome of AL65 and AL22
#qry.fasta : AL65.fasta (in one run) and AL22.fasta (second time run); LH3.fasta (third run) separately


#Just change .snp to .txt and it will work. 

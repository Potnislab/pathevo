

##################        SNV CHANGES          ################# 


#!/bin/bash
module load anaconda/3-2020.07
export PYTHONPATH=$PYTHONPATH:/scratch/aubaxk001/Midas/snps/MIDAS
export PATH=$PATH:/scratch/aubaxk001/Midas/snps/MIDAS/scripts
export PATH=$PATH:/scratch/aubaxk001/Midas/midas_blk_db
export PATH=$PATH:/home/aubaxk001/.local/bin

for fq1 in /scratch/aubaxk001/Midas/atdep_reads/2/*.paired.1.fq
do
echo "working with file $fq1"
base=$(basename $fq1 .paired.1.fq)
echo "base name is $base"
fq1=/scratch/aubaxk001/Midas/atdep_reads/2/${base}.paired.1.fq
fq2=/scratch/aubaxk001/Midas/atdep_reads/2/${base}.paired.1.fq

run_midas.py snps ${base}_Snp_bk --species_id Xanthomonas_perforans \
-1 $fq1   -2 $fq2 -d /scratch/aubaxk001/Midas/midas_blk_db --remove_temp
done




## merging of the result output files based on the different treatments

#Susceptible-Ambient-end season

merge_midas.py snps ./Snps_SUS_A_End --species_id Xanthomonas_perforans -i 1EE_S42_Snp_bk,4EE_S48_Snp_bk,5EE_S50_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any

#Susceptible-Elevated-o3-end season

merge_midas.py snps ./Snps_SUS_E_End --species_id Xanthomonas_perforans -i 2EE_S44_Snp_bk,3EE_S46_Snp_bk,6EE_S52_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any

#Susceptible-Ambient-mid season
merge_midas.py snps ./Snps_SUS_A_Mid --species_id Xanthomonas_perforans -i 1EM_S18_Snp_bk,4EM_S24_Snp_bk,5EM_S26_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any


#Susceptible-Elevated-o3-Mid season
merge_midas.py snps ./Snps_SUS_E_Mid --species_id Xanthomonas_perforans -i 2EM_S20_Snp_bk,3EM_S22_Snp_bk,6EM_S28_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any


#Resistant-Ambient-end season

merge_midas.py snps ./Snps_Resis_A_End --species_id Xanthomonas_perforans -i 1XE_S43_Snp_bk,4XE_S49_Snp_bk,5XE_S51_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any


#Resistant-Elevated-o3-end season

merge_midas.py snps ./Snps_Resis_E_End --species_id Xanthomonas_perforans -i 2XE_S45_Snp_bk,3XE_S47_Snp_bk,6XE_S53_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any

#Resistant-ambient-mid season
merge_midas.py snps ./Snps_Resis_A_Mid --species_id Xanthomonas_perforans -i 1XM_S19_Snp_bk,4XM_S25_Snp_bk,5XM_S27_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any

#Resistant-Elevated-o3-mid season
merge_midas.py snps ./Snps_Resis_E_Mid --species_id Xanthomonas_perforans -i 2XM_S21_Snp_bk,3XM_S23_Snp_bk,6XM_S29_Snp_bk -t list -d midas_blk_db --all_sites --all_samples  --snp_type any





## output file = ./Snps_X10R_E_Mid
## -i : with path details of input files, list names of output files from last step as an input for this step
## species_id: reference id from MIDAS database
## -d: path to reference database




#################################################################################

## Downstream analysis for further details
#convert all the text files into csv

sed 's/\t/,/g'  input.txt > output.csv
# 1. combined the snp depth, snp frequency for each site per replicate within treatment separately using this commnad:


awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' Snp_depth.csv Snp_freq.csv > Snp_freq_depth.csv

#now separate them for each replicate

cut  -f 1,2,5 Snp_freq_depth.csv > 1EE.csv # (for example)

# do it for all replicates

# 2. keep 6 columns from snps_info files

cut  -f 1,2,3,4,5,6 Snp_info.csv > SNps_info.csv


#3. combine also mid and end season snp files using this command: 

awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' snp_info_mid.csv snp_info_end.csv > Snp_info.csv
cut -d, --complement -f 7,8,9,10 Snp_info.csv > Snp_info1.csv

# do for all the replicates

#4. now we have columns like this :

# site_id, ref_id, ref_pos, ref_allele, major_allele_mid, minor_allele_mid, major_allele_end, minor_allele_end

# to find if there are same major from mid to end season: creating new column, by checking if each row of column 5 is equal to column 7, then type true

awk -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "Same_major"; print; next } {if ($5==$7) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, Snp_info1.csv > Snp_info_SM.csv

# to find if there are same minor allele from mid to end season: creating new column
awk  -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "Same_minor"; print; next } {if ($6==$8) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, Snp_info_SM.csv > SMin.csv

# to find allele is shifting from minor to major from mid to end season: creating new column

awk  -F, -vOFS=, 'NR==1 { $0 = $0 OFS  "Min_major"; print; next } {if ($6==$7) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, SMin.csv> MMaj_.csv

# to find allele is shifting from major to minor from mid to end season: creating new column

awk -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "maj_minor"; print; next } {if ($5==$8) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, MMaj_.csv > Final.csv

# to find de novo mutations during mid season; by checking if reference allele is same as mid season major allele
awk  -F, -vOFS=, 'NR==1 { $0 = $0 OFS  "Denovo_Mid"; print; next } {if ($4==$5) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, Final_.csv  >  denovo_Mid.csv


# to find de novo mutations during end season; by checking if reference allele is same as end season major allele
awk  -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "Denovo_End"; print; next } {if ($4==$7) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, denovo_Mid.csv > Denovo_Final.csv


#now we have these column in our final file:
#site_id,ref_id,ref_pos,ref_allele,major_allele_mid,minor_allele_mid,major_allele_end,minor_allele_end,
#9 Same_major,
#10 Same_minor,
#11 Min_major,
#12 maj_minor

# to find for if there is shift in allele for specific sites: if the minor allele is changing or not similar:

awk -F , '$10=="False" { print }' Final.csv > minor.csv

#from above output,  minor alleles which were changed, which ones from those are actual cause of allele shift: meaning which one minor become major from mid minor allele to end season as amjor allele:

awk -F , '$11=="True" { print }' minor.csv > Changed_minor_to_major.csv



# for denovo mutations: we have column in Denovo_Final.csv:
#site_id,ref_id,ref_pos,ref_allele,major_allele_mid,minor_allele_mid,major_allele_end,minor_allele_end,
#9 Same_major,
#10 Same_minor,
#11 Min_major,
#12 maj_minor,
#13 Denovo_Mid,
#14 Denovo_End


# for mid season de novo mutations:

awk -F , '$13=="False" { print }' Denovo_Final.csv > Denovo_Final_mid.csv 

#deleterious might be those when new mutation (major mid) did not pass to end of season
awk -F , '$9=="False" { print }' Denovo_Final_mid.csv > Denovo_Final_deleterious.csv 


# for end season de novo mutations:

awk -F , '$14=="False" { print }' Denovo_Final.csv  > Denovo_Final_end.csv 

#now combine the files of snp_depth and frequency files with these final files for all replicates
# and remove the sites which have read-depth is less than 10 reads
#using R
#then we can merge the files from Mummer to find which genome is belonging to specific alleles
# and find if all of them are truely de novo and also to locate the allele shift.








########################################################################################
##############                  CMH Test                              ###################
###########################################################################################


### CMH test: Where (according to the table and picture shown above): 7 represent 1XM which is X10R ambient mid season replicate 1; 8 represent 2XM which is X10R elevated mid season replicate \
# 10 represent 4XM which is X10R ambient mid season replicate 2; 9 represent 3XM which is X10R elevated mid season replicate 2
# 11 represent 5XM which is X10R ambient mid season replicate 3; 12 represent 6XM which is X10R elevated mid season replicate 3


# Note: The number in the comparison depends on the order of input files when mpileup (see step1 from the Fst test). The table above shows the number for each input files in the 1st column.


#CMH test was done for different comparisons in regard to different hosts and environments. The numbers of comparisons correspond to the order when generating mpileup file.

# To detect allele frequency change in response to ozone environment, first compare it for resistant cultivar X10R in mid season,

#input is from Fst test last step = atDep21.nrpan.ss10.idf.sync

perl ./cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseResis-ambelevM.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 7-8,10-9,11-12 --remove-temp

# then in end season.

perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseResis-ambelevE.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 19-20,22-21,23-24 --remove-temp 

# Next, we compare it for susceptible cultivar in mid season,

perl ./cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseSuscep-ambelevM.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 1-2,4-3,5-6 --remove-temp 

# then in end season.

perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseSuscep-ambelevE.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 13-14,16-15,17-18 --remove-temp

# On the other hand, to detect allele frequency change in response to host defense, first compare it under ambient environment in mid season,

perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseAmb-Resis_Suscep_Mid.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 1-7,4-10,5-11 --remove-temp

# then in end season.

perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseAmb-Resis_Suscep_End.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 13-19,16-22,17-23 --remove-temp

# Next, we compare it under elevated ozone in mid season,

perl cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseElev-Resis_Suscep_Mid.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 2-8,3-9,6-12 --remove-temp

# then in end season.

perl ./cmh-test.pl --input atDep21.nrpan.ss10.idf.sync --output baseElev-Resis_Suscep_End.cmh --min-count 2 --min-coverage 4 --max-coverage 120 --population 14-20,15-21,18-24 --remove-temp





### Convert cmh to GWAS

for file in *.cmh; do tag=${file%.cmh}; perl ./export/cmh2gwas.pl --input $file --output $tag.gwas --min-pvalue 1.0e-20; done

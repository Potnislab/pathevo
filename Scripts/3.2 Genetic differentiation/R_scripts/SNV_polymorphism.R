
---
title: "SNP Polymorphism"

---

#run in terminal on Alabamae Supercomputer
#first I combined  the SNP_frequeny and SNP depth file for each sample

# first make directory and also for each file add a new column for each file
mkdir modified_files && for i in *.csv; do   awk 'NR==1{ sub(/\.csv$/, "", FILENAME) } { $0 = $0 OFS FILENAME }1' OFS=,  "$i"  > "./modified_files/$i"; done

# make directory where all the site will be added to include only the sites which have SNP depth of >= 10

mkdir read_depth && for i in *.csv; do   awk -F, 'NR==1 || $3>=10' $i  > read_depth/$i; done

# make directory where all the site will be added to include only the sites which have greater than 0.2 and less than 0.8 SNP freuqncy 

mkdir Polymorphism_files && for i in *.csv; do   awk -F, 'NR==1 || $2>0.2 &&  $2<0.8' $i  > Polymorphism_files/$i; done

#concatenation of all files in minor and polymorphism files without header

awk 'FNR>2 || NR==2' ./*.csv > ALL-Merged.csv


#to remove the rows based on the first column of other file

awk  -F, -vOFS=, 'NR == FNR {a[$1]; next} !($1 in a)' /scratch/aubaxk001/FInal_midas/blk_sites  ALL-Merged.csv >  ALL-Merged_blk.csv


#to count the rows based on factors. Of columns

awk -F, '$3{array[$2]++;array2[$2" "$3]++;array3[$3]++} END{for(u in array){for(y in array3){if(array2[u" "y]){print u,array[u],y,array2[u" "y]}}}}'   ALL-Merged_blk1.csv > data.csv



################################
#NOW RUN below command on RSTUDIO to make Figure

```{r}
file <- readxl::read_xlsx("/path_to_file_with_normalised_counts/SNV_poly_0.2_0.8.xlsx", sheet=2)

#Note: we normalised the raw counts by dividing the absolute abundance of each sample with respective values


file$Cul_Con <- paste(file$Cultivar, file$Condition)
order <- factor(file$Sample, levels= c("1EM" , "4EM","5EM", "1EE" ,"4EE","5EE", "1XM", "4XM"  , "5XM", "1XE" ,"4XE", "5XE" ,
                                       "2EM", "3EM", "6EM", "2EE" , "3EE" ,"6EE",  "2XM",  "3XM" ,  "6XM","2XE" ,"3XE", "6XE"))

library(RColorBrewer)
library(ggplot2)


colr <- c("tomato", "steelblue", "mediumvioletred", "lightgreen", "slategray")

ggplot(file, aes(fill=site_type, y=Norm_count, x = order)) + 
    geom_bar(position="stack", stat="identity",  width = 0.60 ) +
    facet_wrap( ~ factor(Cul_Con, levels =c( "Susceptible Ambient" , "Susceptible Elevated", "Resistant Ambient"  ,  "Resistant Elevated"  )), scale = "free_x", nrow = 1) +
 # facet_wrap( ~ factor(Cultivar, levels = c("Susceptible", "Resistant")) ) + #, switch = "x") +
 scale_fill_manual(values = colr) +
  theme_bw() +  labs( x = "", y = "SNV counts per Xp density \n(ng of DNA per mg of leaf sample)") + 
  theme(axis.text = element_text(face="bold", size=8), 
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        title =element_text(face="bold", size=8),
        legend.title = element_text(face="bold", size=8),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=8)) +
  #ggtitle("SNV counts with intermediate frequency (range of 0.2: 0.8) \n with atleast 10 reads and without blk genes") +
  theme(legend.position = "bottom"  ) 

#ggsave("/Users/amanpreetkaur/Downloads/SNV_poly_norm.pdf",width = 7, height = 4, units="in", dpi=700)



```
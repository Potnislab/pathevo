---
title: "CMH_plot"
---

---
title: "cmh_simple"
author: "Ranlin"
date: "2023-09-19"
output: html_document
---

```{r setup, include=FALSE}
# packages required

library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(writexl)
library(cowplot)
library(readr)
library(dplyr)
library(purrr)
getwd()

```

```{r}
# Set the directory path to your folder containing gwas extension files from the CMH test last step

file_list <- list.files(path = folder_path, pattern = "\\.gwas$", full.names = TRUE)

# Read all files and combine them into a single data frame
combined_data <- map_dfr(file_list, ~ read_tsv(.x) %>% mutate(file_name_column = basename(.x)))

#remove NA
combined_data1 <- na.omit(combined_data)

is.na(combined_data1$P)

```

```{r cars}
cmh <- combined_data1
cmh$log <- with (cmh, -log10(P))

#extract chr number
# Load the required library
library(stringr)

# Extract numbers between underscores
cmh$Contig <- str_extract(cmh$CHR, "(?<=_)[0-9]+(?=_)")

# Print the extracted numbers
unique(cmh$Contig)

##contig length info
contig_len <- unique(cmh$CHR)
contig_len <- as.data.frame(contig_len)
contig_len$Contig <- str_extract(contig_len$contig_len, "(?<=_)[0-9]+(?=_)")
contig_len$len <- str_extract(contig_len$contig_len, "(?<=\\=)[0-9]+")

#generate continuious contig positions
data_cum <- contig_len %>%
  mutate(bp_add = lag(cumsum(len), default = 0)) %>%
  select(Contig, bp_add)

cmh <- cmh %>%
  inner_join(data_cum, by="Contig") %>%
  mutate(bp_cum=BP+bp_add)
head(cmh)
```

```{r data cleanup}
#change names https://stackoverflow.com/questions/29271549/replace-all-occurrences-of-a-string-in-a-data-frame 
unique(cmh$file_name_column)
cmh <- cmh %>%
  mutate( Population_comparisons= case_when(
    file_name_column == "baseAmb-X10RecwE.gwas" ~ "Amb_HOST_End",
file_name_column == "baseAmb-X10RecwM.gwas" ~ "Amb_HOST_Mid",
file_name_column == "baseElev-X10RecwE.gwas" ~ "Elev_HOST_End",
file_name_column == "baseElev-X10RecwM.gwas" ~ "Elev_HOST_Mid",
 file_name_column == "baseECW-ambelevE.gwas" ~ "S_ENV_End",
file_name_column == "baseECW-ambelevM.gwas" ~ "S_ENV_Mid",
file_name_column == "baseX10R-ambelevE.gwas" ~ "R_ENV_End",
file_name_column == "baseX10R-ambelevM.gwas" ~ "R_ENV_Mid",
TRUE ~ NA_character_  # Keep other values unchanged
  ))

##subset
cmh_h <- subset(cmh, Population_comparisons %in% c("Amb_HOST_End","Amb_HOST_Mid","Elev_HOST_End","Elev_HOST_Mid")) 
unique(cmh_h$Population_comparisons)

cmh_e <- subset(cmh, Population_comparisons %in% c("S_ENV_End","S_ENV_Mid","R_ENV_End","R_ENV_Mid"))
unique(cmh_e$Population_comparisons)

str(cmh_h)
str(cmh_e)

##create subsets with only significant snps
#generate threshold for each subset
cutoff <- as.data.frame(table(cmh$Population_comparisons))
total
cutoff$threshold <- 0.05/cutoff$Freq
cutoff

write_xlsx(cutoff,"cutoff.xlsx")

```



```{r plot}
cmh_h$Population_comparisons<-as.factor(cmh_h$Population_comparisons)

#host_combined
values1<- c("#1F78B4","#A6CEE3","#33A02C", "#B2DF8A")
levels(cmh_h$Population_comparisons)

h <- ggplot(cmh_h, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.25) +
  scale_color_manual(values = values1, 
                     labels =c("Ambient End: Susceptible vs Resistant", 
                               "Ambient Mid: Susceptible vs Resistant",
                               expression(paste("Elevated O"[3], " End: Susceptible vs Resistant")), 
                               expression(paste("Elevated O"[3], " Mid: Susceptible vs Resistant"))))+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_End",3]), color="#1F78B4", linetype="solid")+  #amb_host_end
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_Mid",3]), color="#A6CEE3", linetype="solid")+  #amb_host_mid
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_End",3]), color="#33A02C", linetype="solid")+  #elev_host_end
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_Mid",3]), color="#B2DF8A", linetype="solid")+  #elev_host_mid
  xlab("Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
    theme(legend.text = element_text(size = 3),
        legend.title = element_text(size = 4),
        legend.key.size = unit(0.1, "cm"),
        legend.text.align = 0,
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))

h_nk <-ggplot(cmh_h, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.5) +
  scale_color_manual(values = values1)+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_End",3]), color="#1F78B4", linetype="solid")+  #amb_host_end
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_Mid",3]), color="#A6CEE3", linetype="solid")+  #amb_host_mid
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_End",3]), color="#33A02C", linetype="solid")+  #elev_host_end
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_Mid",3]), color="#B2DF8A", linetype="solid")+  #elev_host_mid
  xlab("Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
    theme(legend.position = "none",
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 5))
ggsave("cmh_h.pdf", width=6, height=3.5, units=c("in"), dpi=700)

```


```{r}
##Amb_host
cmh_h_amb <- subset(cmh_h, Population_comparisons %in% c("Amb_HOST_Mid" ,"Amb_HOST_End"))
unique(cmh_h_amb$Population_comparisons)

ha <- ggplot(cmh_h_amb, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.25) +
  scale_color_manual(values = values1[1:2], 
                     labels =c("Ambient End: Susceptible vs Resistant", 
                               "Ambient Mid: Susceptible vs Resistant"))+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_Mid",3]), color="#A6CEE3", linetype="solid")+#amb_host_mid
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_End",3]), color="#1F78B4", linetype="solid")+ #amb_host_end
  scale_x_continuous(name= "Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
    theme(legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        legend.text.align = 0,
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))

ha_nk <- ggplot(cmh_h_amb, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.5) +
  scale_color_manual(values = values1[1:2])+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_Mid",3]), color="#A6CEE3", linetype="solid")+#amb_host_mid
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Amb_HOST_End",3]), color="#1F78B4", linetype="solid")+ #amb_host_end
  scale_x_continuous(name= "Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
    theme(legend.position = "none",
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 5))

#ggsave("cmh_h_a.pdf", width=6, height=3.5, units=c("in"), dpi=700)

```


```{r}
##Elev_host
cmh_h_elev <- subset(cmh_h, Population_comparisons %in% c("Elev_HOST_Mid" ,"Elev_HOST_End"))

he <- ggplot(cmh_h_elev, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.25) +
  scale_color_manual(values = values1[3:4], 
                     labels =c(
                               expression(paste("Elevated O"[3], " End: Susceptible vs Resistant")), 
                               expression(paste("Elevated O"[3], " Mid: Susceptible vs Resistant"))))+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_Mid",3]), color="#B2DF8A", linetype="solid")+ #elev_host_mid
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_End",3]), color="#33A02C", linetype="solid")+ #elev_host_end
  scale_x_continuous(name= "Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
  theme(legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.2, "cm"),
        legend.text.align = 0,
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))

he_nk <- ggplot(cmh_h_elev, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.5) +
  scale_color_manual(values = values1[3:4])+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_Mid",3]), color="#B2DF8A", linetype="solid")+ #elev_host_mid
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="Elev_HOST_End",3]), color="#33A02C", linetype="solid")+ #elev_host_end
  scale_x_continuous(name= "Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
  theme(legend.position = "none",
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 5))

ggsave("cmh_h_e.pdf", width=6, height=3.5, units=c("in"), dpi=700)

```


```{r}
#env_combined
values2 <- c("#E31A1C","#FB9A99","#FF7F00", "#FDBF6F" )
unique(cmh_e$Population_comparisons)

e <-ggplot(cmh_e, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.25) +
  scale_color_manual(values =values2, 
                     labels =c(expression(paste("Resistant End: Ambient vs Elevated O"[3])), 
                               expression(paste("Resistant Mid: Ambient vs Elevated O"[3])),
                               expression(paste("Susceptible End: Ambient vs Elevated O"[3])),
                               expression(paste("Susceptible Mid: Ambient vs Elevated O"[3]))))+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="S_ENV_Mid",3]), color="#FDBF6F", linetype="solid")+ #"S_ENV_Mid"
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="S_ENV_End",3]), color="#FF7F00", linetype="solid")+ #"S_ENV_End"
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="R_ENV_Mid",3]), color="#FB9A99", linetype="solid")+ #"R_ENV_Mid"
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="R_ENV_End",3]), color="#E31A1C", linetype="solid")+ #"R_ENV_End"
  scale_x_continuous(name= "Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
  theme(legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.key.size = unit(0.05, "cm"),
        legend.text.align = 0,
        axis.title = element_text(size = 5),
        axis.text = element_text(size = 4))

e_nk <-ggplot(cmh_e, aes(x=bp_cum, y= log, color=Population_comparisons)) + 
  geom_point(alpha=1, size=0.5) +
  scale_color_manual(values =values2 )+
  geom_hline(yintercept = -log10(5e-8), color="grey", linetype="dashed", size=0.25)+
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="S_ENV_Mid",3]), color="#FDBF6F", linetype="solid")+ #"S_ENV_Mid"
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="S_ENV_End",3]), color="#FF7F00", linetype="solid")+ #"S_ENV_End"
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="R_ENV_Mid",3]), color="#FB9A99", linetype="solid")+ #"R_ENV_Mid"
  geom_hline(yintercept = -log10(cutoff[cutoff$Var1=="R_ENV_End",3]), color="#E31A1C", linetype="solid")+ #"R_ENV_End"
  scale_x_continuous(name= "Contig Position") +
  ylab("-log10(p-value)")+
  theme_classic()+
  labs(color="Population Comparisons")+
  theme(legend.position = "none",
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 5))

ggsave("cmh_env.pdf", width=6, height=3.5, units=c("in"), dpi=700)

```



```{r combine figs}

a <-plot_grid(h,e, labels=c("i", "ii"), nrow=2, align = "h", label_size = 6)
h_ae <-plot_grid(ha,he, labels=c("iii", "iv"), nrow=2, align = "h", label_size = 6)
all <- plot_grid(a, h_ae,labels = "",ncol=2)

h_e_nk <-plot_grid(h_nk,e_nk, labels="AUTO", nrow=2, align = "h")
h_ae_nk <-plot_grid(ha_nk,he_nk, labels=c("C", "D"), nrow=2, align = "h")
all_nk <- plot_grid(h_e_nk, h_ae_nk,labels = "",ncol=2)

```

```{r}
ggsave("cmh_v2.pdf", width=6, height=3.5, units=c("in"), dpi=700)
```

```{r}
```


### Annotation of SNPs

```{r setup, include=FALSE}
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) #need to first run this line to increase memory then load libraries.
library(xlsx)
library(writexl)
library(dplyr)
library(tidyverse)
```

```{r}

filtered_cmh <- cmh %>%
  filter((Population_comparisons == "Amb_HOST_Mid" & P < cutoff[cutoff$Var1=="Amb_HOST_Mid",3]) |  # Subset A
         (Population_comparisons == "Amb_HOST_End" & P < cutoff[cutoff$Var1 == "Amb_HOST_End",3]) |  
         (Population_comparisons == "Amb_HOST_Mid" & P < cutoff[cutoff$Var1== "Amb_HOST_Mid",3]) |
        (Population_comparisons == "Elev_HOST_End" & P < cutoff[cutoff$Var1== "Elev_HOST_End",3]) |
          (Population_comparisons == "Elev_HOST_Mid" & P < cutoff[cutoff$Var1== "Elev_HOST_Mid",3]) |
          (Population_comparisons == "S_ENV_Mid" & P < cutoff[cutoff$Var1== "S_ENV_Mid",3]) |
          (Population_comparisons == "S_ENV_End" & P < cutoff[cutoff$Var1== "S_ENV_End",3]) |
          (Population_comparisons == "R_ENV_Mid" & P < cutoff[cutoff$Var1== "R_ENV_Mid",3]) |
          (Population_comparisons == "R_ENV_End" & P < cutoff[cutoff$Var1== "R_ENV_End",3]))  

total.no<-as.data.frame(table(filtered_cmh$Population_comparisons))
total.no$Freq<- as.numeric(total.no$Freq)
total.no$Var1<- as.character(total.no$Var1)
total.no[5:8,"Var1"]<-c("R_ENV_End","R_ENV_Mid","S_ENV_End","S_ENV_Mid")
total.no[5:8,"Freq"]<-c(0,3,0,0)
str(total.no)
total.no$Gene_No <- c(19,69,15,59,0,2,0,0)
total.no$Freq2 <- c(total.no$Freq-total.no$Gene_No)
no <- gather(total.no,key = "Type", value = "No",-c("Var1"))
no1 <- no[-c(which(no$Type=="Freq")),]
no2 <- no[-c(which(no$Type=="Freq2")),]
str(no1)
color1<-c("#1F78B4","#A6CEE3","#33A02C", "#B2DF8A","#E31A1C","#FB9A99","#FF7F00", "#FDBF6F")
color2<- c("black", "grey")
color3<-c("#1F78B4","#A6CEE3","#33A02C", "#B2DF8A","#E31A1C","#FB9A99","#FF7F00", "#FDBF6F","black", "grey")
library("ggpattern")   
library("gridpattern") 
library(RColorBrewer)

color4 <-brewer.pal(12,"Paired")
no2$Type <- gsub("Gene_No","Gene",no2$Type)
no2$Type <- gsub("Freq","SNP",no2$Type)

ggplot(no2,aes(x=Var1,y=No,fill=Type))+
  geom_bar(stat = "identity",position = "dodge")+
  labs(y="Number of significant SNPs/Genes",
       x="Comparisons")+
  geom_text(aes(label=No),stat = "identity",position = position_dodge(width = 0.9), vjust=-0.5, size=3)+
  theme_classic()+
  scale_fill_manual(values = color2)+
  theme(legend.position = "right",
        axis.text.x = element_text(size = 5))
ggsave("no.snps_genes.png", width=6, height=3.5, units=c("in"), dpi=700)

retain <- read.xlsx("Counts.xlsx",sheetIndex = 1, header=T)
retain.long <- gather(retain, key = "Type", value = Freq, -c("Comparison"))
ggplot(retain.long,aes(x=Comparison,y=Freq,fill=Type))+
  geom_bar(stat = "identity",position = "dodge")+
  labs(y="Number of significant SNPs/Genes",
       x="Comparisons")+
  geom_text(aes(label=Freq),stat = "identity",position = position_dodge(width = 0.9), vjust=-0.5, size=3)+
  theme_classic()+
  scale_fill_manual(values = color4[9:12])+
  theme(legend.position = "right",
        axis.text.x = element_text(size = 8))
ggsave("no.snps_genes_retained.png", width=6, height=3.5, units=c("in"), dpi=700)

cmh_s_mid<-cmh[cmh$Population_comparisons=="S_ENV_Mid",]
cmh_s_mid[cmh_s_mid$P<cutoff[8,3],]

gwas<- filtered_cmh[, c(1,2,4,7,10)]
str(gwas)
gwas$Contig <- as.numeric(gwas$Contig)

```

```{r}
gff <- read.xlsx("./Annotations.xlsx", sheetIndex = 3, header = T)
gff$CHR <- str_extract_all(gff$Contig, "(?<=_)[0-9]+(?=_)")
gff$CHR<-as.numeric(gff$CHR)
str(gff)


newdata<-matrix(0,nrow = 0,ncol = 8)
colnames(newdata)<-c("CHR","BP","P","Contig","Population_comparisons","Start","End","Annotation")
newdata<-as.data.frame(newdata)

for (gw in 1:nrow(gwas)) {
  temp <- subset(gff, CHR==gwas$Contig[gw] & Start<=gwas$BP[gw] & End >= gwas$BP[gw])
  if (nrow(temp) == 0) {
    rows <- cbind(gwas[gw,], NA, NA, NA)  # Adjust the number of NA columns
    colnames(rows) <- colnames(newdata)  # Match column names
  } else if(nrow(temp) == 1){
  rows <- cbind(gwas[gw,], temp[,c(2,3,4)])} else{
    temp2<-rbind(gwas[gw,],gwas[gw,])
    rows <- cbind(temp2, temp[,c(2,3,4)])
  }
  newdata <- rbind(newdata, rows)
}

newdata<-matrix(0,nrow = 0,ncol = 8)
colnames(newdata)<-c("CHR","BP","P","Contig","Population_comparisons","Start","End","Annotation")
newdata<-as.data.frame(newdata)
for (gw in 1:nrow(gwas)) {
  temp <- subset(gff, CHR==gwas$Contig[gw] & Start<=gwas$BP[gw] & End >= gwas$BP[gw])
  if (nrow(temp) == 0) {
    rows <- cbind(gwas[gw,], NA, NA, NA)  # Adjust the number of NA columns
    colnames(rows) <- colnames(newdata)  # Match column names
  } else if(nrow(temp) == 1){
  rows <- cbind(gwas[gw,], temp[,c(2,3,4)])} else{
    temp2<-rbind(gwas[gw,],gwas[gw,])
    rows <- cbind(temp2, temp[,c(2,3,4)])
  }
  newdata <- rbind(newdata, rows)
}

ahm_s <- subset(newdata, Population_comparisons=="Amb_HOST_Mid")
ahe_s <- subset(newdata, Population_comparisons=="Amb_HOST_End")
ehm_s <- subset(newdata, Population_comparisons=="Elev_HOST_Mid")
ehe_s <- subset(newdata, Population_comparisons=="Elev_HOST_End")
sem_s <- subset(newdata, Population_comparisons=="S_ENV_Mid")
see_s <- subset(newdata, Population_comparisons=="S_ENV_End")
rem_s <- subset(newdata, Population_comparisons=="R_ENV_Mid")
ree_s <- subset(newdata, Population_comparisons=="R_ENV_End")


```

```{r}
write_xlsx(newdata, "signif_ann_raw.xlsx")
write_xlsx(ahm_s, "signif_ann_ahm.xlsx")
write_xlsx(ahe_s, "signif_ann_ahe.xlsx")
write_xlsx(ehm_s, "signif_ann_ehm.xlsx")
write_xlsx(ehe_s, "signif_ann_ehe.xlsx")
write_xlsx(rem_s, "signif_ann_rem.xlsx")
```

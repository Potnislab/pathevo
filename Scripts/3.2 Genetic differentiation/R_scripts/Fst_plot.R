---
title: "Fst_plots"
output: html_document
---

```{r}
#Packages_required
library(ggplot2)
library(readxl)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(writexl)
library(RColorBrewer)
#getwd()
```

```{r cars}
#import file and do some cleanup:
fstraw <- readxl::read_xlsx("/path_to_fst/atDep21.nrpan.ss10.idf_fst.xlsx", sheet=2)

fst <- gather(fstraw, key = "Samples", value = "Fst", -c("CHR","Contig","Window"))

fst <- fst %>%
  mutate(Comparisons=case_when(
    Samples %in% c("7:8","9:10","11:12")~"R_ENV_Mid",
    Samples %in% c("19:20","21:22","23:24")~"R_ENV_End",
    Samples %in% c("1:2","3:4","5:6")~"S_ENV_Mid",
    Samples %in% c("13:14","15:16","17:18")~"S_ENV_End",
    Samples %in% c("1:7","4:10","5:11")~"Amb_HOST_Mid",
    Samples %in% c("13:19","16:22","17:23")~"Amb_HOST_End",
    Samples %in% c("2:8","3:9","6:12")~"Elev_HOST_Mid",
    Samples %in% c("14:20","15:21","18:24")~"Elev_HOST_End",
    TRUE ~ NA_character_
  ))

fst.nona <- na.omit(fst)
fst.nona1 <-fst[complete.cases(fst), ]
```

```{r}
orderlevel<-c("Amb_HOST_Mid", "Amb_HOST_End", "Elev_HOST_Mid", "Elev_HOST_End","S_ENV_Mid" , "S_ENV_End", "R_ENV_Mid", "R_ENV_End")
my_palette <- brewer.pal(8, "Paired")
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
clevel<-c("1:7","4:10","5:11","13:19","16:22","17:23", "2:8","3:9","6:12","14:20","15:21", "18:24","1:2", "3:4","5:6","13:14","15:16","17:18","7:8","9:10","11:12","19:20","21:22","23:24")
rclevel<- rev(clevel)

 ggplot(fst.nona, aes(x=Fst, y=factor(Samples, levels = rclevel), fill=Comparisons)) + 
  geom_boxplot(outlier.size = 0.25 ,lwd=0.25,outlier.alpha=0.5) + 
  scale_fill_manual(values = my_palette,
                    breaks = orderlevel)+
  stat_summary(fun=mean, geom="point", color="red",
     size=0.5, position=position_dodge(),show.legend = FALSE)+
  #geom_text(data = means, aes(label=value, y= value+0.08))+
  xlab("Fst")+
  ylab("Comparisons")+scale_x_continuous(breaks=seq(0,1, 0.15)) + 
  guides(fill = guide_legend(byrow = TRUE, title = NULL)) + geom_vline(xintercept = 0.15, col= "red", lty= "dotted") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.spacing.y = unit(0.3, 'cm'), legend.position='right', legend.justification='top',
    axis.text.x  = element_text(face="bold", size=12), 
       # axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
 # ggtitle("Selective pressures on Pathogen population") +
  theme(legend.position = "right"  ) + theme_bw()
ggsave("/Users/amanpreetkaur/Downloads/fst_boxplo.pdf",width = 11, height = 6, units="in", dpi=700)

#ggsave("fst_boxplot_mean.pdf", width=6, height=3.5, units=c("in"), dpi=700)

```
```{r}
#import file and do some cleanup:
fstraw <- readxl::read_xlsx("/Users/amanpreetkaur/Downloads/fst-1.xlsx", sheet=2)

fst <- gather(fstraw, key = "Samples", value = "Fst", -c("CHR","Contig","Window"))

fst <- fst %>%
  mutate(Comparisons=case_when(
    Samples %in% c("7:8")~"R_ENV_Mid_1",
     Samples %in% c("9:10")~"R_ENV_Mid_2",
     Samples %in% c("11:12")~"R_ENV_Mid_3",

    Samples %in% c("19:20")~"R_ENV_End1",
     Samples %in% c("21:22")~"R_ENV_End2",
     Samples %in% c("23:24")~"R_ENV_End3",
    
    Samples %in% c("1:2")~"S_ENV_Mid1",
    Samples %in% c("5:6")~"S_ENV_Mid2",
    Samples %in% c("3:4")~"S_ENV_Mid3",
    
    Samples %in% c("13:14")~"S_ENV_End1",
    Samples %in% c("15:16")~"S_ENV_End2",
    Samples %in% c("17:18")~"S_ENV_End3",
    
    Samples %in% c("1:7")~"Amb_HOST_Mid1",
        Samples %in% c("4:10")~"Amb_HOST_Mid2",
        Samples %in% c("5:11")~"Amb_HOST_Mid3",
    
    Samples %in% c("13:19")~"Amb_HOST_End1",
        Samples %in% c("17:23")~"Amb_HOST_End2",
        Samples %in% c("16:22")~"Amb_HOST_End3",
    
    
    Samples %in% c("3:9")~"Elev_HOST_Mid1",
    Samples %in% c("6:12")~"Elev_HOST_Mid2",
    Samples %in% c("2:8")~"Elev_HOST_Mid3",
    
    
    Samples %in% c("18:24")~"Elev_HOST_End1",
     Samples %in% c("14:20")~"Elev_HOST_End2",
     Samples %in% c("15:21")~"Elev_HOST_End3",
    
    TRUE ~ NA_character_
  ))

fst.nona <- na.omit(fst)
fst.nona1 <-fst[complete.cases(fst), ]
Fst_mean <- fst.nona1 %>% group_by(Comparisons) %>% 
  summarise(mean_fst=mean(Fst),
            median= median(Fst))
```



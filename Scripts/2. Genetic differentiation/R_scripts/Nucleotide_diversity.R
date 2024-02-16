
---
title: "Nucleotide Diversity"
output: html_notebook
---

```{r}
library(ggplot2)
library(ggpval)
library(ggpubr)
library(tidyverse) 
library(ggthemes)
```
```{r}
local_micro <- read.csv("/path_to_metapop_output/METAPOP/10.Microdiversity/local_contig_micro.csv")

meta <- read.csv("/path_to_metadata_about_samples/atdep.csv") 

local_micro_merge <- merge(local_micro, meta, by = "source", all= FALSE)

Mean_pi <- local_micro_merge %>% group_by(Sample) %>% 
  summarise(mean_pi=mean(pi),
            .groups = 'drop')

Mean_pi <- merge(Mean_pi, meta, by = "Sample")

```
```{r}
################### NUCLEOTIDE DIVERSITY  ################### 

A <- ggplot(local_micro_merge, aes(x = Condition, y = pi, fill=factor(Season, level=c('Mid','End')))) +   
  #stat_boxplot(geom ='errorbar', width = 0.5)
  geom_boxplot() + #geom_point() +
  labs( x = "", y = "Within-host Nucleotide Diversity\n for Xanthomonas perforans (pi)") +
  #scale_fill_hue(name = "Time of Sampling") + 
  scale_fill_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500")) +
  facet_wrap(~ factor(Cultivar, levels = c("Susceptible", "Resistant"))) +
  #scale_x_discrete(labels=c("Ambient"="Ambient", "Elevated O3"= "Elevated O_3"))+
  theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12),
        ) + ylim(0, 0.014) +
  #ggtitle("") +
theme_clean()
 A + theme(legend.position = "bottom")



ggsave("/Users/amanpreetkaur/Downloads/pi.pdf",width = 7, height = 7, units="in", dpi=700)  #save the plot

```{r}
my_chambers <- factor(Mean_pi$Chamber, level=c( "X1",  "X4", "X5", "X2", "X3", "X6", "E1" , "E4" , "E5", "E2" , "E3" ,  "E6")) 
library(ggh4x)
 ggplot(data = Mean_pi, aes(y = my_chambers , x= mean_pi, color = factor(Season, level=c('Mid','End')))) + 
    geom_line(aes(group = Chamber), color = "darkgrey") +
    geom_point() + 
  scale_color_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500"))  +
  facet_nested_wrap(
    vars(Cultivar, Condition),  strip.position = "right", scales = "free_y", ncol = 1,
    #axes = "all", remove_labels = "x"
  ) + theme_bw() +
    labs( x = "Average Within-host Nucleotide Diversity\n for Xanthomonas perforans (pi)", y = "Chambers") +
   theme( axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)
        ) +
   theme(legend.position = "bottom")
 #axis.text.y=element_blank(),
          #axis.ticks.y=element_blank() 
 
 ggsave("/Users/amanpreetkaur/Downloads/Average_pi.pdf",width = 7, height = 7, units="in", dpi=700)  #save the plot

  ```

   
  
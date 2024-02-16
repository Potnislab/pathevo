---
title: "Gene gain and loss plot"
output: html_document
date: "2023-10-05"
---

# We manually combined mid and end season files; based on presence and absence of genes over season.
# then we taken into account only those genes which were absent or present in  two replicates or in all the three replicates
# we also counted the genes based on each replicates in each treatment
```{r}
#packages required

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggthemes)
```

```{r}

parallel_counts <-  readxl::read_xlsx("Gene_gain_loss_counts.xlsx", sheet = 1)

  parallel_counts %>% 
    ggplot(aes(x=gain_loss, y=counts, fill = gain_loss)) + 
    geom_bar(stat = "identity", width = 0.50 )+
    #geom_col_pattern(aes(pattern_angle=genome, color=genome, fill=genome),
              #       pattern_fill="white", pattern_color="white", pattern_spacing=.03, width = 0.60 ) +
  #  scale_color_manual(values = colr) +
    scale_fill_manual(values = colr) +
   # scale_pattern_angle_manual(values = c(0, -45, 45, 80)) + #adding angle for each strip 
    labs(x="", y="Gene Count") +
    facet_grid(cultivar ~ environment) + #, switch = "x")  +
  theme_bw() + theme(legend.position="bottom")
  
```



```{r}
file <-  readxl::read_xlsx("Gene_gain_loss_counts.xlsx", sheet = 3)

file$merger <- paste(file$cultivar, file$condition)
file$genome_gain_loss <- paste(file$genome, file$gain_loss)



  file %>% 
    ggplot(aes(x=genome_gain_loss, y= counts, fill = common)) + 
    geom_bar(position = position_dodge2(preserve = "single") , stat = "identity", width = 0.50 )+
    #geom_col_pattern(aes(pattern_angle=genome, color=genome, fill=genome),
              #       pattern_fill="white", pattern_color="white", pattern_spacing=.03, width = 0.60 ) +
  #  scale_color_manual(values = colr) +
   # scale_fill_manual() +
   # scale_pattern_angle_manual(values = c(0, -45, 45, 80)) + #adding angle for each strip 
    labs(x="", y="Gene Count") +
    facet_grid(merger ~ gain_loss,scales = "free_x")  +
  theme_bw() + theme(legend.position="bottom") +
    theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
  #ggtitle("") +
  theme(legend.position = "bottom"  ) 
  
#ggsave("/Users/amanpreetkaur/Downloads/genes_counts.pdf",width = 7, height = 7, units="in", dpi=700)
  
```
```{r}


file <-  readxl::read_xlsx("Gene_gain_loss_counts.xlsx", sheet = 4)
file$merge <- paste(file$Cultivar, file$Condition)

file %>% 
  ggplot(aes(x= factor(cog, levels = c("NA", "S", "U", "NU", "G", "GM", "IQ", "M", "L", "P", "K", "C", "F", "E")), y= value, fill = gain_loss)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.75 ) +
  #geom_col_pattern(aes(pattern_angle=genome, color=genome, fill=genome),
  #       pattern_fill="white", pattern_color="white", pattern_spacing=.03, width = 0.60 ) +
  #  scale_color_manual(values = colr) +
  # scale_fill_manual() +
  # scale_pattern_angle_manual(values = c(0, -45, 45, 80)) + #adding angle for each strip 
  labs(x="", y="Gene Count") +
  scale_fill_manual(name = "Gain/Loss", values=c(  "#00823c","#CE74A7"))+
  theme_bw() +  
  facet_grid( ~ merge, cols = NULL,  scales = "free_x")  +
  theme_bw() + theme(legend.position="bottom") +
  theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) + coord_flip() +
  #ggtitle("") +
  theme(legend.position = "bottom") 

#ggsave("/Users/amanpreetkaur/Downloads/genes_gain_loss.pdf",width = 9, height = 7, units="in", dpi=700)

```
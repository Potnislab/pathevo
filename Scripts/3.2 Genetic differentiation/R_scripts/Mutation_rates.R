---
title: "Mutation rates"
output: html_document
date: "2023-06-11"
---

```{r}
setwd("/path_to_rhometa_output/Theta_Est_nrpan_Atdep_blk_v2/RHOMETA/theta_estimate")

files <- list.files(pattern="\\.csv$")
 
```

```{r}

#compling all the files together in one file

library(tidyverse)

data <- map_dfr(
  .x = files,
  ~ read.table(.x) %>%
    mutate(
      filename = stringr::str_replace(
        .x,
        "\\.csv$",
        ""
      )
    )
)
library(stringr)
 
# Split name column into firstname and last name

data[c('Theta', 'estimate')] <- str_split_fixed(data$V1, ',', 2)

data[c('Sample', 'extra')] <- str_split_fixed(data$filename, "_Theta_estimate_stats", 2 )
 
data <- data[, c(5,4 ,3)]

is.remove <- c("mean_depth", "median_depth", "tps_median_depth")

data_file <- data[!(data$Theta %in% is.remove),]

data_file <- data_file[, -c(3)]
colnames(data_file) <- c("source", "Theta")

write.csv(data_file, file = "/Save_to_here/Atdep/RHOMETA/Theta_blk_estimate.csv")


```

```{r}
data_file <- read.csv("/path_to_input/RHOMETA/Theta_blk_estimate.csv")

meta <- read.csv("/Path_to_metdata/AtDep/atdep.csv")

theta <- merge(data_file, meta, by = "source", all= FALSE)
theta$Theta <- as.numeric((theta$Theta))

```


```{r}
###################  BoXplot of  MUTATION RATES      ########################

ggplot(theta, aes(x = Condition, y = Theta, fill=factor(Season, level=c('Mid','End')))) +   
  #stat_boxplot(geom ='errorbar', width = 0.5)
  geom_boxplot() + #geom_point() +
  labs( x = "", y = "Mean Population Mutation Rate per site") +  ylim(0, 0.0018) +
  #scale_fill_hue(name = "Time of Sampling") + 
  scale_fill_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500"))+
  theme_clean() + 
   facet_wrap(~ factor(Cultivar, levels = c("Susceptible", "Resistant"))) +
  #scale_x_discrete(labels=c("Ambient"="Ambient", "Elevated O3"= "Elevated O_3"))+
  theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
  #ggtitle("") +
  theme(legend.position = "bottom"  ) 



ggsave("/Users/amanpreetkaur/Downloads/theta.pdf",width = 7, height = 7, units="in", dpi=700)
```


```{r}
Mean_theta <- theta %>% group_by(Sample) %>% 
  summarise(mean_theta=mean(Theta),
            .groups = 'drop')

Mean_theta <- merge(Mean_theta, meta, by = "Sample")

```

```{r}
###################  lineplot of  MUTATION RATES      ########################

my_chambers <- factor(Mean_theta$Chamber, level=c( "X1",  "X4", "X5", "X2", "X3", "X6", "E1" , "E4" , "E5", "E2" , "E3" ,  "E6")) 

 ggplot(data = Mean_theta, aes(y = my_chambers , x= mean_theta, color = factor(Season, level=c('Mid','End')))) + 
    geom_line(aes(group = Chamber), color = "darkgrey") +
    geom_point() + 
  scale_color_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500"))  +
  facet_nested_wrap(
    vars(Cultivar, Condition),  strip.position = "right", scales = "free_y", ncol = 1,
    #axes = "all", remove_labels = "x"
  ) + theme_bw() +
    labs( x = "Average mutation rates\n for Xanthomonas perforans (theta)", y = "Chambers") +
   theme( axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)
        ) +
   theme(legend.position = "bottom")
 
 ggsave("/Users/amanpreetkaur/Downloads/theta_average.pdf",width = 7, height = 7, units="in", dpi=700)

```

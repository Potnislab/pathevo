# Nucleotide DIversity plot

# libraried required

    install.packges("ggplot2")
    library(ggplot2)
    install.packges("ggpval")
    library(ggpval)
    install.packges("ggpubr")
    library(ggpubr)
    install.packges("tidyverse")
    library(tidyverse) 

import or read the input files into Rstudio

    local_micro <- read.csv("/path/to/metapop_pi_diversity_local.csv")
    theta <- read.csv("/path/to/rhometa/THETA.csv")
    #separating the theta file according to the ambient and elevated ozone conditions
     theta1 <- theta[theta$Environment != "Ambient", ]
     theta2 <- theta[theta$Environment != "Elevated O_3", ]

# code for plot

      ggplot(local_micro, aes(x = Cultivar, y = pi, fill=factor(Time, level=c('Mid-season','End-season')))) +   
      geom_boxplot() + #geom_point() +
      labs( x = "", y = "Within-host Nucleotide Diversity") +
      #scale_fill_hue(name = "Time of Sampling") + 
      scale_fill_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500")) +
      theme_bw() + 
      facet_wrap(~Environment) +
      theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
       #ggtitle("") +
      theme(legend.position = "bottom" ) 
      
# save the plot 
     ggsave("/Users/amanpreetkaur/Downloads/pi.pdf",width = 7, height = 7, units="in", dpi=700)  #save the plot



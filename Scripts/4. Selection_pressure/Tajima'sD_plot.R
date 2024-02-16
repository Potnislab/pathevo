---
title: "Tajima's D Plot"

---

## set up data
#We can read in the annotation file in R using the following code:
annotate <- read.csv("path_to_/EGGNOG_annotations.csv")
remove_last_n <- 2
annotate$query_1 <-  substr(annotate$query, 1, nchar(annotate$query) - remove_last_n)
annotate <- annotate[, c(22,7,8,9,21)]
colnames(annotate) <- c("contig_gene"  ,    "COG_category"  , "Description"  ,  "Preferred_name", "PFAMs"   )
#annotate
#rownames(annotate) <- annotate[, c(1)]



# Load data
local_gene <- read.table("/path_to/local_gene_microdiversity.tsv ", header = TRUE)
metadata <- read.csv("/path_to/Atdep_metadata.csv")
metadata$source <- factor(metadata$source) #combining the metadata with the orginal file
local_gene <- merge(local_gene, metadata, by = "source")
local_gene$Treatment <- paste(local_gene$Merger,"-",local_gene$Season)

local_gene_all <- merge(local_gene, annotate, by= "contig_gene", all.x = TRUE)

local_gene <- local_gene_all[local_gene_all$pNpS_ratio != "Inf", ] 

#write.csv(local_gene, file = "/Users/amanpreetkaur/Downloads/Localgene.csv")
```


#tajima's_D
B <- ggplot(data = local_gene_all, aes(x = Condition , y = taj_D, color = factor(Season, level=c('Mid','End')))) + 
    #geom_line(aes(group = Chamber), color = "darkgrey") +
    geom_jitter(position = position_jitter(height = .25, width = .25), size = 1,  na.rm = TRUE) + 
  #geom_boxplot() +
  scale_color_manual(name = "Time of Sampling", values=c( "#4169E1", "#FFA500"))  +
   facet_wrap(~ factor(Cultivar, levels = c("Susceptible", "Resistant"))) +
  #scale_x_discrete(labels=c("Ambient"="Ambient", "Elevated O3"= "Elevated O_3"))+
  theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12),
        ) + 
  #ggtitle("") +
#theme_clean() + 
  theme_bw()+
  theme(legend.position = "bottom")
B
ggsave("/Users/amanpreetkaur/Downloads/tajima'D.pdf", width=7, height=4, units=c("in"), dpi=700)

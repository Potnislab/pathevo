# plot for pN/pS ratio and Tajima's D 
    
    File <- read.csv("/path/to/Metapop/10.Microdiversity/local_gene.csv") #read the file
    metadata <- read.csv("/path/to/metadata.csv")
    metadata$source <- factor(metadata$source) #combining the metadata with the orginal file
    local_gene2 <- merge(File, metadata, by = "source")


#  Libraries  required

    library(ggplot2)
    
# import the input files
removing all the genes have pN/pS value =inf and choosing the genes those have greater than 1 pN/pS ratio and less than -1 Tajima's D values

    local_gene <- local_gene2[local_gene2$pNpS_ratio != "Inf", ] 
    pNpS_taj <- subset(local_gene, pNpS_ratio > 1 & taj_D < 0, select = c(2:19))
    #write.csv(/save/to/pNpS_taj, file= "local_Pn_taj.csv") #to save the csv file

# plot for pN/pS ratio

saving the plot in "A" vector

    A <- ggplot(selection, aes( x= gene_no.,y = pNpS_ratio, group=season)) +  
    geom_jitter(aes(color = factor(season))) + facet_grid(cultivar ~ condition) +
    labs( x = "", y = "pNpS") +
    scale_color_manual(name = "Time of Sampling", values =   c( "#FFA500", "#4169E1"))+
    theme_bw() +
    theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
    #ggtitle("") +
    theme(
      legend.position = "none",
    ) 

# plot for tajima'D value

saving the plot in "B" vector

    B <- ggplot(selection, aes( x= gene_no.,y = taj_D, group=season)) + 
    geom_jitter(aes(color = factor(season))) + facet_grid(cultivar ~ condition) +
    labs( x = "gene_id", y = "taj_D") +
    scale_color_manual(name = "Time of Sampling", values =  c( "#FFA500", "#4169E1"))+
    theme_bw() + 
    theme(axis.text = element_text(face="bold", size=12), 
        axis.title.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=12),
        title =element_text(face="bold", size=12),
        legend.title = element_text(face="bold", size=12),
        legend.text=element_text(face = "bold", size=12),
        strip.text.x = element_text(face="bold", size=12)) +
    #ggtitle("") +
    theme(legend.position = "bottom"  ) 

Arranging both the plots (A & B) in one figure as two rows and one column
    
    ggarrange(A, B, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

Save the plot 

    ggsave("/save/to/Selection.pdf",width = 6, height = 8, units="in", dpi=700) #to save the file



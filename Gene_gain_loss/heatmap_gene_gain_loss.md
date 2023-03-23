#  libraries required
        library("circlize")
        library("RColorBrewer")
        library("RColorBrewer")
        library(dendextend)
        library(ComplexHeatmap)
        library(ggplot2)
        #remotes::install_github("coolbutuseless/ggpattern")
        library(ggpattern)

#  import all files
        gene_presabs <- read.csv("/path/to/gene_gain_loss.csv", row.names = "Type") #read the csv file of gene gain and loss with row_names of annotations or gene_ids
        atdep <- read.csv("/path/to/metadata.csv", row.names = "Sample") #file containing metadata information about all samples

        gain_loss <- read.csv("/path/to/gain_loss_counts.csv") #the file containing all the number of genes counts gained or lost from the samples

```

```{r}

#heatmap

gene_filter <- gene_presabs[, -c(1:5)] #removed the column that was already imported as row names 
mat <- as.matrix(gene_filter[,c(1:24)]) #c is specifying the column numbers 

#to divide the heatmap acording to the different conditions

level_order <- factor(atdep$Merger, level =c("ECW_A", "X10R_A", "ECW_E", "X10R_E")) #ordering the column for splicing

col_fun = colorRamp2(c( 0, 1), c("royalblue", "red2")) #color for the cells

row_ha = rowAnnotation(foo2 = gene_presabs$old) #annotating the rows according to the genome of different strains

pdf("/path/to/save/heatmap.pdf")
Heatmap(mat, name = "key", 
        rect_gp = gpar(col = "black", lwd = 1),
        row_split = gene_presabs$old,
        row_names_gp = gpar(fontsize = 3),
      #row_order = order(as.numeric(gsub("row", "", rownames(gene_filter)))), 
        column_order = order(as.numeric(gsub("column", "", colnames(gene_filter)))),
      # column_split = factor(rep(LETTERS[1:4], 6), levels = LETTERS[1:4]),
      column_split = level_order,
        show_column_dend = FALSE, row_names_side = "left", show_row_dend = FALSE,
       cluster_row_slices = FALSE, 
       cluster_column_slices = FALSE,
      right_annotation =  row_ha,  col =  col_fun
                                              )

dev.off()

```
```{r}
#figure for gene gain/loss counts

colr <- c( "#00823c","#CE74A7", "steelblue", "#9e9e9e")  #selecting colours
  gain_loss %>% 
    ggplot(aes(x=Gain_loss, y=negative), pattern="stripe") + 
    geom_col_pattern(aes(pattern_angle=genome, color=genome, fill=genome),
                     pattern_fill="white", pattern_color="white", pattern_spacing=.03, width = 0.60 ) +
    scale_color_manual(values = colr) +
    scale_fill_manual(values = colr) +
    scale_pattern_angle_manual(values = c(0, -45, 45, 80)) + #adding angle for each strip 
    labs(x="", y="Gene Count") +
    facet_grid(Cultivar ~ Condition) + #, switch = "x")  +
  theme_bw() + theme(legend.position="bottom")
  
ggsave("/path/to/genepresabs_counts_strips.pdf", width = 5, height = 5, units="in", dpi=800) #to save the plot


```





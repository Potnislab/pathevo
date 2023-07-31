---
title: "Figure 1"
output: html_notebook
---

```{r}
library(readxl)
library(dplyr)
library(ggplot2)
library(circlize)
```


```{r}
#Figure 1.B 

ozone_values <- read_excel("Ozone Values.xlsx")
ozone_values$environment <- ifelse(ozone_values$chamber== 2 | ozone_values$chamber== 3 | ozone_values$chamber== 6 | ozone_values$chamber== 7 | ozone_values$chamber== 10 | ozone_values$chamber== 11,"Elevated Ozone", "Ambient")
ozone_values$cycle_date <- as.Date(ozone_values$cycle_date, origin = "1899-12-30")
ozone_values$season <- if_else(ozone_values$cycle_date < '2021-06-18', "Mid", "End")
ozone_values <- subset(ozone_values, ozone_values$treatment == "Inoculated")
ozone_values_growing_season <- subset(ozone_values, ozone_values$cycle_date < '2021-07-29')


ggplot(ozone_values_growing_season, aes(x = as.Date(cycle_date), y = cycle_o3_ave, color = environment)) + geom_point(size = 0.1) + geom_smooth(method = "lm", se = FALSE) + scale_color_manual(values = c("Ambient" = "lightskyblue", "Elevated Ozone" = "red2")) + xlab("") + ylab("Ozone (O3) Concentration (ppm)") + geom_hline(yintercept=40, linetype='dotted') + labs(color ="Environment")
```

```{r}
#Figure 1.C

abund_severity <- read_excel("disease_severity.xlsx")
abund_severity$genotype <- factor(abund_severity$genotype, levels = c("AL65", "AL22", "SC4"))
abund_severity$season <- factor(abund_severity$season, levels = c("Mid", "End"))
abund_severity_resistant <- subset(abund_severity, abund_severity$cultivar == "Resistant")

data_vline <- data.frame(environment = c("Ambient", "Ambient", "Elevated", "Elevated"), vline = c(0.291, 0.375, 2.125, 12.614), season = c("End", "Mid", "End", "Mid"))
data_vline$season <- factor(data_vline$season, levels = c("Mid", "End"))

ggplot(abund_severity_resistant, aes(x = disease_severity, y = genotype_absolute_abundance, color = genotype, shape = season)) + geom_point()  + scale_color_manual(values = c("#CE74A7", "#00823c", "#9e9e9e")) + scale_shape_manual(values = c(1, 2)) + facet_grid(environment~cultivar) + theme_bw() + geom_vline(data = data_vline, aes(xintercept = vline, linetype = season), size = 0.2) + scale_linetype_manual(values = c(4, 1)) + xlab("Disease Severity") + ylab("Genotype Absolute Abundance") + labs(shape = "Season") + labs(linetype = "Season") + labs(color = "Genotype")
```

```{r}
#Figure 1.D

##Circular Chord Diagram

relative_abundance_for_chord <- read_excel("relative_abundance_for_chord.xlsx")

grid.col = c(SC4 = "#9e9e9e", AL22 = "#00823c", AL65 = "#CE74A7", AEM = "#4169E1", AEE = "#FFA500", AXM = "#4169E1", AXE = "#FFA500", OEM = "#4169E1", OEE = "#FFA500", OXM = "#4169E1", OXE = "#FFA500")

chord <- chordDiagram(relative_abundance_for_chord, order = c("SC4", "AL22", "AL65", "AEM", "AEE", "AXM", "AXE", "OEM", "OEE", "OXM", "OXE"), small.gap = 5, big.gap = 10, grid.col = grid.col, row.col = c("#FF000080", "#00FF0010", "#0000FF10"), annotationTrack = c("grid"), preAllocateTracks = list(list(tAbundanceck.height = 0.1), list(track.height = 0.1)), annotationTrackHeight = mm_h(3))

for(si in get.all.sector.index()) {
  if (endsWith(si, "M")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "Mid", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)}
  if (endsWith(si, "E")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "End", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)}
  if (endsWith(si, "5")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "SC6", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)
  }
    if (endsWith(si, "2")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "SC3", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)
    }
    if (endsWith(si, "4")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "SC4", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)
  }
  }

highlight.sector(c("AEM", "AEE"), track.index = 2, col = "#bf6236", text = "Susceptible", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("AXM", "AXE"), track.index = 2, col = "#9d4ee7", text = "Resistant", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("OEM", "OEE"), track.index = 2, col = "#bf6236", text = "Susceptible", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("OXM", "OXE"), track.index = 2, col = "#9d4ee7", text = "Resistant", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("AEM", "AEE", "AXM", "AXE"), track.index = 1, col = "#87CEFA", text = "Ambient", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("OEM", "OEE", "OXM", "OXE"), track.index = 1, col = "#FF4e4e", text = "Elevated O3", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

circos.clear()



##Boxplot

relative_abundance_for_boxplot <- read_excel("relative_abundance_for_boxplot.xlsx")
relative_abundance_for_boxplot$genotype <- factor(relative_abundance_for_boxplot$genotype, levels = c("AL65", "AL22", "SC4"))
relative_abundance_for_boxplot$index <- factor(relative_abundance_for_boxplot$index, levels = c("ASM", "ASE", "ARM", "ARE", "ESM", "ESE", "ERM", "ERE"))
relative_abundance_for_boxplot$season <- factor(relative_abundance_for_boxplot$season, levels = c("Mid", "End"))

ggplot(relative_abundance_for_boxplot, aes(x = season, y = relative_abundance, fill = genotype)) + geom_boxplot(lwd=0.2) + scale_fill_manual(values = c("#CE74A7", "#00823c", "#9e9e9e")) + facet_grid(cultivar~environment) + theme_bw() + xlab("Season") + ylab("Relative Abundance") + labs(fill = "Genotype")
```

```{r}
#Figure 1.E

##Circular Chord Diagram

absolute_abundance_for_chord <- read_excel("absolute_abundance_for_chord.xlsx")

grid.col = c(SC4 = "#9e9e9e", AL22 = "#00823c", AL65 = "#CE74A7", AEM = "#4169E1", AEE = "#FFA500", AXM = "#4169E1", AXE = "#FFA500", OEM = "#4169E1", OEE = "#FFA500", OXM = "#4169E1", OXE = "#FFA500")

chord <- chordDiagram(absolute_abundance_for_chord, order = c("SC4", "AL22", "AL65", "AEM", "AEE", "AXM", "AXE", "OEM", "OEE", "OXM", "OXE"), small.gap = 5, big.gap = 10, grid.col = grid.col, row.col = c("#FF000080", "#00FF0010", "#0000FF10"), annotationTrack = c("grid"), preAllocateTracks = list(list(tAbundanceck.height = 0.1), list(track.height = 0.1)), annotationTrackHeight = mm_h(3))

for(si in get.all.sector.index()) {
  if (endsWith(si, "M")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "Mid", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)}
  if (endsWith(si, "E")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "End", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)}
  if (endsWith(si, "5")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "SC6", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)
  }
    if (endsWith(si, "2")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "SC3", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)
    }
    if (endsWith(si, "4")) {
    xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
    ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circos.text(mean(xlim), mean(ylim), "SC4", sector.index = si, track.index = 3, facing = "bending.inside", niceFacing = TRUE, col = "black", cex = 0.9)
  }
  }

highlight.sector(c("AEM", "AEE"), track.index = 2, col = "#bf6236", text = "Susceptible", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("AXM", "AXE"), track.index = 2, col = "#9d4ee7", text = "Resistant", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("OEM", "OEE"), track.index = 2, col = "#bf6236", text = "Susceptible", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("OXM", "OXE"), track.index = 2, col = "#9d4ee7", text = "Resistant", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("AEM", "AEE", "AXM", "AXE"), track.index = 1, col = "#87CEFA", text = "Ambient", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

highlight.sector(c("OEM", "OEE", "OXM", "OXE"), track.index = 1, col = "#FF4e4e", text = "Elevated O3", niceFacing = TRUE, facing="bending.inside", font = 1, cex = 0.9)

circos.clear()



##Boxplot

absolute_abundance_for_boxplot <- read_excel("absolute_abundance_for_boxplot.xlsx")
absolute_abundance_for_boxplot$genotype <- factor(absolute_abundance_for_boxplot$genotype, levels = c("AL65", "AL22", "SC4"))
absolute_abundance_for_boxplot$index <- factor(absolute_abundance_for_boxplot$index, levels = c("ASM", "ASE", "ARM", "ARE", "ESM", "ESE", "ERM", "ERE"))
absolute_abundance_for_boxplot$season <- factor(absolute_abundance_for_boxplot$season, levels = c("Mid", "End"))

ggplot(absolute_abundance_for_boxplot, aes(x = season, y = genotype_absolute_abundance, fill = genotype)) + geom_boxplot(lwd=0.2) + scale_fill_manual(values = c("#CE74A7", "#00823c", "#9e9e9e")) + facet_grid(cultivar~environment) + theme_bw() + xlab("Season") + ylab("Absolute Abundance") + labs(fill = "Genotype")
```


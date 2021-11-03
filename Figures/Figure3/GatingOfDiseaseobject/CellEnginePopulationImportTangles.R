##############################################################################
# Script:CellEnginePopulationImport.R
# Project: HiRes AD Scanned Data
# Author: Dmitry Tebaykin
# Date : 19th August 2020
# This program reads in a folder of tsv files, adds a filename column, merges the files, writes out a csv & makes heatmaps

# Script: merge_table.R
# Project: brain
# Author: David Glass
# Date: 8-17-20
# From program that reads in a folder of tsv files, adds a filename column, merges the files, and writes out a csv
##############################################################################

### LIBRARIES ###
libs <- c('ggplot2','RColorBrewer','reshape2','devtools',
          'pheatmap','flowCore','plyr','scales','coin','grid',
          'gridExtra','tidyr','randomcoloR','viridis','tidyverse', 'ComplexHeatmap')
for (L in libs){
  if(!require(L, character.only = T)) {
    install.packages(L, dep = T, quiet = T)
    
    if (!require(L, character.only = T)) {
      if (L %in% c('flowCore', 'flowWorkspace')) {
        install_github("RGLab/flowCore", ref="trunk")
      } else {
        suppressWarnings(BiocManager::install(L))
      }
    }
  }
}

library(readr)
library(readr)
library(scales)
require(data.table)
setwd("/Volumes/KausaliaHD/R4R")
source("/Volumes/KausaliaHD/R4R/General/readWriteCytofFCS.R")

### MAIN ####

## read in files
mainPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/CellEngineCurated/EzDeepcellPooled/" 

## read in Pooleddata "EzDeepcellPooled/EzDeepCellDataNFT/TotalTauPHF1Tau/" 
## read in EzSegdata  "/" 
dataPath = paste0(mainPath, "alldataDiseasePooled_TotalTauvsPHF1_GatedPopulations/PHF1Tau_Grades/") 

outputFolder = paste0(dataPath, "TangleThread_TotalTauvsPHF1Tau/")
dir.create(outputFolder, showWarnings= F, recursive = T)

files <- list.files(dataPath, pattern="tsv", full.names=FALSE, recursive=FALSE)
df.list <- list()
for (file in files) {
  df.list[[file]] <- fread(paste0(dataPath, file)) %>%
    .[, file:=file]
}
dat <- rbindlist(df.list)
rm(df.list)

## correct the names
replaceNames <- function(vec) {
  # replaces vec format "Y89Di (CD45)" with "CD45"
  if (!grepl("\\(", vec)) return(vec)
  start <- gregexpr("\\(", vec)[[1]][[1]] + 1
  stop <- gregexpr("\\)", vec)[[1]][[1]] - 1
  return(substr(vec, start, stop))
}

setnames(dat, colnames(dat), unlist(lapply(colnames(dat), replaceNames)))

## if you want to write out the merged data as a csv:
fwrite(dat, paste0(outputFolder, "TotalTauvsPHF1TauGrades.csv"))

## if you want to convert the merged data to data.frame
df <- as.data.frame(dat)
rm(dat)

## to read in the csv of CellEngine 
#df = read_csv (paste0(dataPath,""TotalTauvsPHF1TauGrades.csv"")) # "EzSegHiLoPHF1Tau.csv" or "EZSegDeepCellHiLoPHF1Tau.csv"

#to make plot Scatterplot (No breaklist)
df$file = as.factor(df$file)
table(df$file)
plot_data=df
plot_data = subset(plot_data, obj_size<10000)
plot_data = plot_data[order(plot_data$obj_size),]
ggplot(plot_data, aes(x=asinh(PHF1Tau), y=asinh(TotalTau))) +
  #stat_density_2d(bins = 8000, alpha = 0.005, geom = "polygon") + # after_stat(level)) (estimates density of point in area), geom = "polygon"  bins = 50,
  #scale_colour_gradientn(colours = heat.colors(5), values = breakList) +
  #color-blind friendly
  scale_colour_gradientn(colours = c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  #scale_color_viridis_c(direction = 1, option = "B", values = breakList) +
  geom_point(aes(color = obj_size), size = 1.5, alpha = 1) +
  #ylim(0,1.5) +
  #xlim(0,1.6)+
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10") +
  theme_classic()
ggsave(paste0(outputFolder, "/ScatterPlotTotalTauvsPHF1Tauvsobjsize(colorblind_arcshinh)A2.pdf") , height = 10, width = 15)

#to make plot Scatterplot (breaklist)
df$file = as.factor(df$file)
table(df$file)
plot_data=df
plot_data = subset(plot_data, obj_size<10000)
plot_data = plot_data[order(plot_data$obj_size),]
breakList = c(seq(0.0,0.08,0.01),seq(0.081,1.5,0.5))
ggplot(plot_data, aes(x=asinh(PHF1Tau), y=asinh(TotalTau))) +
  #stat_density_2d(bins = 8000, alpha = 0.005, geom = "polygon") + # after_stat(level)) (estimates density of point in area), geom = "polygon"  bins = 50,
  #scale_colour_gradientn(colours = heat.colors(5), values = breakList) +
  #color-blind friendly
  scale_colour_gradientn(colours = c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"), values = breakList) +
  #scale_color_viridis_c(direction = 1, option = "B", values = breakList) +
  geom_point(aes(color = obj_size), size = 3, alpha = 1) +
  #ylim(0,1.5) +
  #xlim(0,1.6)+
  #scale_x_continuous(trans = "log10") +
  #scale_y_continuous(trans = "log10") +
  theme_classic()
ggsave(paste0(outputFolder, "/ScatterPlotTotalTauvsPHF1Tauvsobjsize(colorblind_arcshinh)B3.pdf") , height = 10, width = 15)


## to convert file id to gatenames
df$file = as.factor(df$file)
table(df$file)
levels(df$file) = c("HiPHF1Tau", "LoPHF1Tau", "MedPHF1Tau") 
table(df$file)
df$file = ordered(df$file, levels = c("LoPHF1Tau","MedPHF1Tau","HiPHF1Tau"))  
table(df$file)

## to make HeatMaps with Rocket or Viridis Colours
callback=function(x,data){
  hc=hclust(dist(data), "ward.D2")
  dend=reorder(as.dendrogram(hc),wts=data)
  as.hclust(dend,"ward.D2")
}
set.seed(240214)

## Basic heatmap. Cluster_cols = F to not cluster the brain regions

## to plot All markers
allMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","Reelin","MAP2","X8OHGuano",
               "PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau",
               "Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1NTF","pTDP43","EEA1","CD56Lyo","CD47", "TotalTau", "PanGAD6567","CD33Lyo")
commonCols = c(allMarkers, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_AllMarkers1.pdf"), width = 8, height = 5, onefile = FALSE)
    breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
    colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
    p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "All Population, Mean Intensities, scaled",
               col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
    dev.off()
    dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_AllMarkers2.pdf"), width = 8, height = 5, onefile = FALSE)
    colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
    p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "All Population, Mean Intensities, scaled",
               col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
    dev.off()
    dev.off()

## to plot Gated markers
markersGated = c("PHF1Tau") #,"PHF1Tau","ApoE4"
commonCols = c(markersGated, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  = plot_data[, match(c("PHF1Tau"),colnames(plot_data))]  

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_GatesMarkers1.pdf"), width = 6, height = 5, onefile = FALSE)
breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
p=pheatmap((plot_data), cluster_rows = F, cluster_cols = F, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Gated Population, Mean Intensities, scaled",
           col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_GatedMarkers2.pdf"), width = 6, height = 5, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
p=pheatmap((plot_data), cluster_rows = F, cluster_cols = F, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Gated Population, Mean Intensities, scaled",
           col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
dev.off()
dev.off()

## to plot Injury markers
injuryMarkers = c("MFN2","PolyubiK48","PolyubiK63","pTDP43","EEA1","X8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4")
commonCols = c(injuryMarkers, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  =plot_data[, match(c("Amyloidbeta140","Amyloidbeta142","PanAmyloidbeta1724","ApoE4","X8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63","EEA1"),colnames(plot_data))] 

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_InjuryMarkers1.pdf"), width = 6, height = 5, onefile = FALSE)
    breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
    colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
    p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Injury Population, Mean Intensities, scaled",
              col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
    dev.off()
    dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_InjuryMarkers2.pdf"), width = 6, height = 5, onefile = FALSE)
    colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
    p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Injury Population, Mean Intensities, Scaled",
                col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
    dev.off()
    dev.off()

## to plot Neuronal markers
neuronalMarkers = c("MAP2","VGAT","VGLUT1","VGLUT2","Reelin","Calbindin","Calretinin","Parvalbumin","Presenilin1NTF","PanGAD6567","CD47","PSD95", "Synaptophysin")
commonCols = c(neuronalMarkers, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  = plot_data[, match(c("MAP2","VGLUT1","VGLUT2","VGAT","PanGAD6567","Reelin","Presenilin1NTF","Parvalbumin","Calbindin","Calretinin","CD47"),colnames(plot_data))]

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_NeuronalMarkers1.pdf"), width = 6, height = 5, onefile = FALSE)
    breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
    colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
    p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Neuronal Population, Mean Intensities, Scaled",
               col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
    dev.off()
    dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_NeuronalMarkers2.pdf"), width = 6, height = 5, onefile = FALSE)
    colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
    p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Neuronal Population, Mean Intensities, Scaled",
               col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
    dev.off()
    dev.off()

## to plot Neuronal Structure (soma/dendrites vs axonal)
neuronalStructure = c("MAP2","MBP","MAG")
commonCols = c(neuronalStructure, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  = plot_data[, match(c("MAP2","MBP","MAG"),colnames(plot_data))]

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_NeuronalStructure1.pdf"), width = 6, height = 5, onefile = FALSE)
  breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
  colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Neuronal Neuronal Structure (soma/dendrites vs axonal), Mean Intensities, scaled",
             col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_NeuronalStructure2.pdf"), width = 6, height = 5, onefile = FALSE)
  colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "Neuronal Neuronal Structure (soma/dendrites vs axonal), Mean Intensities, scaled",
             col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()

## to plot NonNeuronal markers
nonNeuronalMarkers = c("Iba1","CD45","GFAP","CD31","MCT1","CD105","CD33Lyo")
commonCols = c(nonNeuronalMarkers, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  = plot_data[, match(c("Iba1","CD45","CD33Lyo","GFAP","CD31","CD105","MCT1"),colnames(plot_data))]

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_NonNeuronalMarkers1.pdf"), width = 6, height = 5, onefile = FALSE)
  breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
  colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "NonNeuronal Population, Mean Intensities, scaled",
             col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap_scaled_NonNeuronalMarkers2.pdf"), width = 6, height = 5, onefile = FALSE)
  colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "NonNeuronal Population, Mean Intensities,scaled",
             col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()
  
  ## to plot Glia markers
  nonNeuronalMarkers = c("Iba1","CD45","GFAP","CD33Lyo")
  commonCols = c(nonNeuronalMarkers, "file")
  df_subset = df[,which(colnames(df) %in% commonCols)]
  plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
  plot_data$file = NULL
  rownames(plot_data) =  plot_data$Group.1
  plot_data$Group.1 = NULL
  plot_data  = plot_data[, match(c("Iba1","CD45","CD33Lyo","GFAP"),colnames(plot_data))]
  
  # this maps  Viridis color palette
  pdf(paste0(outputFolder,"/Heatmap_scaled_GliaMarkers1.pdf"), width = 6, height = 5, onefile = FALSE)
  breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
  colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "NonNeuronal Population, Mean Intensities, scaled",
             col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()
  
  # this maps Rocket color palette
  pdf(paste0(outputFolder,"/Heatmap_scaled_GliaMarkers2.pdf"), width = 6, height = 5, onefile = FALSE)
  colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "NonNeuronal Population, Mean Intensities,scaled",
             col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()
  
  ## to plot Vessel markers
  nonNeuronalMarkers = c("CD31","MCT1","CD105")#"CD45",
  commonCols = c(nonNeuronalMarkers, "file")
  df_subset = df[,which(colnames(df) %in% commonCols)]
  plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
  plot_data$file = NULL
  rownames(plot_data) =  plot_data$Group.1
  plot_data$Group.1 = NULL
  plot_data  = plot_data[, match(c("CD31","CD105","MCT1"),colnames(plot_data))]#"CD45",
  
  # this maps  Viridis color palette
  pdf(paste0(outputFolder,"/Heatmap_scaled_VesselOnlyMarkers1.pdf"), width = 6, height = 5, onefile = FALSE)
  breakList = c(0.0,0.05,0.10,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39, 0.4,0.44,0.46)
  colfunc <- colorRampPalette(c("#440154FF","#471164FF","#481F70FF","#472D7BFF","#443A83FF","#404688FF","#3B528BFF","#365D8DFF","#31688EFF","#2C728EFF","#287C8EFF","#24868EFF","#21908CFF","#1F9A8AFF","#20A486FF","#27AD81FF","#35B779FF","#47C16EFF","#5DC863FF","#75D054FF","#8FD744FF","#AADC32FF","#C7E020FF","#E3E418FF","#FDE725FF"))
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "NonNeuronal Population, Mean Intensities, scaled",
             col= colfunc(length(breakList)),cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()
  
  # this maps Rocket color palette
  pdf(paste0(outputFolder,"/Heatmap_scaled_VesselonlyMarkers2.pdf"), width = 6, height = 5, onefile = FALSE)
  colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#FFFFFF" White,"#C62828" Red
  p=pheatmap((plot_data), cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "column", main = "NonNeuronal Population, Mean Intensities,scaled",
             col = colfunc(100), direction = 1, cellwidth = 10, cellheight = 10)
  dev.off()
  dev.off()




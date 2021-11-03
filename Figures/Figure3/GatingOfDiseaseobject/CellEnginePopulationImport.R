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
          'gridExtra','tidyr','randomcoloR','viridis','tidyverse')
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
mainPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/CellEngineCurated/"  # change this to the path to the folder with your tsv files

## read in Pooleddata "EzDeepcellPooled/EzDeepCellDataNFT/TotalTauPHF1Tau/" or "EzDeepcellPooled/EzDeepCellDataNFT/EzDeepCellDataAmyloidPlaques/Abeta40Abeta42/" 
## read in EzSegdata  "EzDataSegOnly/EzDataOnlyAmyloidPlaques/TotalTauPHF1Tau/" or "EzDataSegOnly/EzDataOnlyAmyloidPlaques/Abeta40Abeta42/" 
dataPath = paste0(mainPath,"EzDeepcellPooled/EzDeepCellDataNFT/EzDeepCellDataAmyloidPlaques/Abeta40Abeta42/") 

outputFolder = paste0(dataPath, "heatmapsAmyloidPlaques/")
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
fwrite(dat, paste0(outputFolder, "EZSegDeepCellHiLoPHF1Tau.csv"))

## if you want to convert the merged data to data.frame
df <- as.data.frame(dat)
rm(dat)

## to read in the csv of CellEngine 
#df = read_csv (paste0(dataPath,"EzSegAmyloidPlaques.csv")) # "EzSegAmyloidPlaques.csv" or "EzSegHiLoPHF1Tau.csv" or "EZSegDeepCellHiLoPHF1Tau.csv"

## to convert file id to gatenames
df$file = as.factor(df$file)
table(df$file)
levels(df$file) =c("Abeta42pos","Abeta40pos", "Abeta4240neg") #c("HiPHF1Tau", "LoPHF1Tau") or c("Abeta42pos","Abeta40pos", "Abeta4240neg")
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
markersGated = c("PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142") #,"PHF1Tau","ApoE4","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142"
commonCols = c(markersGated, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  = plot_data[, match(c("Amyloidbeta140","Amyloidbeta142","PanAmyloidbeta1724"),colnames(plot_data))] #plot_data[, match(c("PHF1Tau"),colnames(plot_data))]  or plot_data[, match(c("Amyloidbeta140","Amyloidbeta142","PanAmyloidbeta1724"),colnames(plot_data))]

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
plot_data  =plot_data[, match(c("PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","X8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63","EEA1"),colnames(plot_data))] #plot_data[, match(c("PHF1Tau","ApoE4","X8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63","EEA1"),colnames(plot_data))]

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
neuronalMarkers = c("MAP2","VGAT","VGLUT1","VGLUT2","Reelin","Calbindin","Calretinin","Parvalbumin","Presenilin1NTF","PanGAD6567","CD47")
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
nonNeuronalMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","CD33Lyo")
commonCols = c(nonNeuronalMarkers, "file")
df_subset = df[,which(colnames(df) %in% commonCols)]
plot_data = aggregate(df_subset, by=list(df_subset$file), FUN = mean)
plot_data$file = NULL
rownames(plot_data) =  plot_data$Group.1
plot_data$Group.1 = NULL
plot_data  = plot_data[, match(c("Iba1","CD45","CD33Lyo","GFAP","CD31","CD105","MCT1","MBP","MAG"),colnames(plot_data))]

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




libs <- c('ggplot2','RColorBrewer','reshape2','devtools',
          'pheatmap','flowCore','plyr','scales','coin','grid',
          'gridExtra','tidyr','Rtsne','FlowSOM','randomcoloR','viridis','tidyverse')

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

library(flowCore)
library(readr)
setwd("/Volumes/KausaliaHD 1/R4R")
source("/Volumes/KausaliaHD 1/R4R/General/readWriteCytofFCS.R")

libs <- c('ggplot2','ggplot','RColorBrewer','reshape2','devtools',
          'pheatmap','flowCore','plyr','scales','coin','grid',
          'gridExtra','tidyr','Rtsne','FlowSOM','randomcoloR','viridis','tidyverse')
library(umap)
library(tsne)
library(readr)
library(flowCore)
library(readr)
library(scales)

mantisPath = "/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/"
ezPath ="/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/"

outputFolder = paste0(mantisPath,"CombineGatedMantisEz/UMAP10/")

dir.create(outputFolder, showWarnings= F, recursive = T)

##NewMatlabSettings: Deepcell
#all_data_mod <- read_csv(paste0(mantisPath,"GatedPopWithMantisCuration/NewMatlabSettings/UMAP/umap_data_sampled.csv"))
all_data_mod <- read_csv(paste0(mantisPath,"GatedPopWithMantisCuration/NewMatlabSettings/UMAP/umap_data_sampledUMAPCurated.csv"))
all_data_mod$point_id = all_data_mod$Point
all_data_mod$obj_id = all_data_mod$cellLabelInImage

##OldMatlabSettings: Deepcell
#all_data_mod <- read_csv(paste0(mantisPath,"/GatedPopWithMantisCuration/OldMatlabSettings/umap_data_sampledUMAPCurated.csv"))

all_data_EzSeg = read_csv (paste0(ezPath, "EzSegData_uci2712J/all_dataEzSeg.csv"))

allMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","Reelin","MAP2","X8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau","Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1NTF","pTDP43","EEA1","CD56Lyo","CD47", "TotalTau", "PanGAD6567","CD33Lyo","HistoneH3Lyo", "Synaptophysin", "PSD95", "PanApoE2E3E4")
commonCols = c(allMarkers, "Meta","point_id","obj_id")

all_data_EzSeg$Meta = all_data_EzSeg$obj_type_id
all_data_mod$Meta = all_data_mod$MantisPopulation
all_data_mod$X8OHGuano = all_data_mod$`8OHGuano` 


#For Natural log transformation. IN USE
applyAsinh <- function(value) {
  value <- value + 1 
  value <- log(value)
  return(value)
}

#to downSample
set.seed(140214)

all_data_EzSeg_sampled = data.frame()
for (i in unique(all_data_EzSeg$Meta)){
  all_data_EzSeg_subset = subset(all_data_EzSeg, Meta == i)
  fcs.expr_transform <- all_data_EzSeg_subset[,allMarkers]
  fcs.expr_transform <- apply(fcs.expr_transform, 2, applyAsinh)
  
  # 2.1  99th Percentile narmalisation for normalisation individual marker expression
  fcsquantile <- apply(fcs.expr_transform, 2, function(x) {
    quantile(x, 0.99, na.rm = T)
  })
  
  # 2.2 If all values are zero, set quantile to 1 to avoid Undefined values
  fcsquantile[fcsquantile == 0] <- 1
  fcsquantile[is.na(fcsquantile)] <- 1
  
  # between zero and max signal/99th percentile
  fcs.expr_transform <- data.frame(t(apply(fcs.expr_transform, 1, function(x) {
    x / fcsquantile
  })))
 
  all_data_EzSeg_subset = cbind(fcs.expr_transform, all_data_EzSeg_subset[,-which(names(all_data_EzSeg_subset) %in% allMarkers)])
 
  if(i %in% c("astrocyte_process")){
    all_data_EzSeg_subset = all_data_EzSeg_subset[sample(1:nrow(all_data_EzSeg_subset),ceiling(nrow(all_data_EzSeg_subset)/10)),] # to give a sample size
    all_data_EzSeg_sampled = rbind(all_data_EzSeg_sampled, all_data_EzSeg_subset)
  } else {
  if(i %in% c("microglia_process")){
    all_data_EzSeg_subset = all_data_EzSeg_subset[sample(1:nrow(all_data_EzSeg_subset),ceiling(nrow(all_data_EzSeg_subset)/5.56)),] # to give a sample size
    all_data_EzSeg_sampled = rbind(all_data_EzSeg_sampled, all_data_EzSeg_subset)
    } else {
   
       all_data_EzSeg_sampled = rbind(all_data_EzSeg_sampled, all_data_EzSeg_subset)
       
  }
}
}

#this is to show region and point number
plot_data = rbind(all_data_mod[,commonCols], all_data_EzSeg_sampled[,commonCols])
plot_data$Meta <- as.factor(plot_data$Meta)


# To remove double counting of assign EZseg and Deepcell objects (any mglia or endothelial cells near projections (107) / vessels (108) respectively to null)
# for (cell in seq_along(plot_data)) {
#   if (plot_data[cell,]$Meta == "microglia") {
#     x = plot_data[cell,]$x_centroid
#     y = plot_data[cell,]$y_centroid
#     #filter projections by -10 & + 10 in x and y direction (3 for mglia-processes, 4 for vessels)
#     data_indexes = which(with(plot_data, (Meta == "microglia_process" & (x_centroid < abs(x-MajorAxisLength/2) | x_centroid < abs(x+MajorAxisLength/2) & (y_centroid < abs(x-MajorAxisLength/2) | y_centroid < abs(x+MajorAxisLength/2))))))
#     plot_data[data_indexes, "Meta"] = "dumpbucket" #
#   }
# }
# filter(plot_data, Meta == 'cells') %>% select(Meta) %>% table() # check plot_data$Meta
# plot_data %>% select(Meta) %>% table()

umap_data = plot_data
#umap_data = subset(umap_data, Meta %in% c("neurons", "endothelial","microglia","astrocytes"))
#umap_data$Meta = ordered(umap_data$Meta, levels = c("microglia","astrocytes","endothelial","neurons")) # in console to check if ordered : table(umap_data_sampled$Meta)

umap_data_sampled = umap_data

##for plotting Specific Markers
#umapLineageMarkers= c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","MAP2","Calbindin","Calretinin","Parvalbumin","Presenilin1") # Phenotypic Markers
#umapADSpecificMarkers = c("8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","MFN2","PHF1Tau","PolyubiK48","PolyubiK63","pTDP43")
umapGCMarkers = c("MBP","MAG","MAP2","Iba1","CD45","GFAP","CD31","MCT1","CD105") # Gated markers

umap_data = umap_data[,umapGCMarkers]
#umap_data = umap_data[,umapAllMarkers]
#umap_data = umap_data[,umapLineageMarkers]
#umap_data = umap_data[,umapADSpecificMarkers]

# to find order of populations
table(umap_data_sampled$Meta) 

#Custom umap settings
customSettings = umap.defaults # Run umap.defaults to see original settings
customSettings$n_neighbors = 30 # knn
customSettings$min_dist = 0.01 # distance betweeen 2 closet cells.  default 0.1, at 0.05 farway, at 0.5 close together
customSettings$local_connectivity = 1 # 
customSettings$negative_sample_rate = 15 # choose 5 random points 
customSettings$knn_repeats = 3 # 3 repeats for each iteration or epoch of 200 attempts
# customSettings$a = 0.6 # higher numbers brings similar clusters closer
# customSettings$b = 0.7 # higher numbers brings different clusters closer

#Scale data
umap_data = apply(umap_data, 2, scale)#scale by col=2, by row=1
umap_data = rescale(umap_data, to = c(0, 1))

#for t-sne plot
# tsne_out= tsne(umap_data)
# tsnecoord = tsne_out$Y
# colnames(tsnecoord) = c("tsne1", "tsne2")
# umap_data_sampled$umapX = tsnecoord[,"tsne1"]
# umap_data_sampled$umapY = tsnecoord[,"tsne2"]
#for t-sne plot

seed = 211264
set.seed(seed)
#for umap plot
umap_out = umap(umap_data, config = customSettings,  verbose = T) 
umap_data_sampled$umapX = umap_out$layout[,1] 
umap_data_sampled$umapY = umap_out$layout[,2]

## To Plot Umaps

## Umap of old
# p = ggplot(umap_data_sampled, aes(umapX, umapY, color = as.factor(Meta))) + 
#   geom_point(pch = 19, cex =2.5 , alpha=0.8) +
#   theme_bw() + labs(x=NULL, y=NULL) +
#   scale_colour_manual(values = palette) +
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         panel.border=element_blank(),
#         axis.ticks=element_blank(),
#         axis.text=element_blank()) + 
#   guides(color = guide_legend(override.aes = list(size=5)))
# ggsave(paste0(outputFolderUmap, "/HiADCluster_umapOriginal.pdf"), plot = p, width = 12, height = 10)

## Umap of Mantis Population

# To get coordinates, run from ggplot .... up to but not including theme(panel.....)
p = ggplot(umap_data_sampled, aes(umapX, umapY, color = as.factor(Meta))) + 
  geom_point(pch = 19, cex =0.01 , alpha=0.8) +
  theme_bw() + labs(x=NULL, y=NULL) +
  guides(color = guide_legend(override.aes = list(size=length(palette))))
ggsave(paste0(outputFolder, "/HiADCluster_CombineMantisEz_UmapScaling.pdf"), plot = p, width = 7, height = 5)

## To Plot All Populations (coloured by unique Ez and GatedMantis population)
plot_data = umap_data_sampled
table(plot_data$Meta) #to check number of cell types per Lineage

## amyloidopathy Fuchia "#FF007F", astrocytes_process LightGreen "#99FF33", astrocytes Deepcell DarkGreen "#106900",  endothelial Deepcell DarkBlue "#0028D4",  microglia Deepcell Mango "#FFBF00", microglia_process Yellow "#FFFF00",neurons DarkRed "#980A0A",tangle-threads "#6335A1", vessel_CD31_CD105 Baby Blue "#2D9BCA",vessel_MCT1 royal Blue "#0000FF"
uniquePalette = c("#FF007F","#99FF33","#106900","#0028D4","#FFBF00","#FFFF00","#980A0A","#6335A1","#2D9BCA","#0000FF")

table(umap_data_sampled$Meta)
palette = c("#FFBF00","#106900","#0028D4","#980A0A")

p = ggplot(umap_data_sampled, aes(umapX, umapY, color = as.factor(Meta))) + 
  geom_point(pch = 19, cex =0.5, alpha=0.8) +
  theme_bw() + labs(x=NULL, y=NULL) +
  scale_colour_manual(values = uniquePalette) +
  xlim(-16,16) +
  ylim(-15,15) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank()) + 
  guides(color = guide_legend(override.aes = list(size=length(uniquePalette))))
ggsave(paste0(outputFolder, "/CombineMantisEzUnique.pdf"), plot = p, width = 10, height = 10)

## To Plot All Populations (coloured by common Ez and GatedMantis population)
plot_data = umap_data_sampled
table(plot_data$Meta) #to check number of cell types per Lineage

## amyloidopathy Fuchia "#FF007F", astrocytes DarkGreen "#106900",astrocytes_process DarkGreen "#106900", endothelial DarkBlue "#0028D4",  microglia Mango "#FFBF00", microglia_process Mango "#FFBF00",neurons DarkRed "#980A0A" ,tangle-threads "#6335A1", vessel_CD31_CD105 DarkBlue "#0028D4",vessel_MCT1 DarkBlue "#0028D4"
#commonPalette = c("#FF007F","#106900","#106900","#0028D4","#FFBF00","#FFBF00","#980A0A","#6335A1","#0028D4","#0028D4") # option 1
commonPalette = c("#FF007F","#106900","#106900","#0028D4","#FFBF00","#FFBF00","#980A0A","#00ffff","#0028D4","#0028D4") # option 2

p = ggplot(umap_data_sampled, aes(umapX, umapY, color = as.factor(Meta))) + 
  geom_density2d(color = "#3c3b3f", alpha = 0.5) +
  geom_point(pch = 19, cex =0.5, alpha=0.8) +
  theme_bw() + labs(x=NULL, y=NULL) +
  scale_colour_manual(values = commonPalette) +
  xlim(-16,16) +
  ylim(-15,15) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank()) + 
  guides(color = guide_legend(override.aes = list(size=length(commonPalette))))
ggsave(paste0(outputFolder, "/CombineMantisEzCommon.pdf"), plot = p, width = 10, height = 10)

#To Plot Subsets

#Lineage Markers
plot_data_subsetLineage = subset(plot_data, Meta %in% c("astrocyte_process", "astrocytes","endothelial","microglia","microglia_process","neurons","vessel_CD31_CD105","vessel_MCT1"))
table(plot_data_subsetLineage$Meta) #to check number of cell types per Lineage

## astrocytes DarkGreen "#106900",astrocytes_process DarkGreen "#106900", endothelial DarkBlue "#0028D4",  microglia Mango "#FFBF00", microglia_process Mango "#FFBF00",neurons DarkRed "#980A0A", vessel_CD31_CD105 DarkBlue "#0028D4",vessel_MCT1 DarkBlue "#0028D4"
lineagePalette = c("#106900","#106900","#0028D4","#FFBF00","#FFBF00","#980A0A","#0028D4","#0028D4")

t=ggplot(plot_data_subsetLineage, aes(umapX, umapY, color = Meta)) + 
  geom_density2d(data = plot_data, color = "#3c3b3f", alpha = 0.5) +
  geom_point(pch = 19, cex = 0.5, alpha=0.8) +
  xlim(-16,16) +
  ylim(-15,15) +
  scale_colour_manual(values = lineagePalette) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme_bw() + labs(x=NULL, y=NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())
ggsave(paste0(outputFolder, "/CombineMantisEzLineage",".pdf"), plot = t, width = 10, height = 10)

## DeepCell Objects
plot_data_subsetDeepCell = subset(plot_data, Meta %in% c("astrocytes","endothelial","microglia","neurons"))
table(plot_data_subsetDeepCell$Meta) #to check number of cell types per Lineage

## astrocytes DarkGreen "#106900",astrocytes_process DarkGreen "#106900", endothelial DarkBlue "#0028D4",  microglia Mango "#FFBF00", microglia_process Mango "#FFBF00",neurons DarkRed "#980A0A", vessel_CD31_CD105 DarkBlue "#0028D4",vessel_MCT1 DarkBlue "#0028D4"
#deepCellPalette = c("#106900","#0028D4","#FFBF00","#980A0A")
#deepCellPalette = c("#000000","#000000","#000000","#000000")
deepCellPalette = c("#603813","#603813","#603813","#603813")

t=ggplot(plot_data_subsetDeepCell, aes(umapX, umapY, color = Meta)) + 
  geom_density2d(data = plot_data, color = "#3c3b3f", alpha = 0.5) +
  geom_point(pch = 19, cex = 0.2, alpha = 0.3) +
  xlim(-16,16) +
  ylim(-15,15) +
  scale_colour_manual(values = deepCellPalette) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme_bw() + labs(x=NULL, y=NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())
ggsave(paste0(outputFolder, "/CombineMantisDeepCellObjects",".pdf"), plot = t, width = 10, height = 10)

## EzSeg Lineage Objects
plot_data_subsetEzSeg = subset(plot_data, Meta %in% c("astrocyte_process","microglia_process","vessel_CD31_CD105","vessel_MCT1"))
table(plot_data_subsetEzSeg$Meta) #to check number of cell types per Lineage

##  astrocytes_process LightGreen "#99FF33",  microglia_process Yellow "#FFFF00", vessel_CD31_CD105 Baby Blue "#2D9BCA",vessel_MCT1 royal Blue "#0000FF"
#ezSegPalette = c("#99FF33","#FFFF00","#2D9BCA","#0000FF")
#ezSegPalette = c("#000000","#000000","#000000","#000000")
ezSegPalette = c("#605c3c","#605c3c","#605c3c","#605c3c")
  
t=ggplot(plot_data_subsetEzSeg, aes(umapX, umapY, color = Meta)) + 
  geom_density2d(data = plot_data, color = "#3c3b3f", alpha = 0.5) +
  geom_point(pch = 19, cex = 0.2, alpha = 0.3) +
  xlim(-16,16) +
  ylim(-15,15) +
  scale_colour_manual(values = ezSegPalette) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme_bw() + labs(x=NULL, y=NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())
ggsave(paste0(outputFolder, "/CombineMantisEzSegObjects",".pdf"), plot = t, width = 10, height = 10)

## Disease Structures (I)
plot_data_DiseaseI = umap_data_sampled
table(plot_data_DiseaseI$Meta) #to check number of cell types per Lineage

## amyloidopathy Fuchia "#FF007F", astrocytes gray #A0A0A0", astrocytes_process gray #A0A0A0", endothelial gray #A0A0A0"",microgliagray #A0A0A0", microglia_process gray #A0A0A0",neurons gray #A0A0A0",tanglesthreads depppurple "#6335A1", vessel_CD31_CD105 gray #A0A0A0",vessel_MCT1 gray #A0A0A0"
diseasePaletteI = c("#FF007F","#E0E0E0","#E0E0E0","#E0E0E0","#E0E0E0","#E0E0E0","#E0E0E0","#00ffff","#E0E0E0","#E0E0E0")

p = ggplot(plot_data_DiseaseI, aes(umapX, umapY, color = as.factor(Meta))) + 
  geom_density2d(data = plot_data, color = "#3c3b3f", alpha = 0.5) +
  geom_point(pch = 19, cex =0.5, alpha= 1) +
  theme_bw() + labs(x=NULL, y=NULL) +
  scale_colour_manual(values = diseasePaletteI) +
  xlim(-16,16) +
  ylim(-15,15) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank()) + 
  guides(color = guide_legend(override.aes = list(size=length(diseasePaletteI))))
ggsave(paste0(outputFolder, "/CombineMantisEzDiseaseI.pdf"), plot = p, width = 10, height = 10)

## Disease Structures (II)
plot_data_subsetDiseaseII = subset(plot_data, Meta %in% c("amyloidPlaques","tauTangles"))
table(plot_data_subsetDiseaseII$Meta) #to check number of cell types per Lineage

## amyloidopathy Fuchia "#FF007F",tanglesthreads cyan "#00ffff"
diseasePaletteII = c("#FF007F","#00ffff")

t=ggplot(plot_data_subsetDiseaseII, aes(umapX, umapY, color = Meta)) +
  geom_density2d(data = plot_data, color = "#3c3b3f", alpha = 0.5) +
  geom_point(pch = 19, cex = 0.5, alpha = 2) +
  theme_bw() + labs(x=NULL, y=NULL) +
  scale_colour_manual(values = diseasePaletteII) +
  xlim(-16,16) +
  ylim(-15,15) +
  guides(color = guide_legend(override.aes = list(size=10))) +
  theme_bw() + labs(x=NULL, y=NULL) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())
ggsave(paste0(outputFolder, "/CombineMantisEzDiseaseII",".pdf"), plot = t, width = 10, height = 10)

##Umap by channel
plot_data_out = umap_data_sampled[, c(allMarkers)]
plot_data_out = apply(plot_data_out, 2, function(x){
  x/max(x)
})
plot_data_out = cbind(plot_data_out, umap_data_sampled[, c("Meta","umapX", "umapY")])
##Settings for NonNeuronNeuron
for (channel in allMarkers) {
  maxlimit = 1
  if (channel %in%  c("PanApoE2E3E4"))
    maxlimit = 0.0002
  
  if (channel %in%  c("EEA1"))
    maxlimit = 0.001
  
  if (channel %in%  c("ApoE4","PolyubiK63"))
    maxlimit = 0.003
  
  if (channel %in%  c("CD31","Amyloidbeta140","MFN2"))
    maxlimit = 0.005
  
  if (channel %in%  c("PHF1Tau"))
    maxlimit = 0.008
  
  if (channel %in%  c("PanAmyloidbeta1724","PolyubiK48"
                      ,"MCT1","Calretinin"))
    maxlimit = 0.01
  
  if (channel %in%  c("Presenilin1NTF"))
    maxlimit = 0.015
  
  if (channel %in%  c())
    maxlimit = 0.15
  
  if (channel %in%  c("pTDP43","CD105","Parvalbumin"))
    maxlimit = 0.03
  
  if (channel %in%  c("MAP2"))
    maxlimit = 0.06
  
  if (channel %in%  c("Reelin","PanGAD6567","Synaptophysin","X8OHGuano","CD33Lyo","Amyloidbeta142","CD47","VGLUT2"))
    maxlimit = 0.05
  
  if (channel %in%  c("TotalTau"))
    maxlimit = 0.07
  
  if (channel %in%  c("CD45"))
    maxlimit = 0.1
  
  if (channel %in%  c("CD56Lyo"))
    maxlimit = 0.5
  
  if (channel %in%  c("GFAP"))
    maxlimit = 0.19
  
  if (channel %in%  c("VGAT","Calbindin","Iba1","VGLUT1"))
    maxlimit = 0.05
  
  if (channel %in%  c("PSD95"))
    maxlimit = 0.15
  
  if (channel %in%  c())
    maxlimit = 0.25
  
  if (channel %in%  c("MBP","MAG","TH","HistoneH3Lyo"))
    maxlimit = 0.3
  
  if (channel %in%  c("TotalTau"))
    maxlimit = 0.2
  
  plot_data_out[which(plot_data_out[, channel] >maxlimit), channel] = maxlimit
  
  if (channel == "8OHGuano")
    channel = "`8OHGuano`"
  breakList = c(0,0,0,0.1,0.1,0.1,0.1,0.2,0.3,1)  # adds colour break in every value add more 0 to compress the blues. VeryStrong c(0,0.1,0.2,0.3,1), Strong c(0,0.1,0.1,0.1,0.1,0.2,0.3,1), Weak c(0,0,0,0.1,0.1,0.1,0.1,0.2,0.3,1), VeryWeak c(0,0,0,0,0,0.1,0.2,0.3,1)
  p = ggplot(plot_data_out, aes_string("umapX", "umapY", color = channel)) + 
    #geom_density2d(color = "#808080", alpha = 0.2) +
    geom_point(pch = 10, cex = 0.3, alpha=0.5) +
    xlim(-16,16) +
    ylim(-15,15) +
    theme_bw() + labs(x=NULL, y=NULL) +
    #scale_color_manual(values  =  rev(colorRampPalette(brewer.pal(11,"RdBu")) (100))) +
   # scale_color_manual(values = colorRampPalette(c("#202020","#00CCCC"))(100)) +
   # scale_color_brewer(palette = "Greens") +
    #scale_color_viridis_c(direction = 1, option = "D", values = breakList) +
     scale_color_gradient2(
     low = ("#26094d"),
     mid = ("#752c40"),
     high = ("#cb5334"),
     midpoint = maxlimit/2,
     space = "Lab",
     na.value = "gray50",
     guide = "colourbar",
     aesthetics = c("fill","colour")) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank())
  ggsave(paste0(outputFolder, "/CombineMantisEz_", channel, ".pdf"), plot = p, width = 15, height = 15)
}
## To find umapY and UmapY coordinates to center plot and remover outliers
range(plot_data_out$umapY)
range(plot_data_out$umapX)

##Generate csv 
write.csv(umap_data_sampled, paste0(outputFolder, "/umap_data_sampledCombineMantisEz1.csv"))

##Generate fsc for both Ezseg and Deepcelldata
plot_data = umap_data_sampled
plot_data$Region = factor(plot_data$Region)
levels(plot_data$Region) = c('DG','CA2')
plot_data$brainRegion = paste0('p', plot_data$Point, '_',plot_data$Region)

# plot_data$Region = as.numeric(plot_data$Region)
# plot_data$SampleNames = as.numeric(plot_data$SampleNames)
# plot_data$Point = as.numeric(plot_data$Point)
# plot_data$MantisPopulation = as.numeric(as.factor(plot_data$MantisPopulation))
plot_data$Meta = as.numeric(as.factor(plot_data$Meta))
plot_data$brainRegion = as.numeric(as.factor(plot_data$brainRegion))


#to remove NA columns
plot_data = plot_data[ , colSums(is.na(plot_data)) == 0]

fcsExprs <- flowCore::flowFrame(as.matrix(plot_data))
suppressWarnings(flowCore::write.FCS(fcsExprs, paste0(outputFolder, 'alldataMantisCurated.fcs'), what = "numeric"))

##Generate fsc for Ezseg Deepcelldata
plot_data = all_data_EzSeg

plot_data$obj_type_id = as.numeric(as.factor(plot_data$obj_type_id))
plot_data$run_type_id = as.numeric(as.factor(plot_data$run_type_id))
fcsExprs <- flowCore::flowFrame(as.matrix(plot_data))
suppressWarnings(flowCore::write.FCS(fcsExprs, paste0(outputFolder, 'alldataEzSeg.fcs'), what = "numeric"))

#______________________________________________________________________________________________________________________________________________________________

#barplot for number of cell per population
plot_data = umap_data_sampled
table(plot_data$Meta) #to check number of cell types per Lineage
plot_data$Meta = ordered(plot_data$Meta, levels = c("amyloidPlaques","tauTangles","neurons","microglia","microglia_process","astrocyte_process","astrocytes","endothelial","vessel_CD31_CD105","vessel_MCT1"))
population = c("amyloidPlaques","tauTangles")
plot_data_subset = droplevels(subset(plot_data, Meta %in% population))
table(plot_data_subset$Meta)

##to remove 0 counts population
droplevels(plot_data_subset)

## amyloidopathy Fuchia "#FF007F", astrocytes_process LightGreen "#99FF33", astrocytes Deepcell DarkGreen "#106900",  endothelial Deepcell DarkBlue "#0028D4",  microglia Deepcell Mango "#FFBF00", microglia_process Yellow "#FFFF00",neurons DarkRed "#980A0A",tauopathy Cyan "#00CCCC", vessel_CD31_CD105 Baby Blue "#2D9BCA",vessel_MCT1 royal Blue "#0000FF"
uniquePalette = c("#FF007F","#99FF33") #,"#106900","#0028D4","#FFBF00","#FFFF00","#980A0A","#00CCCC","#2D9BCA","#0000FF")

numberOfCellperPopulation = table(plot_data_subset$Meta)
write.csv(numberOfCellperPopulation, paste0(outputFolder,'/numberOfCellPerPopulation.csv'))
numberOfCellperPopulation = data.frame(Meta=numberOfCellperPopulation)
colnames(numberOfCellperPopulation) = c('Meta', 'cellCounts')
h=ggplot(numberOfCellperPopulation, aes(x=Meta, y=cellCounts, fill = Meta)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = uniquePalette) +
  theme_classic(15)
ggsave(paste0(outputFolder, "/CellNumPerPopulationD",".pdf"), plot = h, width = 12, height = 10)









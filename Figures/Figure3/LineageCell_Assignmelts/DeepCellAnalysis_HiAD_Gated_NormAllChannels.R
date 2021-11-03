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
setwd("/Volumes/KausaliaHD/R4R")
source("/Volumes/KausaliaHD/R4R/General/readWriteCytofFCS.R")

#For Natural log transformation. IN USE
applyAsinh <- function(value) {
  value <- value + 1 
  value <- log(value)
  return(value)
}

#Make vector of Markers Names across all the input files
allMarkerNames = c ()
concatFCStoDF <- function(file.names, region) {
  cat("Reading FCS files...")
  cell.data <- data.frame()
  files.FCS <- read.flowSet(file.names, transformation = FALSE, emptyValue = FALSE, truncate_max_range = FALSE)
  cat("Complete\n")
  
  #colnames(files.FCS) <- channel.metals
  colnames(files.FCS) <- gsub("-|_|\\s+", "", colnames(files.FCS), ignore.case = T)
  print(colnames(files.FCS))
  
  allMarkerNames <<-c(allMarkerNames, colnames(files.FCS))
}

##Paper uci2717J DG AND CA2 Scans
outputFolder = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/GatedPoP/NewMatlabSettings/4Gated0FlowSOMClustersNoOligo1/"
dir.create(outputFolder, showWarnings= F, recursive = T)

fcsPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSegData_uci2712J/"
fcsRuns = c("190505HiResDG/", "190604HiResCA2/") 
fcsPixel = "fcs_single_cell_dynamic_expansion_Smaller_no_fftnoiseDist"


allDataFcs = data.frame()
allFcsFileNames =c()

for (runNumber in 1:length(fcsRuns)) {
  print(paste0(fcsPath,fcsRuns[runNumber]))
  #TMAPanel <- read_csv(list.files(path=paste0(fcsPath, fcsRuns[runNumber],"/info/"),pattern ="*.csv", full.names = T))
  fcsNames = list.files(path = paste0(fcsPath,fcsRuns[runNumber], fcsPixel), pattern = paste0("^dataScaleSizeFCS.*fcs$"), full.names = T)
  allFcsFileNames = c(allFcsFileNames, fcsNames)
  
  #use this if relying info.csvpanel for channel names
  #allDataFcs = rbind(allDataFcs, concatFCStoDF(fcsNames, fcsRuns[runNumber], TMAPanel$Label))
  
  #use this if relying on tifs folder for channel names
  #channelnames = list.files(path = paste0(fcsPath,fcsRuns[runNumber], '/no_fftnoise/Point1/TIFs/'), pattern = paste0("*.tif"), full.names = F)
  #print(channelnames)
  allDataFcs = rbind(allDataFcs, concatFCStoDF(fcsNames, fcsRuns[runNumber]))
}

allMarkerNames = unique(allMarkerNames)
concatFCStoDF <- function(file.names, region) {
  cat("Reading FCS files...")
  cell.data <- data.frame()
  files.FCS <- read.flowSet(file.names, transformation = FALSE, emptyValue = FALSE, truncate_max_range = FALSE)
  cat("Complete\n")
  
  colnames(files.FCS) <- gsub("-|_|\\s+", "", colnames(files.FCS), ignore.case = T)
  
  missingMarkerNames = setdiff(allMarkerNames,colnames(files.FCS))
  
  files.FCS_transform <- files.FCS[, -which(colnames(files.FCS) %in% c("cellLabelInImage","cellSize"))]
  print(colnames(files.FCS_transform))
  
  files.FCS_nontransform <- files.FCS[, c("cellLabelInImage","cellSize")]
  #print (files.FCS_nontransform)
  # files.FCS_nontransform$HistoneH3 <- files.FCS_transform$HistoneH3
  # files.FCS_transform$HistoneH3 <- NULL
  #print(colnames(files.FCS_transform))
  #files.FCS_transform <- files.FCS[, colnames(files.FCS) %in% channel.metals]
  #files.FCS_nontransform <- files.FCS[, !(colnames(files.FCS) %in% channel.metals)]
  
  for (i in 1:length(files.FCS)) {
    cat(paste("Concatenating:", file.names[i], "\n", sep = " "))
    
    file.fcs_transform <- files.FCS_transform[[i]]
    fcs.expr_transform <- as.data.frame(exprs(file.fcs_transform))
    fcs.expr_transform <- apply(fcs.expr_transform, 2, applyAsinh)
    
    # # Percentile narmalisation across runs (100th percentile)
    # fcsquantile <- apply(fcs.expr_transform, 2, function(x) {
    #   quantile(x, 1)
    # })
    # 
    # # If all values are zero, set quantile to 1 to avoid Undefined values
    # fcsquantile[fcsquantile == 0] <- 1
    # 
    # # Rescale the data to [0:5] interval to mimic CyTOF data scale
    # fcs.expr_transform <- data.frame(t(apply(fcs.expr_transform, 1, function(x) {
    #   x / fcsquantile * 5
    # })))
    
    # Normalise signal intensites before summing by:
    
    # 1. Unit vector calculation for narmalisation across different runs (Batch effect). need to change the way data is imported. wont work in this format
    # fcs.expr_transform <- data.frame(t(apply(fcs.expr_transform, 1, function(x) {
    #   x / sqrt(sum(x^2))
    # })))
    
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
    
    # 3.0 Quantile for narmalisation across different runs (Batch effect)
    #to fill in use nonsignal channels like C12 to calculate quantile normalisation
    
    colnames(fcs.expr_transform) = gsub("X8OHGuano", "8OHGuano", colnames(fcs.expr_transform))
    
    file.fcs_nontransform <- files.FCS_nontransform[[i]]
    fcs.expr_nontransform <- as.data.frame(exprs(file.fcs_nontransform))
    
    fcs.expr <- cbind(fcs.expr_transform, fcs.expr_nontransform)
    
    fcs.expr[missingMarkerNames] = NA
    
    fcs.expr$Region = region
    fcs.expr$Point = gsub('\\.fcs', '', gsub('dataScaleSizeFCS_p', '', basename(file.names[i])))
    
    cell.data <- rbind(cell.data, fcs.expr)
  }
  cat("Concatenation complete.\n")
  
  return(cell.data)
}

#To import data
allDataFcs = data.frame()
allFcsFileNames =c()

for (runNumber in 1:length(fcsRuns)) {
  print(paste0(fcsPath,fcsRuns[runNumber]))
  fcsNames = list.files(path = paste0(fcsPath,fcsRuns[runNumber], fcsPixel), pattern = paste0("^dataScaleSizeFCS.*fcs$"), full.names = T)
  allFcsFileNames = c(allFcsFileNames, fcsNames)
  
  #use this if relying info.csvpanel for channel names
  #allDataFcs = rbind(allDataFcs, concatFCStoDF(fcsNames, fcsRuns[runNumber], TMAPanel$Label))
  
  #use this if relying on tifs folder for channel names
  #channelnames = list.files(path = paste0(fcsPath,fcsRuns[runNumber], '/no_fftnoise/Point1/TIFs/'), pattern = paste0("*.tif"), full.names = F)
  #print(channelnames)
  allDataFcs = rbind(allDataFcs, concatFCStoDF(fcsNames, fcsRuns[runNumber]))
}

markers = setdiff(allMarkerNames, c("cellLabelInImage" ,"cellSize","Au197","Background", "C12" , "Ca40", "Si28" ,"Ta181", "Na23" ,"empty113" ,"empty139"))

plot_histograms = function(populationName ) {
  #plot histograms of all channels before outlier removal
  plot_data = reshape2::melt(allDataFcs[,markers], ids = c())
  colnames(plot_data) = c("Markers","LogCellIntensities")
  ggplot(plot_data, aes(x=LogCellIntensities)) + 
    geom_histogram(bins = 30)+ 
    facet_wrap(Markers~., scales="free")
  ggsave(filename = paste0(outputFolder,"/MarkersHistograms_",populationName,".pdf"), width = 20, height = 20)
}

allDataFcs =  allDataFcs[!is.na(allDataFcs$CD45),]

#to save the full ungated data set.
allDataFcsFull = allDataFcs 

#____________________________________________________________________________________________________________________________________________________________________
# Make a heatmap for channels vs all cell with expression values from allDATAFcsFull
# plot_data = allDataFcsFull # for all ungated cells
# plot_data = subset(all_data_full, Meta == 0) # for Meta gate 1 and if set to 0 see remain in ungated cells
# plot_data = plot_data[, c(markers)]
# 
# # to fing out number of cell in each condition
# sum(!is.na(plot_data$Amyloidbeta140)) # number of AD cells
# sum(!is.na(plot_data$pASyn)) # number of PD cells
# 
# # plot_data = apply(plot_data, 2, scale)#marker scaled
# # plot_data = rescale(plot_data, to = c(0, 1))#adjust between 0 and 2
# 
# # MarkersGating = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","TH","VGAT","VGLUT1","VGLUT2", "MAP2")
# # plot_data_subset = plot_data[,which(colnames(plot_data) %in% MarkersGating)]
# # dev.off()
# # pdf(paste0(outputFolder, "/HeatMap_GCMarkers.pdf"), width = 6, height = 6, onefile = F)
# # 
# # MarkersCommon = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","TH","VGAT","VGLUT1","VGLUT2","MAP2","CD33Lyo","TotalTau","CD56Lyo","SERT", "Synaptophysin","PSD95")
# # plot_data_subset = plot_data[,which(colnames(plot_data) %in% MarkersCommon)]
# # dev.off()
# # pdf(paste0(outputFolder, "/HeatMap_CommonMarkers.pdf"), width = 6, height = 6, onefile = F)
# 
# MarkersADSpecific = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","MAP2","8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau","PanAmyloidbeta1724","Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1","Reelin","pTDP43","CD33Lyo")
# plot_data_subset = plot_data[,which(colnames(plot_data) %in% MarkersADSpecific)]
# dev.off()
# pdf(paste0(outputFolder, "/HeatMap_MarkersADSpecificMeta0.pdf"), width = 6, height = 6, onefile = F)
# 
# # MarkersPDSpecific = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","TH","VGAT","VGLUT1","VGLUT2","MAP2","8OHGuano","PanAmyloidbeta1724","Amyloidbeta142","ApoE4","CD47","DAT","EEA1","LRRK2Parkin8","PARK7DJ1","PHF1Tau","PINK1PARK6","PanASyn","PanGAD6567","ParkinPARK2","pASyn","pTDP43","MFN2")
# # plot_data_subset = plot_data[,which(colnames(plot_data) %in% MarkersPDSpecific)]
# # dev.off()
# # pdf(paste0(outputFolder, "/HeatMap_MarkersPDSpecific.pdf"), width = 6, height = 6, onefile = F)
# 
# breakList = c(0.0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
# colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
# 
# pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "All Cells Marker Expression Scaled",
#          show_rownames = F,
#          #  color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
#          col= colfunc(length(breakList)),breaks = breakList, cellwidth = 10, cellheight = 0.2)
# #col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)
# dev.off()

#__________________________________________________________________________________________________________
#plot Scatterplot. change names
# ggplot(allDataFcsFull, aes(x=CD31, y=MAP2)) + #__`8OHGuano` only for this marker_______________________________________________________________________________________
#   #stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", bins = 50, alpha = 0.3) +
#   #scale_fill_viridis_c(direction = -1, option = "E") +
#   geom_point(aes(colour = MCT1), size=1.5, alpha = 0.5) +
#   ylim(0,2) +
#   xlim(0,2)+
#   #scale_x_continuous(trans = "log10") +
#   #scale_y_continuous(trans = "log10") +
#   theme_classic()
# ggsave(paste0(outputFolder, "/ScatterPlotCD31vsCD105.pdf") , height = 6, width = 6) #______________________________________________________________________
#__________________________________________________________________________________________________________

#START OF GATING
all_data_full = allDataFcsFull[, -which(colnames(allDataFcsFull) %in% c("Au197","Background", "C12" , "Ca40", "Si28" ,"Ta181", "Na23" ,"empty113" ,"empty139"))] # all gated data set together

all_data_full$Meta = 0

gatedIndexes = c()
gateNames = c("neurons","endothelial","microglia","astrocytes")
#gateNames = c("neurons","endothelial","microglia","astrocytes","oligodendrocytes")
#gateNames = c("ungated","microglia","astrocytes","endothelial","excitatory","oligodendrocytes","inhibitory","dopaminergic")

## To Gate Population.
#data_indexes = which(with(allDataFcsFull, CD31 > 0.1 & MCT1 > 0.1 & CD105 > 0.1)) # endothelial Strictergates
data_indexes = which(with(allDataFcsFull, CD31 > 0.07 & MCT1 > 0.07 & CD105 > 0.07)) # endothelial
data_indexes = setdiff(data_indexes, gatedIndexes)
all_data_full[data_indexes, "Meta"] = 1 # to assign cluster ID for endothelial
plot_histograms(gateNames[2])
gatedIndexes = c(gatedIndexes, data_indexes)

#data_indexes = which(with(allDataFcsFull, (CD45 > 0.3 & Iba1 > 0.5))) # microglia Strictergates
data_indexes = which(with(allDataFcsFull, (CD45 > 0.08 & Iba1 > 0.06) | (Iba1 > 0.4) )) # microglia 
data_indexes = setdiff(data_indexes, gatedIndexes)
all_data_full[data_indexes, "Meta"] = 2 # to assign cluster ID for microglia good
plot_histograms(gateNames[3])
gatedIndexes = c(gatedIndexes, data_indexes)

#data_indexes = which(with(allDataFcsFull, GFAP > 0.6))  # astrocytes Strictergates
data_indexes = which(with(allDataFcsFull, GFAP > 0.15)) # astrocytes
data_indexes = setdiff(data_indexes, gatedIndexes)
all_data_full[data_indexes, "Meta"] = 3 # to assign cluster ID for astrocytes good
plot_histograms(gateNames[4])
gatedIndexes = c(gatedIndexes, data_indexes)

# data_indexes = which(with(allDataFcsFull, MAP2 > 0.75)) # neurons
# data_indexes = setdiff(data_indexes, gatedIndexes)
# all_data_full[data_indexes, "Meta"] = 4 # to assign cluster ID for VGLUT2 neurons
#plot_histograms(gateNames[5])
#gatedIndexes = c(gatedIndexes, data_indexes)


# data_indexes = which(with(allDataFcsFull, MAP2 > 0.001 & VGLUT1 > 0.001 & VGLUT2 > 0.002 & VGAT > 0.005 & TH > 0.04)) # neurons
# data_indexes = setdiff(data_indexes, gatedIndexes)
# all_data_full[data_indexes, "Meta"] = 4 # to assign cluster ID for VGLUT2 neurons
# plot_histograms(gateNames[5])
# gatedIndexes = c(gatedIndexes, data_indexes)

# data_indexes = which(with(allDataFcsFull,VGLUT1 > 0.3 & VGLUT2 > 0.2)) # excitatory
# data_indexes = setdiff(data_indexes, gatedIndexes)
# all_data_full[data_indexes, "Meta"] = 4 # to assign cluster ID for VGLUT2 neurons
# plot_histograms(gateNames[5])
# gatedIndexes = c(gatedIndexes, data_indexes)

# data_indexes = which(with(allDataFcsFull, MBP > 0.42 & MAG > 0.42)) # oligodendrocytes
# data_indexes = setdiff(data_indexes, gatedIndexes)
# all_data_full[data_indexes, "Meta"] = 4 # to assign cluster ID for oligodendrocytes
# plot_histograms(gateNames[5])
# gatedIndexes = c(gatedIndexes, data_indexes)

# data_indexes = which(with(allDataFcsFull, VGAT > 0.3)) # VGAT
# data_indexes = setdiff(data_indexes, gatedIndexes)
# all_data_full[data_indexes, "Meta"] = 6 # to assign cluster ID for VGAT neurons
# plot_histograms(gateNames[7])
# gatedIndexes = c(gatedIndexes, data_indexes)
# 
# data_indexes = which(with(allDataFcsFull, TH > 0.3)) # dopaminergic
# data_indexes = setdiff(data_indexes, gatedIndexes)
# all_data_full[data_indexes, "Meta"] = 7 # to assign cluster ID for dopaminergic neurons
# plot_histograms(gateNames[8])
# gatedIndexes = c(gatedIndexes, data_indexes)

table(all_data_full$Meta)

## left with ungated
allDataFcs = subset(all_data_full, Meta == 0) # ungated
allDataFcs$Meta = NULL
plot_histograms(gateNames[0])
## left with ungated so run flow SOM

# numClusters = 1 # dont change cos original flowSom cluster
numGates = 4

#_____________________________________________________________________________________________________________________________________________________________________________________
all_data = all_data_full
all_data$Meta = as.factor (all_data$Meta)
levels(all_data$Meta) = gateNames
#numClusters = numClusters + numGates

sampleNames = c("190505HiResDG","190604HiResCA2")   

all_data$Region = factor(all_data$Region)
levels(all_data$Region) = c('DG','CA2')
all_data$brainRegion = paste0('p', all_data$Point, '_',all_data$Region)
all_data$SampleNames = as.factor(all_data$brainRegion)
# levels(all_data$SampleNames) = sampleNames   

table(all_data[,c("Meta","SampleNames")])

#to save into CSV file
write.csv(all_data,  paste0(outputFolder,file = "all_dataGatedPop.csv"))

#to save into FSC file
plot_data = all_data_full
plot_data$Region = factor(plot_data$Region)
levels(plot_data$Region) = c('DG','CA2')
plot_data$brainRegion = paste0('p', plot_data$Point, '_',plot_data$Region)

plot_data$Region = as.numeric(plot_data$Region)
plot_data$Point = as.numeric(plot_data$Point)
plot_data$Meta = as.numeric(plot_data$Meta)
plot_data$brainRegion = as.numeric(as.factor(plot_data$brainRegion))

#to remove NA columns
plot_data = plot_data[ , colSums(is.na(plot_data)) == 0]

fcsExprs <- flowCore::flowFrame(as.matrix(plot_data))
suppressWarnings(flowCore::write.FCS(fcsExprs, paste0(outputFolder, 'alldataGatedPop.fcs'), what = "numeric"))
#_______________________________________________________________________________________________________________________________________#
#Read in csv
all_data_full = read_csv("/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/GatedPoP/NewMatlabSettings/4Gated0FlowSOMClustersNoOligo2/all_dataGatedPop.csv")

# START OF Make a barplot of cluster sizes
plot_data = all_data_full
plot_data = data.frame(table(plot_data$Meta))
#plot_data = subset(plot_data, Var1!="OtherNeurons")#____________________________________________________________________________________________________________

plot_data$Var1 = as.factor(plot_data$Var1)
levels(plot_data$Var1) = gateNames
levels(plot_data$Var1)

# check order of population : levels(plot_data$Var1)
plot_data$Var1 = ordered(plot_data$Var1, levels = c("microglia","astrocytes","endothelial","neurons"))
levels(plot_data$Var1)

## astrocytes Green "#106900", endothelial Blue"#0028D4",  microglia mango "#FFBF00",neurons deep Red "#E80009"...Color order follows order BEFORE custom ordering which is ("neurons""endothelial" "microglia" "astrocytes")
palette = c("#E80009","#0028D4","#FFBF00","#106900")

ggplot(plot_data, aes(x = Var1, y = Freq)) +
  geom_bar(aes(fill = Var1), stat = "identity", colour = "black", fill = palette) + 
  xlab("metaclusters") + ylab("Cluster size") +
  theme_classic(15) +
  theme(legend.position="none")
ggsave(filename = paste0(outputFolder,"/MetaClusterSizeGatedPoP3.pdf"), width = 20, height = 10)
#_______________________________________________________________________________________________________________________________________


#_______________________________________________________________________________________________________________________________________#
# Make a heatmap for channels vs clusters with mean expression values from all_data. Scaled by marker (col) version of the above
plot_data = all_data[, c(markers, "Meta")]

# to finD out number of cell in each condition
sum(!is.na(plot_data$Amyloidbeta140)) # number of AD cells

## to remove NA columns to make AD specific HeatMap
# plot_data = plot_data[!is.na(plot_data$Amyloidbeta140),] 
table(plot_data$Meta)

plot_data = melt(plot_data, id = c("Meta"))
plot_data = dcast(plot_data, variable ~ Meta, function(x){
  mean(x, na.rm = T )
})
plot_data$variable = NULL
plot_data = t(plot_data)
plot_data_rowNames = rownames(plot_data)

plot_data = apply(plot_data, 2, scale)#marker scaled
plot_data = rescale(plot_data, to = c(0, 1))#adjust between 0 and 2
colnames(plot_data) = markers
rownames(plot_data) = plot_data_rowNames
plot_data = plot_data[c(3,4,2,1),] # IN console to check order rownames(plot_data). Which should be #("microglia"   "astrocytes"  "endothelial" "neurons"))_____________________________________________________________________________________________________________________________

GatingMarkers = c("MAP2","MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% GatingMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_GatingMarkers2.pdf"), width = 7, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "Gated, Mean Marker Expression, \nGated Markers", 
           col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)

#breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

AllMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","MAP2","8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau","PanAmyloidbeta1724","Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1","pTDP43", "Synaptophysin","PSD95")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% AllMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_AllMarkers2.pdf"), width = 13, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, \nAll All Markers", 
           col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)          

#breakList = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

LineageMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","MAP2","Calbindin","Calretinin","Parvalbumin","Presenilin1")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% LineageMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_LineageMarkers2.pdf"), width = 13, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, \nPhenotypic Markers", 
           col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)          

#breakList = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

ADSpecificMarkers = c("8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","MFN2","PHF1Tau","PolyubiK48","PolyubiK63","pTDP43")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% ADSpecificMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_ADSpecificMarkers2.pdf"), width = 13, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, \nAD Specific Markers", 
           col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)          

#breakList = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

#_______________________________________________________________________________________________________________________________________#

#_______________________________________________________________________________________________________________________________________#
# Make a heatmap for channels vs clusters with SD (Z-scores) from all_data. 
plot_data = all_data[, c(markers, "Meta")]
plot_data = melt(plot_data, id = c("Meta"))
plot_data = dcast(plot_data, variable ~ Meta, sd)
plot_data$variable = NULL
plot_data = t(plot_data)
plot_data_rowNames = rownames(plot_data)
plot_data = apply(plot_data, 2, scale)#marker scaled
plot_data = rescale(plot_data, to = c(0, 1))#adjust between 0 and 2
colnames(plot_data) = markers
rownames(plot_data) = plot_data_rowNames
#plot_data = plot_data[c(1,3,2,4,5,6,7,8,9,10,11),]

dev.off()
pdf(paste0(outputFolder, "/FlowSOMClusters_ZScores.pdf"), width = 10, height = 10, onefile = F)

#breakList = c(0.0, 0.1, 0.11, 0.12, 0.14, 0.16, 0.18, 0.2,0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4,0.42, 0.44, 0.46, 0.48, 0.5,0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0)
#colfunc <- colorRampPalette(c("white","red"))

pheatmap(plot_data, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Z Scores",
         #color = viridis(length(plot_data),alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
         #col=colfunc(length(breakList)),breaks = breakList, cellwidth = 10, cellheight = 10)
         col= (colorRampPalette(brewer.pal(9,"Reds"))(50)[1:50]) , cellwidth = 10, cellheight = 10)

dev.off()
#_______________________________________________________________________________________________________________________________________#

#_______________________________________________________________________________________________________________________________________#
# DATA Prep for Violin plots
all_data_full = all_data
all_data_full$Meta = as.factor(all_data_full$Meta)
all_data_full$Meta = ordered(all_data_full$Meta, levels = c("microglia","astrocytes","endothelial","oligodendrocytes","neurons"))
all_data_sampled <- data.frame()
unifSize = 10000000
palette = c("#FF007F","#A634E4","#00CC66","#FF95FF","#078AFD") # For "neurons","endothelial","microglia","astrocytes","oligodendrocytes"

for (i in 1:numGates) {
  all_data_subset = all_data_full
  all_data_subset = subset(all_data_subset, Meta == levels(all_data_full$Meta)[i])
  
  if (unifSize < nrow(all_data_subset)) {
    unif = sample(1:nrow(all_data_subset), size = unifSize)
  } else {
    unif = 1:nrow(all_data_subset)
  }
  unif = all_data_subset[unif,]
  unif = melt(unif, id = c("Meta"))
  colnames(unif) = c("Meta", "Channel", "Values")
  
  all_data_sampled <<- rbind(all_data_sampled, unif)
}

#_______________________________________________________________________________________________________________________________________#

#_______________________________________________________________________________________________________________________________________#
# Violin plots per fSOM cluster for each channel
plot_list = list()
for (channel in markers) {
  plot_data = subset(all_data_sampled, Channel == channel)
  g = ggplot(plot_data, aes(x = as.factor(Meta), y = as.numeric(Values))) +
    geom_violin(aes(fill = Meta), scale = "width", width = 0.6) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
    # scale_y_continuous(trans=log2_trans(),
    #                    breaks = trans_breaks("log2", function(x) 2^x),
    #                    labels = trans_format("log2", math_format(2^.x)))+
    scale_fill_manual(values = palette) + 
    xlab("") + ylab("") +
    ggtitle(channel) +
    ylim(c(0,0.5)) + # watch for warnings to see how many outliers were removed in each plot.
    theme_classic(12) +
    theme(legend.position = "none")
  plot_list[[channel]] <- g
}

# Plot violin + boxplots for each channel separately. Show all clusters on each plot
picName = paste0(outputFolder, "/Violins_perChannelPHF1Tau.pdf")
cat(paste0("Plotting: ", picName, "\n"))
pdf(picName, width = 40, height = 60, onefile = T)
do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
dev.off()
# #_______________________________________________________________________________________________________________________________________#


# #_______________________________________________________________________________________________________________________________________#
# Violin plots per each channel for all fSOM clusters
plot_list = list()
violinMarkers = c("8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63") # c("8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63") c("PHF1Tau","ApoE4","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142")
#violinMarkers = markers #  = markers ## for all markers

for (meta in levels(all_data_full$Meta)) {
  #meta = c("microglia") #run if single population "microglia","astrocytes","endothelial","oligodendrocytes","inhibitory","excitatory","nonneuronal"
  plot_data = subset(all_data_sampled, Meta == meta)
  plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
  plot_data = subset(plot_data, Channel %in% violinMarkers)
  g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
    geom_violin(aes(fill = Channel), scale = "width", width = 0.6) +
    geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), outlier.size = -10) +
    scale_y_continuous(trans=log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)))+
    #limits = c(2^-10,2^10) +
    #ylim(c(2^-10,2^5)) +
    scale_fill_manual(values = colorRampPalette(c("#000000"))(10)[1:10]) + # colorRampPalette(c("deepskyblue4", "darkorange")), colorRampPalette(brewer.pal(9,"Blue"))
    #xlab("") + ylab("") +
    #ylim(0,2) +
    ggtitle(paste0("Cluster: ", meta)) +
    theme_classic(12) +
    theme(legend.position = "none") +
    theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
          axis.title.x = ggplot2::element_text(face = "bold", size = 10),
          axis.title.y = ggplot2::element_text(face = "bold", size = 10),
          axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
  
  # for single plot
  #ggsave(paste0(outputFolder, "/Violins_", meta,".pdf"),g, width = 6, height = 6)
  
  plot_list[[meta]] <- g
}

# Plot violin + boxplots for each channel separately. Show all clusters on each plot
picName = paste0(outputFolder, "/Violins_perCluster_InjuryMarksAllLog10.pdf")
cat(paste0("Plotting: ", picName, "\n"))
pdf(picName, width = 25, height = 5, onefile = T)
do.call(grid.arrange, list(grobs = plot_list, ncol = 5))
dev.off()
#_______________________________________________________________________________________________________________________________________#












#_______________________________________________________________________________________________________________________________________#
# # Make a stackedbarplot of cluster sizes for subset sampletypes 
# palette = c("#EE1289","#00008B","#FCF534","#FF7F00","#7850FF","#9932CC","#00FF33","#276C05","#4D1805","#FF0000", "#00FFFF")#ALL____________________________
# #this is to show region and point number
# plot_data = all_data_full
# plot_data$Region = factor(plot_data$Region)
# levels(plot_data$Region) = c('AD1','AD3')
# plot_data$brainRegion = paste0('p', plot_data$Point, '_',plot_data$Region)
# 
# ## plot_data = subset(plot_data, brainRegion %in% c("p1_PD", "p11_AD3" ,"p2_AD1" , "p2_PD"  , "p5_AD1" , "p5_AD3" , "p6_AD3" , "p7_AD1" , "p8_AD3" , "p9_AD3" ))
# ## sampleNames =c("SN_PD","HIP_AD","SN_H", "EC_PD","Str_H","MO_H","Crb_H","MTG_AD","LC_H","Hip_H" )   
# plot_data = subset(plot_data, brainRegion %in% c("p1_PD", "p11_AD3" ,"p2_AD1" , "p2_PD" , "p7_AD1","p9_AD3" ))
# sampleNames =c("SN_PD","HIP_AD","SN_H", "EC_PD","MTG_AD","Hip_H" )
# plot_data$brainRegion = as.factor(plot_data$brainRegion)
# levels(plot_data$brainRegion) = sampleNames [1:6] 
# plot_data$brainRegion = ordered(plot_data$brainRegion, levels = c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) #  in console to check if levels(plot_data$brainRegion)
# 
# plot_data = subset(plot_data, Meta %in% c("oligodendrocytes","astrocytes","microglia","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5"))          
# PopulationNames =c("oligodendrocytes","astrocytes","microglia","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5")
# plot_data$Meta = as.factor(plot_data$Meta)
# levels(plot_data$Meta) = PopulationNames [1:11] 
# plot_data$Meta = ordered(plot_data$Meta, levels = c("oligodendrocytes","microglia","astrocytes","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5")) # in console to check if levels(plot_data$Meta)
# 
# plot_data = data.frame(table(plot_data$brainRegion, plot_data$Meta))
# 
# #plot_data = subset(plot_data, Var2!="OtherNeurons")# TO EXCLUDE adjust color palette_________________________________________________________________________
# plot_data = subset(plot_data, Var2 %in% c("oligodendrocytes","astrocytes","microglia","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5"))# TO INCLUDE adjust color palette
# #plot_data = subset(plot_data, Var2 %in% c("TH","VGAT","FC1","FC2","FC3","FC4","FC5"))# TO INCLUDE adjust color palette
# #plot_data = subset(plot_data, Var2 %in% c("oligodendrocytes","astrocytes","microglia","endothelial")) # TO INCLUDE adjust color palette
# tempPlotData  = data.frame()
# 
# for (i in unique(plot_data$Var1)){
#   tempdata = subset(plot_data, Var1 == i)
#   sumCellCounts = sum(tempdata$Freq) # total number of cells in the brain region
#   tempdata$Freq = (tempdata$Freq/sumCellCounts)*100 #convert to fractions
#   tempPlotData  = rbind(tempPlotData , tempdata)
# }
# plot_data = tempPlotData 
# 
# ggplot(plot_data, aes(x=Var1, y=Freq, fill = Var2)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = palette[1:11]) +
#   labs(x="Brain Region", y="Percentage") +
#   theme_bw(15)
# ggsave(filename= paste0(outputFolder,"/", "StackedBarPlot_Neurons.pdf"), width = 4, height = 6)
# #_______________________________________________________________________________________________________________________________________#
# 
# 
# 
# 
# #_______________________________________________________________________________________________________________________________________#
# # DATA prep for Violin plots per brain region for specific channel
# plot_data = all_data_full
# plot_data$Region = factor(plot_data$Region)
# levels(plot_data$Region) = c('AD1','AD3', 'PD')
# plot_data$brainRegion = paste0('p', plot_data$Point, '_',plot_data$Region)
# 
# # plot_data = subset(plot_data, brainRegion %in% c("p1_PD", "p11_AD3" ,"p2_AD1" , "p2_PD"  , "p5_AD1" , "p5_AD3" , "p6_AD3" , "p7_AD1" , "p8_AD3" , "p9_AD3" ))
# # sampleNames =c("SN_PD","HIP_AD","SN_H", "EC_PD","Str_H","MO_H","Crb_H","MTG_AD","LC_H","Hip_H" )   
# plot_data = subset(plot_data, brainRegion %in% c("p1_PD", "p11_AD3" ,"p2_AD1" , "p2_PD" , "p7_AD1","p9_AD3" ))
# sampleNames =c("SN_PD","HIP_AD","SN_H", "EC_PD","MTG_AD","Hip_H" )
# plot_data$brainRegion = as.factor(plot_data$brainRegion)
# levels(plot_data$brainRegion) = sampleNames [1:6] 
# plot_data$brainRegion = ordered(plot_data$brainRegion, levels = c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) # levels(plot_data$brainRegion) in console to check if 
# 
# plot_data = subset(plot_data, Meta %in% c("oligodendrocytes","astrocytes","microglia","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5"))          
# PopulationNames =c("oligodendrocytes","astrocytes","microglia","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5")
# plot_data$Meta = as.factor(plot_data$Meta)
# levels(plot_data$Meta) = PopulationNames [1:11] 
# plot_data$Meta = ordered(plot_data$Meta, levels = c("oligodendrocytes","microglia","astrocytes","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5")) # in console to check if levels(plot_data$Meta)
# 
# #plot_data = subset(plot_data, Meta!="OtherNeurons")# TO EXCLUDE adjust color palette_________________________________________________________________________
# plot_data = subset(plot_data, Meta %in% c("oligodendrocytes","astrocytes","microglia","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5"))# TO INCLUDE adjust color palette
# #plot_data = subset(plot_data, Meta %in% c("TH","VGAT","FC1","FC2","FC3","FC4","FC5"))# TO INCLUDE adjust color palette
# #plot_data = subset(plot_data, Meta %in% c("oligodendrocytes","astrocytes","microglia","endothelial")) # TO INCLUDE adjust color palette
# plot_data$Region = NULL
# plot_data$Point = NULL
# 
# all_data_sampled = plot_data
# all_data_sampled = melt(all_data_sampled, id = c("Meta", "brainRegion"))
# colnames(all_data_sampled) = c("Meta", "brainRegion", "Channel", "Values") 
# 
# 
# #_______________________________________________________________________________________________________________________________________#
# # Violin plots per brain region for specific channel.......NONNEURONAL
# specificChannel = "PHF" # "pTDP43","8OHGuano", "beta142", "ApoE", "PolyubiK63", "PanAbeta","PHF"__________________________________________
# plot_list = list()
# for (BR in c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) {
#   plot_data = subset(all_data_sampled, Channel == specificChannel & brainRegion == BR)
#   g = ggplot(plot_data, aes(x = Meta, y = Values)) +
#     geom_violin(aes(fill = Meta), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0)
#   if (BR != "Hip_H" & BR != "EC_PD")#
#     g = g + scale_fill_manual(values=palette)
#   if (BR == "Hip_H")#
#     g = g + scale_fill_manual(values=palette[2:9])  # this is custom for the Hip_H and this TMA
#   if (BR == "EC_PD")#
#     g = g + scale_fill_manual(values=palette[c(1:4,6:11)])  # this is custom for the EC_PD and this TMA #
#   g = g + xlab("") + ylab("") +
#     ggtitle(BR) +
#     ylim(c(0, 1.5)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic() +
#     theme(legend.position = "none")
#   plot_list[[BR]] <- g
# }
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_NonNeuronsSingle_perBR_", specificChannel, ".pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 15, height = 8, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 6))
# dev.off()
# 
# ##Violin plots per brain region for specific channel.......NEURONAL.MULTIPLE CHANNELS Disease Marks
# specificChannel = c("PHF","PanAbeta") # "pTDP43","8OHGuano", "beta142", "ApoE", "PolyubiK63", "PanAbeta","PHF"__________________________________________
# plot_list = list()
# paletteChannelDisease = c("#FFDE17", "#00CCCC")
# for (BR in c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) {
#   plot_data = subset(all_data_sampled, Channel %in% specificChannel & brainRegion == BR)
#   g = ggplot(plot_data, aes(x = Meta, y = Values)) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) 
#   #geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0, outlier.shape = NA) 
#   g = g + scale_fill_manual(values=paletteChannelDisease)
#   g = g + xlab("") + ylab("") +
#     ggtitle(BR) +
#     ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic() +
#     theme(legend.position = "none")
#   plot_list[[BR]] <- g
# }
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_NeuronsMultipleDisease_perBR_", specificChannel, ".pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 20, height = 2, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 6))
# dev.off()
# 
# ##Violin plots per brain region for specific channel.......NEURONAL.MULTIPLE CHANNELS Injury Marks
# specificChannel = c("8OHGuano","pTDP43") # "pTDP43","8OHGuano", "beta142", "ApoE", "PolyubiK63", "PanAbeta","PHF"__________________________________________
# plot_list = list()
# paletteChannelInjury = c("#D7DF23", "#BE1E2D")
# for (BR in c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) {
#   plot_data = subset(all_data_sampled, Channel %in% specificChannel & brainRegion == BR)
#   g = ggplot(plot_data, aes(x = Meta, y = Values)) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) 
#   #geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0, outlier.shape = NA) 
#   g = g + scale_fill_manual(values=paletteChannelInjury)
#   g = g + xlab("") + ylab("") +
#     ggtitle(BR) +
#     ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic() +
#     theme(legend.position = "none")
#   plot_list[[BR]] <- g
# }
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_NeuronsMultipleInjury_perBR_", specificChannel, ".pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 20, height = 2, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 6))
# dev.off()
# #_______________________________________________________________________________________________________________________________________#
# 
# 
# #_______________________________________________________________________________________________________________________________________#
# ##Violin plots per brain region for specific channel.......NEURONAL....INDIVIDUAL CHANNELS
# specificChannel = "PHF" # "pTDP43","8OHGuano", "beta142", "ApoE", "PolyubiK63", "PanAbeta","PHF"__________________________________________
# plot_list = list()
# for (BR in c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) {
#   plot_data = subset(all_data_sampled, Channel == specificChannel & brainRegion == BR)
#   g = ggplot(plot_data, aes(x = Meta, y = Values)) +
#     geom_violin(aes(fill = Meta), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0)
#   if (BR != "Hip_H" & BR != "EC_PD" & BR != "SN_PD")# general case where it works
#     g = g + scale_fill_manual(values=palette[5:9])
#   if (BR == "Hip_H")#
#     g = g + scale_fill_manual(values=palette[c(5,6,7,8,9)])  # this is custom for those that dont work
#   if (BR == "EC_PD")#
#     g = g + scale_fill_manual(values=palette[c(6,7,8,9)])  #  this is custom for those that dont work
#   if (BR == "SN_PD")#
#     g = g + scale_fill_manual(values=palette[c(6,7,8,9)])  #  this is custom for those that dont work
#   g = g + xlab("") + ylab("") +
#     ggtitle(BR) +
#     ylim(c(0, 1.2)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_bw(15) +
#     theme(legend.position = "none")
#   plot_list[[BR]] <- g
# }
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_NeuronsSingle_perBR_", specificChannel, ".pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 15, height = 8, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 6))
# dev.off()
# 
# ##Violin plots per brain region for specific channel.......NEURONAL.MULTIPLE CHANNELS Disease Marks
# specificChannel = c("PHF","PanAbeta") # "pTDP43","8OHGuano", "beta142", "ApoE", "PolyubiK63", "PanAbeta","PHF"__________________________________________
# plot_list = list()
# paletteChannelDisease = c("#FFDE17", "#00FFFF")
# for (BR in c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) {
#   plot_data = subset(all_data_sampled, Channel %in% specificChannel & brainRegion == BR)
#   g = ggplot(plot_data, aes(x = Meta, y = Values)) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) 
#   #geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0, outlier.shape = NA) 
#   g = g + scale_fill_manual(values=paletteChannelDisease)
#   g = g + xlab("") + ylab("") +
#     ggtitle(BR) +
#     ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic() +
#     theme(legend.position = "none")
#   plot_list[[BR]] <- g
# }
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_NeuronsMultipleDisease_perBR_", specificChannel, ".pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 20, height = 2, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 6))
# dev.off()
# 
# ##Violin plots per brain region for specific channel.......NEURONAL.MULTIPLE CHANNELS Injury Marks
# specificChannel = c("8OHGuano","pTDP43") # "pTDP43","8OHGuano", "beta142", "ApoE", "PolyubiK63", "PanAbeta","PHF"__________________________________________
# plot_list = list()
# paletteChannelInjury = c("#D7DF23", "#BE1E2D")
# for (BR in c("SN_H","SN_PD","EC_PD","Hip_H","HIP_AD","MTG_AD")) {
#   plot_data = subset(all_data_sampled, Channel %in% specificChannel & brainRegion == BR)
#   g = ggplot(plot_data, aes(x = Meta, y = Values)) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) 
#   #geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0, outlier.shape = NA) 
#   g = g + scale_fill_manual(values=paletteChannelInjury)
#   g = g + xlab("") + ylab("") +
#     ggtitle(BR) +
#     ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic() +
#     theme(legend.position = "none")
#   plot_list[[BR]] <- g
# }
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_NeuronsMultipleInjury_perBR_", specificChannel, ".pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 20, height = 2, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 6))
# dev.off()
# #_______________________________________________________________________________________________________________________________________#
# 
# 
# 
# #_______________________________________________________________________________________________________________________________________#
# # DATA Prep for Violin plots per fSOM cluster for each channel
# 
# all_data_sampled <- data.frame()
# unifSize = 10000000
# for (i in 1:numClusters) {
#   all_data_subset = all_data_full
#   all_data_subset = subset(all_data_subset, Meta == levels(all_data_full$Meta)[i])
#   
#   if (unifSize < nrow(all_data_subset)) {
#     unif = sample(1:nrow(all_data_subset), size = unifSize)
#   } else {
#     unif = 1:nrow(all_data_subset)
#   }
#   unif = all_data_subset[unif,]
#   unif = melt(unif, id = c("Meta"))
#   colnames(unif) = c("Meta", "Channel", "Values")
#   
#   all_data_sampled <<- rbind(all_data_sampled, unif)
# }
# 
# # Violin plots per fSOM cluster for each channel
# plot_list = list()
# for (channel in markers) {
#   plot_data = subset(all_data_sampled, Channel == channel)
#   plot_data$Meta = ordered(plot_data$Meta, levels = c("oligodendrocytes","microglia","astrocytes","endothelial","TH","VGAT","FC1","FC2","FC3","FC4","FC5")) # in console to check if levels(plot_data$Meta)
#   g = ggplot(plot_data, aes(x = as.factor(Meta), y = as.numeric(Values))) +
#     geom_violin(aes(fill = Meta), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
#     scale_fill_manual(values=palette[1:11]) +#_____________________________________________________________________________________________________________________
#     xlab("") + ylab("") +
#     ggtitle(channel) +
#     ylim(c(0, 2)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_bw(15) +
#     theme(legend.position = "none")
#   plot_list[[channel]] <- g
# }
# 
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_perChannel.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 40, height = 60, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
# dev.off()
# #_______________________________________________________________________________________________________________________________________#
# 
# 
# 
# #_______________________________________________________________________________________________________________________________________#
# # DATA Prep for Violin plots per each channel for all fSOM clusters
# all_data_sampled <- data.frame()
# unifSize = 10000000
# for (i in 1:numClusters) {
#   all_data_subset = all_data_full
#   all_data_subset = subset(all_data_subset, Meta == levels(all_data_full$Meta)[i])
#   
#   if (unifSize < nrow(all_data_subset)) {
#     unif = sample(1:nrow(all_data_subset), size = unifSize)
#   } else {
#     unif = 1:nrow(all_data_subset)
#   }
#   unif = all_data_subset[unif,]
#   unif = melt(unif, id = c("Meta"))
#   colnames(unif) = c("Meta", "Channel", "Values")
#   
#   all_data_sampled <<- rbind(all_data_sampled, unif)
# }
# # Violin plots per each channel for all fSOM clusters
# plot_list = list()
# violinMarkers = c("TH", "VGAT", "VGLUT1", "VGLUT2", "PSD95", "Synaptophysin")
# plot_data = subset(plot_data, Channel %in% violinMarkers)
# 
# for (meta in levels(all_data_full$Meta)) {
#   meta = "OtherNeurons" #run if single population
#   plot_data = subset(all_data_sampled, Meta == meta)
#   plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
#   #plot_data = subset(plot_data, Channel %in% markers)
#   plot_data = subset(plot_data, Channel %in% violinMarkers)
#   g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
#     scale_fill_manual(values = colorRampPalette(c("#4D1805"))(10)[1:10]) + # colorRampPalette(c("deepskyblue4", "darkorange")), colorRampPalette(brewer.pal(9,"Blue"))
#     xlab("") + ylab("") +
#     ylim(0,1.0) +
#     ggtitle(paste0("Cluster: ", meta)) +
#     theme_bw(12) +
#     theme(legend.position = "none") +
#     theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.x = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.y = ggplot2::element_text(face = "bold", size = 10),
#           axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
#   # for single plot 
#   ggsave(paste0(outputFolder, "/Violins2_", meta,".pdf"),g, width = 6, height = 6)
#   plot_list[[meta]] <- g
# }
# 
# # Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_perCluster.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 40, height = 2.5*numClusters, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
# dev.off()
#_______________________________________________________________________________________________________________________________________#






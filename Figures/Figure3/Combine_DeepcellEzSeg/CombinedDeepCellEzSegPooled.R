libs <- c('ggplot2','RColorBrewer','reshape2','devtools','readr',
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
require(data.table)
setwd("/Volumes/KausaliaHD/R4R")
source("/Volumes/KausaliaHD/R4R/General/readWriteCytofFCS.R")


mantisPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/"
ezPath ="/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/"

outputFolder = paste0(mantisPath,"CombineGatedMantisEz/ForDiseaseObjectCellEngineOverlay/")
dir.create(outputFolder, showWarnings= F, recursive = T)

##NewMatlabSettings: Deepcell
all_data_mod <- read_csv(paste0(mantisPath,"GatedPopWithMantisCuration/NewMatlabSettings/UMAP/umap_data_sampled.csv"))

##Import csv deepcell 
table(all_data_mod$Point)
all_data_mod$point_id = all_data_mod$Point
all_data_mod$obj_id = all_data_mod$cellLabelInImage
all_data_mod$obj_size = all_data_mod$cellSize


##Pool csv EzSEG
all_data_EzSeg = read_csv (paste0(ezPath, "EzSegData_uci2712J/all_dataEzSeg.csv"))

table(all_data_EzSeg$point_id)

allMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","Reelin","MAP2","X8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau","Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1NTF","pTDP43","EEA1","CD56Lyo","CD47", "TotalTau", "PanGAD6567","CD33Lyo","HistoneH3Lyo", "Synaptophysin", "PSD95", "PanApoE2E3E4")
commonCols = c(allMarkers, "Meta","point_id","obj_id", "obj_size")

all_data_EzSeg$Meta = all_data_EzSeg$obj_type_id
all_data_mod$Meta = all_data_mod$MantisPopulation
all_data_mod$X8OHGuano = all_data_mod$`8OHGuano` 
#all_data_EzSeg$X8OHGuano = all_data_EzSeg$`8OHGuano` 

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
    ##all_data_EzSeg_subset = all_data_EzSeg_subset[sample(1:nrow(all_data_EzSeg_subset),ceiling(nrow(all_data_EzSeg_subset)/10)),] # to give a sample size
    all_data_EzSeg_sampled = rbind(all_data_EzSeg_sampled, all_data_EzSeg_subset)
  } else {
    if(i %in% c("microglia_process")){
      ##all_data_EzSeg_subset = all_data_EzSeg_subset[sample(1:nrow(all_data_EzSeg_subset),ceiling(nrow(all_data_EzSeg_subset)/5.56)),] # to give a sample size
      all_data_EzSeg_sampled = rbind(all_data_EzSeg_sampled, all_data_EzSeg_subset)
    } else {
      
      all_data_EzSeg_sampled = rbind(all_data_EzSeg_sampled, all_data_EzSeg_subset)
      
    }
  }
}


##this is to save combined fsc file for CellEngine
all_data_mod$Source = 0
all_data_EzSeg_sampled$Source = 1

plot_data = rbind(all_data_mod[,c(commonCols, "Source")], all_data_EzSeg_sampled[,c(commonCols, "Source")])
table(plot_data$Meta)

plot_data$Meta = as.numeric(as.factor(plot_data$Meta))
table(plot_data$Meta)

#to remove NA columns
plot_data = plot_data[ , colSums(is.na(plot_data)) == 0]
fcsExprs <- flowCore::flowFrame(as.matrix(plot_data))
suppressWarnings(flowCore::write.FCS(fcsExprs, paste0(outputFolder, "alldataDiseasePooled.fcs"), what = "numeric"))

##this is to save combined csv file for whatever your heart desires
## if you want to write out the merged data as a csv:
fwrite(plot_data, paste0(outputFolder, "alldataDiseasePooled2.csv"))


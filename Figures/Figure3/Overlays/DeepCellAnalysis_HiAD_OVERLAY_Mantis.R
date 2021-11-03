library(EBImage)
library(nnet)
library(reshape2)
library(readr)
setwd("/Volumes/KausaliaHD 1/R4R/General/")

#Run DeepCellAnalysis.R then DeepCellAnalysisUmapTMAQuntile before using this Rscript to get "gated_UMapData" or "all_data" clusters!!!!!!

#DeepCell Segmentation Paths
fcsPath =  "/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSegData_uci2712J/DGCA2Pooled/"
# fcsRuns = c("190505HiResDG/", "190604HiResCA2/")
# fcsRunsShort = c("DG", "CA2")
fcsPixel = "single_cell_dynamic_expansion_Smaller_no_fftnoiseDist/Point"

# Mantis Paths (Dmitry's debugging)
#all_data_mod <- read_csv(paste0(fcsPath, "MantisData/MantisDataPooledOverlayNoOligoResults/all_data_mod.csv"))
mantisPath = "/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/"

# #all_data_mod <- read_csv(paste0(mantisPath, "GatedPopWithMantisCuration/NewMatlabSettings/umap_data_sampledUMAPCurated.csv"))
all_data_mod <- read_csv(paste0(mantisPath, "GatedPopWithMantisCuration/NewMatlabSettings/AfterMantisCuration(26Oct20)/all_data_modMore.csv"))

outputFolder = paste0(mantisPath,"Curation/OverlayAll_White/DG")
dir.create(outputFolder, showWarnings= F, recursive = T)

outputFolder = paste0(mantisPath,"Curation/OverlayAll_White/CA2")
dir.create(outputFolder, showWarnings= F, recursive = T)

outputFolder = paste0(mantisPath,"Curation/OverlayAll_White/")

all_data = all_data_mod # reads in Mantis Corrected Population
#all_data = umap_data_sampledCurated # reads in UMAP & Mantis Corrected Population
all_data$OriginalMeta = all_data$Meta
all_data$Meta = as.factor(all_data$MantisPopulation)

##To Import CSV  

all_dataDG = subset(all_data, Region == "DG")
all_dataDG$Point = all_dataDG$Point + 34
all_dataCA2 = subset(all_data, Region == "CA2")
all_data =rbind(all_dataCA2 , all_dataDG)
table(all_data$Point)


# Run these 2 lines if SampleNames column does not exist in all_data
# all_data$SampleNames = as.factor(all_data$brainRegion)
# levels(all_data$SampleNames) = sampleNames

table(all_data$Meta) #to order colours
table(all_data$MantisPopulation)
table(all_data[,c("MantisPopulation","brainRegion")])

## astrocytes Green "#106900", endothelial Blue "#0028D4",  microglia mango "#FFBF00",neurons deep red "#980A0A"...Color order follows all_data$Meta order 
palette = c("#106900","#0028D4","#FFBF00","#980A0A") #All astrocytes endothelial   microglia     neurons
# palette = c("#106900","#000000","#000000","#000000") #Astro
# palette = c("#000000","#0028D4","#000000","#000000") #Endo
# palette = c("#000000","#000000","#FFBF00","#000000") #Micro
# palette = c("#106900","#000000","#FFBF00","#000000") #Glia
#  palette = c("#106900","#0028D4","#000000","#000000") #BBB
# palette = c("#000000","#000000","#000000","#980A0A") #Neuro

 ## White backgroud
 # palette = c("#106900","#FFFFFF","#FFFFFF","#FFFFFF") #Astro
 # palette = c("#FFFFFF","#0028D4","#FFFFFF","#FFFFFF") #Endo
 # palette = c("#FFFFFF","#FFFFFF","#FFBF00","#FFFFFF") #Micro
 # palette = c("#106900","#FFFFFF","#FFBF00","#FFFFFF") #Glia
 # palette = c("#106900","#0028D4","#FFFFFF","#FFFFFF") #BBB
# palette = c("#FFFFFF","#FFFFFF","#FFFFFF","#980A0A") #Neuro
 # palette = c("#FFFFFF","#0028D4","#FFFFFF","##980A0A") #Neuro_Endo
 # palette = c("#106900","#0028D4","#FFFFFF","#FFFFFF") #Astro_Endo
 #  palette = c("#FFFFFF","#0028D4","#FFBF00","#FFFFFF") #Micro_Endo


# for (runNumber in 1:length(fcsRunsShort)) {
#   region = fcsRuns[runNumber]
 
pointList = unique(all_data$Point)
for (p in pointList) {  # comment out to run one point 
    #p = c("18")#,"15","30") # comment in to run one point
    print(paste0("Processing Point: ", p))
    newLmod <- read.csv(paste0(fcsPath, fcsPixel, p, "/newLmod.csv"), header = F)
    overlayImage = newLmod
    #cell cluster assignments for each point
    cellClustersAssign = data.frame("cellid" = character(), "population" = character())
    
    # Set background to zero from whatever index it is. Assume: background is the most prevalent object (most pixels) in the image
    bg_finder = table(reshape2::melt(overlayImage, id.vars = NULL)$value)
    bg_index = which.is.max(bg_finder)
    bg_cellLabel = as.numeric(names(bg_finder[bg_index]))
    overlayImage[overlayImage == bg_cellLabel] = 0
    
    overlayImage_r = overlayImage
    overlayImage_g = overlayImage
    overlayImage_b = overlayImage
    
    #for each cell in Point p we replace cell ID with colour based on meta number. Change between all_data or gated_UMapData.
    #pointdata =  subset(gated_UMapData, Region == region & Point == p)
    #pointdata =  subset(all_data, Region == fcsRunsShort[i] & Point == p) 
    pointdata =  subset(all_data, Point == p)
    #print(paste0("p",p, "_",fcsRunsShort[runNumber]))
    #pointdata =  subset(point_data, Meta!= "MIBIFUA")# remove 8
   
    #warning to say cell data foundin alldata or point empty
    if (nrow(pointdata) == 0){
      print("no cell data found")
      next
    }
    
    for (cell in unique(reshape2::melt(overlayImage, id.vars = NULL)$value)) {
      #print(paste0("Processing region: ", region, ". cell: ", cell))
      cellAssign = data.frame("cellid" = character(), "population" = character())
      if (cell %in% unique(pointdata$cellLabelInImage)) {
        cellMetas = subset(pointdata, cellLabelInImage == cell)$Meta
        
        if (length(cellMetas) == 1) {
          # Cell has a single assignment
          overlayImage_r[overlayImage == cell] = col2rgb(palette)[1, cellMetas] 
          overlayImage_g[overlayImage == cell] = col2rgb(palette)[2, cellMetas] 
          overlayImage_b[overlayImage == cell] = col2rgb(palette)[3, cellMetas] 
          
        } else if (length(cellMetas) <= 4) {
          # For now - up to 4 Meta assignments, can do more splits in the future if needed
          if (length(cellMetas) == 2) {
            cellMetas = c(cellMetas[1], cellMetas[2], cellMetas[2], cellMetas[1])
          } 
          
          if(length(cellMetas) == 3) {
            cellMetas = c(cellMetas[1], cellMetas[2], cellMetas[3], cellMetas[1])
          } 
          
          cellLocation = data.frame(which(newLmod == cell, arr.ind = TRUE))
          cellCenterRow = ceiling(mean(cellLocation$row))
          cellCenterCol = ceiling(mean(cellLocation$col))
          
          cellLoc = list()
          
          # First quadrant - top left
          cellLoc[[1]] = subset(cellLocation, row %in% c(min(cellLocation$row) : cellCenterRow) & col %in% c(min(cellLocation$col) : cellCenterCol))
          
          # Second quadrant - top right
          cellLoc[[2]] = subset(cellLocation, row %in% c(min(cellLocation$row) : cellCenterRow) & col %in% c((cellCenterCol + 1) : max(cellLocation$col)))
          
          # Third quadrant - bottom left
          cellLoc[[3]] = subset(cellLocation, row %in% c((cellCenterRow + 1) : max(cellLocation$row)) & col %in% c(min(cellLocation$col) : cellCenterCol))
          
          # Fourth quadrant - bottom right
          cellLoc[[4]] = subset(cellLocation, row %in% c((cellCenterRow + 1) : max(cellLocation$row)) & col %in% c((cellCenterCol + 1) : max(cellLocation$col)))
          
          for (i in 1:4) {
            for (j in 1:nrow(cellLoc[[i]])) {
              overlayImage_r[cellLoc[[i]]$row[j], cellLoc[[i]]$col[j]] = col2rgb(palette)[1, cellMetas[i]] 
              overlayImage_g[cellLoc[[i]]$row[j], cellLoc[[i]]$col[j]] = col2rgb(palette)[2, cellMetas[i]] 
              overlayImage_b[cellLoc[[i]]$row[j], cellLoc[[i]]$col[j]] = col2rgb(palette)[3, cellMetas[i]] 
            }
          }
        } else {
          # If cell has more than 4 assignments - make it white
          overlayImage_r[overlayImage == cell] = 255
          overlayImage_g[overlayImage == cell] = 255
          overlayImage_b[overlayImage == cell] = 255
        }
        cellAssign = data.frame("cellid" = as.character(cell), "population" = as.character(subset(pointdata, cellLabelInImage == cell)$Meta))
      } else {
        # Set background and non-cell objects to black
        overlayImage_r[overlayImage == cell] = 255
        overlayImage_g[overlayImage == cell] = 255
        overlayImage_b[overlayImage == cell] = 255
      }
      
      if (nrow(cellAssign) > 0) {
        cellClustersAssign = rbind(cellClustersAssign, cellAssign) 
      }
    }
    #print(paste0("saving"))
    overlayImage_r = overlayImage_r / 255
    overlayImage_g = overlayImage_g / 255
    overlayImage_b = overlayImage_b / 255
    
    img <- EBImage::transpose(rgbImage(Image(as.matrix(overlayImage_r)), Image(as.matrix(overlayImage_g)), Image(as.matrix(overlayImage_b))))
    if(p <= 34){
      sampleName = paste0("CA2/p", p, "_CA2")
    } else {
      sampleName = paste0("DG/p", (p-34), "_DG")
    }
    writeImage(img, paste0(outputFolder,'/', sampleName,".tif"), bits.per.sample = 16) # If want smaller TIFs, swap to 8, hopefully doesn't break anything else
    write.csv(cellClustersAssign, paste0(outputFolder,'/', sampleName,".csv"), row.names = F)
}




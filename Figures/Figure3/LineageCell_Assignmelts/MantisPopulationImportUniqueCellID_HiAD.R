##RUN DeepCellAnalysis_HiAD_Gated_NormAllChannel.R to generate all_data ###
# Import existing all_data from the enviroment using Import Dataset. Rename simply into all_data
# Need to pool all the new Mantis corrected csv into one folder.

library(readr)

dataPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/"

outputFolder = paste0(dataPath,"/GatedPopWithMantisCuration/NewMatlabSettings/")
dir.create(outputFolder, showWarnings= F, recursive = T)

importPath = paste0(dataPath,"/Curation/OverlayAll/")

listOfFileNames = list.files(path = importPath, pattern = ".*ModO\\.csv$", full.names = F) 
listOfFileNames = gsub("\\.csv", "", listOfFileNames)
listOfFileNames = gsub("ModO", "", listOfFileNames)
table(all_data$Point)
table(all_data$SampleNames)
## check list of listOfFileNames , need to have same data as point names point colmn the spead all_data sheet 

#need to run original scipt for all_data
all_data_mod = data.frame()

for (filename in listOfFileNames) {
  print(filename)
  tempdf =  read_csv(paste0(importPath,filename, "ModO.csv"))
  all_data_mod_subset = subset(all_data, SampleNames == filename) # if doing single run replace SampleNames with Point
  all_data_mod_subset$MantisPopulation = NA
  
  for (cell_id in unique(all_data_mod_subset$cellLabelInImage)){
    cell = subset(tempdf, cellid == cell_id)
    if (nrow(cell) == 0) {
      next
    }
    if (nrow(cell) == 1) {
      all_data_mod_subset[which(all_data_mod_subset$cellLabelInImage == cell_id), "MantisPopulation"] = cell$population
    } else {
      all_data_mod_subset[which(all_data_mod_subset$cellLabelInImage == cell_id), "MantisPopulation"] = cell$population[sample(1:nrow(cell),1)]
    }
  }
  all_data_mod = rbind(all_data_mod, all_data_mod_subset)
}

all_data_mod$NumericCode = as.numeric(as.factor(all_data_mod$MantisPopulation))


write.csv(all_data_mod, paste0(outputFolder, "/all_data_mod_unique_cell_id.csv"))

# START OF Make a barplot of cluster sizes
plot_data = all_data_mod
plot_data = data.frame(table(plot_data$MantisPopulation))

# Cell&Population Distribution Before Mantis Curation
table(all_data_mod$Meta) 
table(all_data_mod[,c("Meta","SampleNames")])

# Cell&Population Distribution After Mantis Curation
table(all_data_mod$MantisPopulation) 
table(all_data_mod[,c("MantisPopulation","SampleNames")])

#plot_data = subset(plot_data, Var1!="Ungated")#____________________________________________________________________________________________________________

plot_data$Var1 = as.factor(plot_data$Var1)
levels(plot_data$Var1) = gateNames
levels(plot_data$Var1)

plot_data$Var1 = ordered(plot_data$Var1, levels = c("microglia","astrocytes","endothelial","neurons"))  # in console to check if levels(plot_data$Var1)
levels(plot_data$Var1)

## astrocytes Cyan"#00FFFF", endothelial Blue"#004C99",  microglia Yellow "#E8D316",neurons Fuchia "#FF007F"...Color order follows order before custome ordering
palette = c("#FF007F","#004C99","#E8D316","#00FFFF")

ggplot(plot_data, aes(x = Var1, y = Freq)) +
  geom_bar(aes(fill = Var1), stat = "identity", colour = "black", fill = palette) + 
  xlab("metaclusters") + ylab("Cluster size") +
  theme_classic(15) +
  theme(legend.position="none")
ggsave(filename = paste0(outputFolder,"/MetaClusterSizeGatedMantisCurated.pdf"), width = 20, height = 10)

# #_______________________________________________________________________________________________________________________________________#
# # Make a heatmap for channels vs clusters with mean expression values from all_data. Scaled by marker (col) version of the above
# plot_data = all_data[, c(markers, "Meta")]
plot_data = all_data_mod
plot_data = plot_data[, c(markers, "MantisPopulation")]

# to finD out number of cell in each condition
sum(!is.na(plot_data$Amyloidbeta140)) # number of AD cells

plot_data = melt(plot_data, id = c("MantisPopulation"))
plot_data = dcast(plot_data, variable ~ MantisPopulation, function(x){
  mean(x, na.rm = T )
})
plot_data$variable = NULL
plot_data = t(plot_data)
plot_data_rowNames = rownames(plot_data)

plot_data = apply(plot_data, 2, scale)#marker scaled
plot_data = rescale(plot_data, to = c(0, 1))#adjust between 0 and 2
colnames(plot_data) = markers
rownames(plot_data) = plot_data_rowNames
plot_data = plot_data[c(3,2,4),]# 5 is TH . IN console to check order rownames(plot_data). Which should be #("microglia","astrocytes","endothelial", neurons"))_____________________________________________________________________________________________________________________________

GatingMarkers = c("MAP2","Iba1","CD45","GFAP","CD31","MCT1","CD105")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% GatingMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_GatingMarkers.pdf"), width = 7, height = 7, onefile = F)
breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 

pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, Gated Markers", 
         #  color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
         col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)
dev.off()

AllMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","MAP2","8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau","PanAmyloidbeta1724","Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1","pTDP43")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% AllMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_AllMarkers.pdf"), width = 13, height = 7, onefile = F)
breakList = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, All All Markers", 
         #  color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
         col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)
dev.off()

LineageMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","VGLUT1","VGLUT2","MAP2","Calbindin","Calretinin","Parvalbumin","Presenilin1NTF")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% LineageMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_LineageMarkers.pdf"), width = 13, height = 7, onefile = F)
breakList = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, Lineage Markers", 
         #  color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
         col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)
dev.off()

NeuralMarkers = c("VGAT","VGLUT1","VGLUT2","MAP2","Calbindin","Calretinin","Parvalbumin","Presenilin1NTF","PanGAD6567")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% NeuralMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_NeuralMarkers4.pdf"), width = 13, height = 7, onefile = F)
breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.32,0.34,0.46,0.38,0.4,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.8,0.9,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, Neural Markers", 
         #  color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
         col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)
dev.off()

ADSpecificMarkers = c("8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","MFN2","PHF1Tau","PolyubiK48","PolyubiK63","pTDP43")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% ADSpecificMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_ADSpecificMarkers.pdf"), width = 13, height = 7, onefile = F)
breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.8,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = F, border_color = FALSE, scale = "none", main = "FlowSOM clusters, Mean Marker Expression, AD Specific Markers", 
         #  color = viridis(length(plot_data),  alpha=1, begin=0, end=1, direction=1, option = "B"), cellwidth = 10, cellheight = 10)
         col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)
dev.off()
#_______________________________________________________________________________________________________________________________________#





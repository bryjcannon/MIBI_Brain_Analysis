# Run DeepCellAnalysisTMAQuantile.R to get all_data
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
library(premessa)
setwd("/Volumes/KausaliaHD 1/R4R")
source("/Volumes/KausaliaHD 1/R4R/General/readWriteCytofFCS.R")


dataPath = "/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/GatedPopWithMantisCuration/NewMatlabSettings/AfterMantisCuration(26Oct20)/"
outputFolder = paste0(dataPath,"FlowSOM30_HiAD_Neurons(26Oct20)")
dir.create(outputFolder, showWarnings= F, recursive = T)

markers = c("HistoneH3Lyo",
            "TotalTau","CD56Lyo","CD47","Synaptophysin","PSD95","MAP2",
            "PanGAD6567","VGAT","VGLUT1","VGLUT2","Reelin",
            "Presenilin1NTF","Calbindin","Calretinin","Parvalbumin",
            "PanApoE2E3E4","ApoE4","Iba1","CD45","CD33Lyo","GFAP","CD31","CD105","MCT1","MAG","MBP", 
            "PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","PHF1Tau", 
            "8OHGuano","pTDP43","PolyubiK48","PolyubiK63","MFN2","EEA1") #  "SERT","TH",     

##Read in csv
all_data_modMore <- read_csv(paste0(dataPath,"all_data_modMore.csv"))

#all_dataSubCluster = subset(all_data_full, Meta == "0") #Cluster ID to extract. for multiple Meta == c("endothelial", "microglia") what ever cluster ID
#all_dataSubCluster = subset(all_data_mod, MantisPopulation == "neurons")
all_dataSubCluster = subset(all_data_modMore, MantisPopulation == "neurons")

colSums(is.na(all_dataSubCluster)) # to show number of NAs...should to be 0
table(all_dataSubCluster$MantisPopulation)

# FOR neurons in HiAD Data
all_dataSubCluster = subset(all_dataSubCluster, !is.na(Calbindin)) # check all_dataSubCluster inconsole
SubclusteringMarkers = c("MAP2","VGLUT1","VGLUT2","VGAT","PanGAD6567", "Presenilin1NTF","Parvalbumin","Calbindin","Calretinin","MAG", "MBP") 

NumSubClusters = 30 # nClus = number of subcluster
#writeCytofFCS(as.matrix(all_dataSubCluster[, SubclusteringMarkers]), paste0(outputFolder, "/all_dataSubCluster.fcs"))

dataFlowframe = as_flowFrame(as.matrix(all_dataSubCluster[, SubclusteringMarkers]), source.frame = NULL)

fSOM <- FlowSOM(dataFlowframe, #(paste0(outputFolder, "/all_dataSubCluster.fcs"),
                compensate = F, transform = F,scale = T, 
                colsToUse = SubclusteringMarkers, nClus = NumSubClusters, seed = 140214)  

table(fSOM$metaclustering)
metaClustering = fSOM$metaclustering

fSOM_clustering = data.frame(fSOM$FlowSOM$map$mapping)
colnames(fSOM_clustering) = c("Cluster", "Value")
fSOM_clustering$Meta = as.numeric(fSOM$metaclustering[fSOM_clustering[,1]])

#to get weights for each marker
flowSOMCodes=fSOM[["FlowSOM"]][["map"]][["codes"]]
flowSOMCodes=data.frame(flowSOMCodes)
flowSOMCodes$Metacluster = fSOM$metaclustering 
flowSOMCodes$Metacluster = as.numeric(flowSOMCodes$Metacluster)

#_____________________________________________________________________________________________________________________________________________________________________________________
### HACK For Neurons, Inhibitory 1, glia 2, Excitatory 3, Undetermined 4
fSOM_clustering[which(fSOM_clustering$Meta == 1), "Meta"] = 1
fSOM_clustering[which(fSOM_clustering$Meta == 2), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 3), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 4), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 5), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 6), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 7), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 8), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 9), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 10), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 11), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 12), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 13), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 14), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 15), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 16), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 17), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 18), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 19), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 20), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 21), "Meta"] = 2
fSOM_clustering[which(fSOM_clustering$Meta == 22), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 23), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 24), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 25), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 26), "Meta"] = 1
fSOM_clustering[which(fSOM_clustering$Meta == 27), "Meta"] = 4
fSOM_clustering[which(fSOM_clustering$Meta == 28), "Meta"] = 3
fSOM_clustering[which(fSOM_clustering$Meta == 29), "Meta"] = 1
fSOM_clustering[which(fSOM_clustering$Meta == 30), "Meta"] = 3

# 
 NumSubClusters = 4
# 
flowSOMCodes[which(flowSOMCodes$Metacluster == 1), "Metacluster"] = 1
flowSOMCodes[which(flowSOMCodes$Metacluster == 2), "Metacluster"] = 2
flowSOMCodes[which(flowSOMCodes$Metacluster == 3), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 4), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 5), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 6), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 7), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 8), "Metacluster"] = 2
flowSOMCodes[which(flowSOMCodes$Metacluster == 9), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 10), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 11), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 12), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 13), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 14), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 15), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 16), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 17), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 18), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 19), "Metacluster"] = 2
flowSOMCodes[which(flowSOMCodes$Metacluster == 20), "Metacluster"] = 2
flowSOMCodes[which(flowSOMCodes$Metacluster == 21), "Metacluster"] = 2
flowSOMCodes[which(flowSOMCodes$Metacluster == 22), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 23), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 24), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 25), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 26), "Metacluster"] = 1
flowSOMCodes[which(flowSOMCodes$Metacluster == 27), "Metacluster"] = 4
flowSOMCodes[which(flowSOMCodes$Metacluster == 28), "Metacluster"] = 3
flowSOMCodes[which(flowSOMCodes$Metacluster == 29), "Metacluster"] = 1
flowSOMCodes[which(flowSOMCodes$Metacluster == 30), "Metacluster"] = 3

### END HACK HERE
#_____________________________________________________________________________________________________________________________________________________________________________________
### HACK For Disease 
# fSOM_clustering[which(fSOM_clustering$Meta == 1), "Meta"] = 1
# fSOM_clustering[which(fSOM_clustering$Meta == 2), "Meta"] = 2
# fSOM_clustering[which(fSOM_clustering$Meta == 3), "Meta"] = 3
# fSOM_clustering[which(fSOM_clustering$Meta == 4), "Meta"] = 4
# fSOM_clustering[which(fSOM_clustering$Meta == 5), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 6), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 7), "Meta"] = 6
# fSOM_clustering[which(fSOM_clustering$Meta == 8), "Meta"] = 7
# fSOM_clustering[which(fSOM_clustering$Meta == 9), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 10), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 11), "Meta"] = 8
# fSOM_clustering[which(fSOM_clustering$Meta == 12), "Meta"] = 9
# fSOM_clustering[which(fSOM_clustering$Meta == 13), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 14), "Meta"] = 10
# fSOM_clustering[which(fSOM_clustering$Meta == 15), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 16), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 17), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 18), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 19), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 20), "Meta"] = 11
# fSOM_clustering[which(fSOM_clustering$Meta == 21), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 22), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 23), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 24), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 25), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 26), "Meta"] = 12
# fSOM_clustering[which(fSOM_clustering$Meta == 27), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 28), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 29), "Meta"] = 5
# fSOM_clustering[which(fSOM_clustering$Meta == 30), "Meta"] = 5
# 
# 
# NumSubClusters = 12
# 
# flowSOMCodes[which(flowSOMCodes$Metacluster == 1), "Metacluster"] = 1 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 2), "Metacluster"] = 2 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 3), "Metacluster"] = 3 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 4), "Metacluster"] = 4 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 5), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 6), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 7), "Metacluster"] = 6
# flowSOMCodes[which(flowSOMCodes$Metacluster == 8), "Metacluster"] = 7
# flowSOMCodes[which(flowSOMCodes$Metacluster == 9), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 10), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 11), "Metacluster"] = 8 ###
# flowSOMCodes[which(flowSOMCodes$Metacluster == 12), "Metacluster"] = 9 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 13), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 14), "Metacluster"] = 10 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 15), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 16), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 17), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 18), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 19), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 20), "Metacluster"] = 11 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 21), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 22), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 23), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 24), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 25), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 26), "Metacluster"] = 12 #
# flowSOMCodes[which(flowSOMCodes$Metacluster == 27), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 28), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 29), "Metacluster"] = 5
# flowSOMCodes[which(flowSOMCodes$Metacluster == 30), "Metacluster"] = 5
### END HACK HERE
#_____________________________________________________________________________________________________________________________________________________________________________________
flowSOMmetaCodes = data.frame()
for (i in 1:NumSubClusters) {
  fscluster = subset(flowSOMCodes,Metacluster == i)
  fscluster = data.frame(t(apply(fscluster, 2, mean)))
  flowSOMmetaCodes = rbind(flowSOMmetaCodes,fscluster)
}
flowSOMmetaCodes$Metacluster = NULL
#_____________________________________________________________________________________________________________________________________________________________________________________

#_____________________________________________________________________________________________________________________________________________________________________________________
dev.off()
dev.off()

pdf(paste0(outputFolder, "/FlowSOMWeightsSubCluster.pdf"), width = 5, height = 5, onefile = F)
colfunc <- colorRampPalette(c("white", "black", "red"))
pheatmap(flowSOMmetaCodes, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "FlowSOM Weights,\n30 MetaClusters into 4 MetaClusters",
         col=colfunc(15),cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

pdf(paste0(outputFolder, "/FlowSOMWeightsColScaledSubCluster.pdf"), width = 5, height = 5, onefile = F)
colfunc <- colorRampPalette(c("white", "black", "red"))
pheatmap(flowSOMmetaCodes, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "FlowSOM Weights Column Scaled,\n30 MetaClusters into 4 MetaClusters",
         col=colfunc(15),cellwidth = 10, cellheight = 10)
dev.off()
dev.off()
#_____________________________________________________________________________________________________________________________________________________________________________________

#_____________________________________________________________________________________________________________________________________________________________________________________
#when using Flowsom
all_dataSubCluster$MetaSubCluster = fSOM_clustering$Meta 

all_dataSubCluster$MetaSubCluster = as.factor(all_dataSubCluster$MetaSubCluster)
#inhibitory Meta1 , oligodendrocytes Meta2 , excitatory Meta 3 , undetermined Meta4 
levels(all_dataSubCluster$MetaSubCluster) = c(paste0(c("inhibitory","oligodendrocytes","excitatory","undetermined"))) 

#when have  less than 100 cells and not runing flowSOM
#all_dataSubCluster$MetaSubCluster = 1:nrow(all_dataSubCluster)

fSOMcodes = fSOM$FlowSOM$map$codes
fSOMcodesMean = apply(fSOMcodes,2,mean)

table(all_dataSubCluster$MetaSubCluster)

write.csv(fSOMcodesMean, paste0(outputFolder,'/fSOMcodesMean_neuralSubCluster.csv'))
#_____________________________________________________________________________________________________________________________________________________________________________________

#_____________________________________________________________________________________________________________________________________________________________________________________
# Make a barplot of cluster sizes

table(all_dataSubCluster$MetaSubCluster)
table(fSOM_clustering$Meta)

#inhibitory Meta1 DeepPurple "#330066", nonneuronal Meta2 PurplePink "#990099", excitatory Meta 3 Pink "#CC0066", undetermined Meta4 Gray "#A5A5A5"
palette = c("#330066","#990099","#CC0066","#A5A5A5") 
#palette = distinctColorPalette(NumSubClusters)

plot_data = data.frame(table(all_dataSubCluster$MetaSubCluster))
ggplot(plot_data, aes(x = Var1, y = Freq)) +
  geom_bar(aes(fill = Var1), stat = "identity", colour = "black", fill = palette) +
  #geom_text(aes(x = Var1, y = (max(plot_data$Freq)/2), label = paste0("Cluster ", Var1, "\n", Freq)), col='black', size = 2.5) +
  xlab("FlowSOM metaclusters") + ylab("Cluster size") +
  theme_bw(15) +
  theme(legend.position="none")
ggsave(filename = paste0(outputFolder,"/FlowSOM_NeuralSubClusterSize.pdf"), width = 12, height = 3)

all_data_sampled <- data.frame()
unifSize = 10000000
for (i in 1:NumSubClusters) {
  all_data_subset = all_dataSubCluster
  all_data_subset = subset(all_data_subset, MetaSubCluster == i)
  
  if (unifSize < nrow(all_data_subset)) {
    unif = sample(1:nrow(all_data_subset), size = unifSize)
  } else {
    unif = 1:nrow(all_data_subset)
  }
  unif = all_data_subset[unif,]
  unif = melt(unif, id = c("MetaSubCluster"))
  colnames(unif) = c("MetaSubCluster", "Channel", "Values")
  
  all_data_sampled <<- rbind(all_data_sampled, unif)
}
###___________________________________________________________________________________________________________________________________________________

###___________________________________________________________________________________________________________________________________________________
## Make a heatmap for channels vs clusters with mean expression values from all_data. Scaled by marker (col) version of the above
plot_data = all_dataSubCluster[, c(markers, "MetaSubCluster")]

## to find out number of cell in each condition
sum(!is.na(plot_data$Amyloidbeta140)) # number of AD cells

## to remove NA columns to make AD specific HeatMap
# plot_data = plot_data[!is.na(plot_data$Amyloidbeta140),] 
table(plot_data$MetaSubCluster)

plot_data = melt(plot_data, id = c("MetaSubCluster"))
plot_data = dcast(plot_data, variable ~ MetaSubCluster, function(x){
  mean(x, na.rm = T )
})
plot_data$variable = NULL
plot_data = t(plot_data)

plot_data_rowNames = rownames(plot_data)

plot_data = apply(plot_data, 2, scale)#marker scaled
plot_data = rescale(plot_data, to = c(0, 1))#adjust between 0 and 2
colnames(plot_data) = markers
#rownames(plot_data) = 1:NumSubClusters
plot_data <- plot_data[,colSums(is.na(plot_data))<nrow(plot_data)]
rownames(plot_data) = plot_data_rowNames
rownames(plot_data)
plot_data = plot_data[c(1,3,2,4),] #Which should be #(inhibitory 1, nonneuronal 2, excitatory 3, undetermined 4))
rownames(plot_data)

SubclusteringMarkers = c("MAP2","VGLUT1","VGLUT2","VGAT","PanGAD6567", "Presenilin1NTF","Parvalbumin","Calbindin","Calretinin","MAG", "MBP")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% SubclusteringMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_SubclusteringMarkers0.pdf"), width = 7, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "30 FlowSOM MetaClusters into 4 MetaClusters, \nMean Marker Expression, \nGated Markers", 
         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)

#breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

AllMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","PanGAD6567","VGLUT1","VGLUT2","MAP2","8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","ApoE4","Calbindin","Calretinin","MFN2","PHF1Tau","PanAmyloidbeta1724","Parvalbumin","PolyubiK48","PolyubiK63","Presenilin1","pTDP43","Synaptophysin","PSD95")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% AllMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_AllMarkers0.pdf"), width = 7, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "30 FlowSOM MetaClusters into 4 MetaClusters, \nMean Marker Expression, \nAll All Markers", 
         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)

#breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

LineageMarkers = c("MBP","MAG","Iba1","CD45","GFAP","CD31","MCT1","CD105","VGAT","PanGAD6567","VGLUT1","VGLUT2","MAP2","Calbindin","Calretinin","Parvalbumin","Presenilin1NTF")
plot_data_subset = plot_data[,which(colnames(plot_data) %in% LineageMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_LineageMarkers0.pdf"), width = 7, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "30 FlowSOM MetaClusters into 4 MetaClusters, \nMean Marker Expression, \nLineage Markers", 
         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)

#breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

NeuralMarkers = c("VGAT","PanGAD6567","PanGAD6567","VGLUT1","VGLUT2","MAP2","Calbindin","Calretinin","Parvalbumin","Presenilin1NTF","MBP","MAG","Synaptophysin","PSD95") 
plot_data_subset = plot_data[,which(colnames(plot_data) %in% NeuralMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_NeuralMarkers0.pdf"), width = 7, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "30 FlowSOM MetaClusters into 4 MetaClusters, \nMean Marker Expression, \nNeural Markers", 
         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)

#breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()

ADSpecificMarkers = c("8OHGuano","PanAmyloidbeta1724","Amyloidbeta140","Amyloidbeta142","PHF1Tau","pTDP43","MFN2","PolyubiK48","PolyubiK63") 
plot_data_subset = plot_data[,which(colnames(plot_data) %in% ADSpecificMarkers)]
dev.off()
pdf(paste0(outputFolder, "/HeatMap_ADSpecificMarkers0.pdf"), width = 7, height = 7, onefile = F)
p=pheatmap(plot_data_subset, cluster_cols = T, cluster_rows = T, border_color = FALSE, scale = "none", main = "30 FlowSOM MetaClusters into 4 MetaClusters, \nMean Marker Expression, \nAD Specific Markers", 
         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)

#breakList = c(0.0,0.1,0.12,0.14,0.16,0.18,0.2,0.3,0.4,0.5,0.6,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)  #this one has 11color,  breakList=c(0,0,0,seq(0,1, by = 0.1))
#colfunc <- colorRampPalette(c("#60BEEE","#FFE5CC","#980808")) 
#color = viridis(length(plot_data_subset),  alpha=1, begin=0, end=1, direction=1, option = "D"), cellwidth = 10, cellheight = 10)
#col= colfunc(length(breakList)),breaks = breakList, cellwidth = 20, cellheight = 20)
#col= rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)[1:100]) , cellwidth = 10, cellheight = 10)

dev.off()
dev.off()
###___________________________________________________________________________________________________________________________________________________

###___________________________________________________________________________________________________________________________________________________

## DATA Prep for Violin plots
all_data_sampled <- data.frame()
unifSize = 10000000
for (i in 1:(NumSubClusters)) {
  all_data_subset = subset(all_dataSubCluster, MetaSubCluster == levels(all_dataSubCluster$MetaSubCluster)[i])

  if (unifSize < nrow(all_data_subset)) {
    unif = sample(1:nrow(all_data_subset), size = unifSize)
  } else {
    unif = 1:nrow(all_data_subset)
  }
  unif = all_data_subset[unif,]
  unif = melt(unif, id = c("MetaSubCluster"))
  colnames(unif) = c("MetaSubCluster", "Channel", "Values")

  all_data_sampled <<- rbind(all_data_sampled, unif)
}
###___________________________________________________________________________________________________________________________________________________

###___________________________________________________________________________________________________________________________________________________

## Violin plots per fSOM cluster for each channel ##

## inhibitory Meta1 DeepPurple "#330066", nonneuronal Meta2 PurplePink "#990099", excitatory Meta 3 Pink "#CC0066", undetermined Meta4 Gray "#A5A5A5"
#palette = c("#330066","#990099","#CC0066","#A5A5A5")
palette = c("#330066","#CC0066","#A5A5A5")

## Violin
plot_list = list()
for (channel in markers) {
  plot_data = subset(all_data_sampled, Channel == channel & MetaSubCluster %in% c("inhibitory","excitatory","undetermined"))
  plot_data$MetaSubCluster = as.factor(plot_data$MetaSubCluster)
  plot_data$MetaSubCluster = ordered(plot_data$MetaSubCluster, levels = c("inhibitory", "excitatory", "undetermined")) 
  g = ggplot(plot_data, aes(x = as.factor(MetaSubCluster), y = as.numeric(Values))) +
    geom_violin(aes(fill = MetaSubCluster), scale = "width", width = 0.6) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
    scale_fill_manual(values = palette) +
    xlab("") + ylab("") +
    ggtitle(channel) +
    ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
    theme_classic(12) +
    theme(legend.position = "none")
  plot_list[[channel]] <- g
}
picName = paste0(outputFolder, "/Violins_perChannel0.pdf")
cat(paste0("Plotting: ", picName, "\n"))
pdf(picName, width = 40, height = 60, onefile = T)
do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
dev.off()

## Scatter
plot_list = list()
for (channel in markers) {
  plot_data = subset(all_data_sampled, Channel == channel & MetaSubCluster %in% c("inhibitory","excitatory","undetermined"))
  plot_data$MetaSubCluster = as.factor(plot_data$MetaSubCluster)
  plot_data$MetaSubCluster = ordered(plot_data$MetaSubCluster, levels = c("inhibitory", "excitatory", "undetermined")) 
  g = ggplot(plot_data, aes(x = as.factor(MetaSubCluster), y = as.numeric(Values))) +
  geom_point(aes(color = MetaSubCluster), cex = 1, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.9)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white", notch = F, outlier.size = 0) +
  scale_color_manual(values = palette) +
    xlab("") + ylab("") +
    ggtitle(channel) +
    ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
    theme_classic(12) +
    theme(legend.position = "none")
  plot_list[[channel]] <- g
}
picName = paste0(outputFolder, "/Scatter_perChannel0.pdf")
cat(paste0("Plotting: ", picName, "\n"))
pdf(picName, width = 40, height = 60, onefile = T)
do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
dev.off()

###_______________________________________________________________________________________________________________________________________#
 
###_______________________________________________________________________________________________________________________________________________

## Plots for Injury ##
plot_data = subset(all_data_sampled, MetaSubCluster %in% c("inhibitory","excitatory","undetermined")) # "inhibitory","excitatory","undetermined" 
plot_data$MetaSubCluster = as.factor(plot_data$MetaSubCluster)
plot_data$MetaSubCluster = ordered(plot_data$MetaSubCluster, levels = c("inhibitory", "excitatory","undetermined")) 
levels(plot_data$MetaSubCluster)

## inhibitory Meta1 DeepPurple "#330066", nonneuronal Meta2 PurplePink "#990099", excitatory Meta 3 Pink "#CC0066", undetermined Meta4 Gray "#A5A5A5"
palette = c("#330066","#CC0066","#A5A5A5")  #,"#A5A5A5"

## ViolinPlot
violinMarkers = c("8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63") 
violin_data = subset(plot_data, Channel %in% violinMarkers)
violin_data$Channel = ordered(violin_data$Channel, levels = violinMarkers)

g = ggplot(violin_data, aes(x = Channel, y = as.numeric(Values), fill = MetaSubCluster)) +
  geom_violin(scale = "width", width = 0.8) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.size = 0) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(0,1) +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Violins_Injury2",".pdf"),g, width = 8, height = 6)

## ScatterPlot
scatterMarkers = c("8OHGuano","pTDP43","MFN2","PolyubiK48","PolyubiK63")
scatter_data = subset(plot_data, Channel %in% scatterMarkers)
scatter_data$Channel = ordered(scatter_data$Channel, levels = violinMarkers)
scatter_data$Values = as.numeric(scatter_data$Values)

## Number of cells and Mean expression
scatter_tablefull = data.frame(table(scatter_data$MetaSubCluster, scatter_data$Channel))
scatter_tablefull$mean = aggregate(scatter_data, by = list(scatter_data$MetaSubCluster, scatter_data$Channel),FUN =  mean)$Values
write.csv(scatter_tablefull, paste0(outputFolder, "/scatter_tablefull_Injury.csv"))

scatter_data = subset(scatter_data, Values>0.001 & Values<1.5)

scatter_tablesubset = data.frame(table(scatter_data$MetaSubCluster, scatter_data$Channel))
scatter_tablesubset$mean = aggregate(scatter_data, by = list(scatter_data$MetaSubCluster, scatter_data$Channel),FUN =  mean)$Values
write.csv(scatter_tablesubset, paste0(outputFolder, "/scatter_tablesubset_Injury.csv"))

g = ggplot(scatter_data, aes(x = Channel, y = Values, color = MetaSubCluster)) +
  geom_point(cex = 0.8, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9),injury
               outlier.size = 0, fill = "white", notch = T, outlier.alpha = 0, 
               ymin = quantile(scatter_data$Values, 0.6), ymax = quantile(scatter_data$Values, 0.6)) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(c(0.001,1.5)) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Scatter_Injury1",".pdf"),g, width = 8, height = 8)

## Shapiro-Wilk normality test for the differences and and NonParametric Test
statTestResults = data.frame(matrix(NA,length(scatterMarkers), 4))
colnames(statTestResults) = c("Normality","inhibitory_excitatory", "excitatory_undetermined", "inhibitory_undetermined")
rownames(statTestResults) = scatterMarkers
for (marker in scatterMarkers){
  marker_data = subset(scatter_data, Channel == marker)
  marker_data$MetaSubCluster = as.character(marker_data$MetaSubCluster)
  res = shapiro.test(marker_data[,"Values"])
  statTestResults [marker,1] = res$p.value
  
  # Compute t-test for inhibitory_excitatory
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "excitatory"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,2] = res$p.value
  
  # Compute t-test for excitatory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("excitatory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,3] = res$p.value
  
  # Compute t-test for inhibitory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,4] = res$p.value
}
write.csv(statTestResults, paste0(outputFolder,'/NonParametricTestInjury.csv'))

###_______________________________________________________________________________________________________________________________________________

###_______________________________________________________________________________________________________________________________________________

## Plots for Amyloid ##
plot_data = subset(all_data_sampled, MetaSubCluster %in% c("inhibitory","excitatory","undetermined")) # "inhibitory","excitatory","undetermined"
plot_data$MetaSubCluster = as.factor(plot_data$MetaSubCluster)
plot_data$MetaSubCluster = ordered(plot_data$MetaSubCluster, levels = c("inhibitory", "excitatory","undetermined")) # levels(plot_data$MetaSubCluster) in console to check if
levels(plot_data$MetaSubCluster)
## inhibitory Meta1 DeepPurple "#330066", glia Meta2 PurplePink "#990099", excitatory Meta 3 Pink "#CC0066", undetermined Meta4 Gray "#A5A5A5"
palette = c("#330066","#CC0066","#A5A5A5") # ,"#A5A5A5"

## ViolinPlot
violinMarkers = c("PanAmyloidbeta1724","Amyloidbeta140", "Amyloidbeta142")
violin_data = subset(plot_data, Channel %in% violinMarkers)

g = ggplot(violin_data, aes(x = Channel, y = as.numeric(Values), fill = MetaSubCluster)) +
  geom_violin(scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.62), outlier.size = 0) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(0,3) +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Violins_Amyloid1",".pdf"),g, width = 8, height = 8)

## ScatterPlot 
scatterMarkers = c("PanAmyloidbeta1724","Amyloidbeta140", "Amyloidbeta142") 
scatter_data = subset(plot_data, Channel %in% scatterMarkers)
scatter_data$Channel = ordered(scatter_data$Channel, levels = violinMarkers)
scatter_data$Values = as.numeric(scatter_data$Values)

## Number of cells and Mean expression
scatter_tablefull = data.frame(table(scatter_data$MetaSubCluster, scatter_data$Channel))
scatter_tablefull$mean = aggregate(scatter_data, by = list(scatter_data$MetaSubCluster, scatter_data$Channel),FUN =  mean)$Values
write.csv(scatter_tablefull, paste0(outputFolder, "/scatter_tablefull_Amyloid.csv"))

scatter_data = subset(scatter_data, Values>0 & Values<5)

scatter_tablesubset = data.frame(table(scatter_data$MetaSubCluster, scatter_data$Channel))
scatter_tablesubset$mean = aggregate(scatter_data, by = list(scatter_data$MetaSubCluster, scatter_data$Channel),FUN =  mean)$Values
write.csv(scatter_tablesubset, paste0(outputFolder, "/scatter_tablesubset_Amyloid.csv"))

g = ggplot(scatter_data, aes(x = Channel, y = Values, color = MetaSubCluster)) +
  geom_point(cex = 0.8, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9)) + 
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9),
               outlier.size = 0, fill = "white", notch = T, outlier.alpha = 0, 
               ymin = quantile(scatter_data$Values, 0.6), ymax = quantile(scatter_data$Values, 0.6)) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(c(0,5)) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Scatter_Amyloid1",".pdf"),g, width = 8, height = 8)

## Shapiro-Wilk normality test for the differences and and NonParametric Test
statTestResults = data.frame(matrix(NA,length(scatterMarkers), 4))
colnames(statTestResults) = c("Normality","inhibitory_excitatory", "excitatory_undetermined", "inhibitory_undetermined")
rownames(statTestResults) = scatterMarkers
for (marker in scatterMarkers){
  marker_data = subset(scatter_data, Channel == marker)
  marker_data$MetaSubCluster = as.character(marker_data$MetaSubCluster)
  res = shapiro.test(marker_data[,"Values"])
  statTestResults [marker,1] = res$p.value
  
  # Compute t-test for inhibitory_excitatory
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "excitatory"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,2] = res$p.value
  
  # Compute t-test for excitatory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("excitatory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,3] = res$p.value
  
  # Compute t-test for inhibitory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,4] = res$p.value
}
write.csv(statTestResults, paste0(outputFolder,'/NonParametricTestAmyloid.csv'))
###_______________________________________________________________________________________________________________________________________________ 

###_______________________________________________________________________________________________________________________________________________  
## Plots for NFT ##
plot_data = subset(all_data_sampled, MetaSubCluster %in% c("inhibitory","excitatory","undetermined")) # "inhibitory","excitatory","undetermined"
plot_data$MetaSubCluster = as.factor(plot_data$MetaSubCluster)
plot_data$MetaSubCluster = ordered(plot_data$MetaSubCluster, levels = c("inhibitory", "excitatory","undetermined")) 
levels(plot_data$MetaSubCluster)
## inhibitory Meta1 DeepPurple "#330066", glia Meta2 PurplePink "#990099", excitatory Meta 3 Pink "#CC0066", undetermined Meta4 Gray "#A5A5A5"
palette = c("#330066","#CC0066","#A5A5A5") #,"#A5A5A5"

## ViolinPlot 
violinMarkers = c("PHF1Tau")
violin_data = subset(plot_data, Channel %in% violinMarkers)

g = ggplot(violin_data, aes(x = Channel, y = as.numeric(Values), fill = MetaSubCluster)) +
  geom_violin(scale = "width", width = 0.6) +
  geom_boxplot(width = 0.07, position = position_dodge(width = 0.6), outlier.size = 0, notch = T,) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(0,5) +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Violins_NFT1",".pdf"),g, width = 8, height = 8)

## ScatterPlot
scatterMarkers = c("PHF1Tau") 
scatter_data = subset(plot_data, Channel %in% scatterMarkers)
scatter_data$Channel = ordered(scatter_data$Channel, levels = violinMarkers)
scatter_data$Values = as.numeric(scatter_data$Values)

## Number of cells and Mean expression
scatter_tablefull = data.frame(table(scatter_data$MetaSubCluster, scatter_data$Channel))
scatter_tablefull$mean = aggregate(scatter_data, by = list(scatter_data$MetaSubCluster, scatter_data$Channel),FUN =  mean)$Values
write.csv(scatter_tablefull, paste0(outputFolder, "/scatter_tablefull_NFT.csv"))

scatter_data = subset(scatter_data, Values>0 & Values<5)

scatter_tablesubset = data.frame(table(scatter_data$MetaSubCluster, scatter_data$Channel))
scatter_tablesubset$mean = aggregate(scatter_data, by = list(scatter_data$MetaSubCluster, scatter_data$Channel),FUN =  mean)$Values
write.csv(scatter_tablesubset, paste0(outputFolder, "/scatter_tablesubset_NFT.csv"))

g = ggplot(scatter_data, aes(x = Channel, y = Values, color = MetaSubCluster)) +
  geom_point(cex = 2, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.4)) + 
  geom_boxplot(width = 0.04, position = position_dodge(width = 0.4),
               outlier.size = 0, fill = "white", notch = T, outlier.alpha = 0, 
               ymin = quantile(scatter_data$Values, 0.6), ymax = quantile(scatter_data$Values, 0.6)) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(c(0,5)) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Scatter_NFT1",".pdf"),g, width = 8, height = 8)

## Shapiro-Wilk normality test for the differences and and NonParametric Test
statTestResults = data.frame(matrix(NA,length(scatterMarkers), 4))
colnames(statTestResults) = c("Normality","inhibitory_excitatory", "excitatory_undetermined", "inhibitory_undetermined")
rownames(statTestResults) = scatterMarkers
for (marker in scatterMarkers){
  marker_data = subset(scatter_data, Channel == marker)
  marker_data$MetaSubCluster = as.character(marker_data$MetaSubCluster)
  #res = shapiro.test(marker_data[,"Values"])
  #statTestResults [marker,1] = res$p.value
  
  # Compute t-test for inhibitory_excitatory
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "excitatory"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,2] = res$p.value
  
  # Compute t-test for excitatory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("excitatory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,3] = res$p.value
  
  # Compute t-test for inhibitory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,4] = res$p.value
}
write.csv(statTestResults, paste0(outputFolder,'/NonParametricTest.csv'))  
###_____________________________________________________________________________________

  
  
  
  



###_____________________________________________________________________________________

##### STATISTICAL TEST ###########

###_____________________________________________________________________________________
## Shapiro-Wilk normality test for the differences and and NonParametric Test
statTestResults = data.frame(matrix(NA,length(scatterMarkers), 4))
colnames(statTestResults) = c("Normality","inhibitory_excitatory", "excitatory_undetermined", "inhibitory_undetermined")
rownames(statTestResults) = scatterMarkers
for (marker in scatterMarkers){
  marker_data = subset(scatter_data, Channel == marker)
  marker_data$MetaSubCluster = as.character(marker_data$MetaSubCluster)
  res = shapiro.test(marker_data[,"Values"])
  statTestResults [marker,1] = res$p.value

  # Compute t-test for inhibitory_excitatory
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "excitatory"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,2] = res$p.value
  
  # Compute t-test for excitatory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("excitatory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,3] = res$p.value
  
  # Compute t-test for inhibitory_undetermined
  test_data = subset(marker_data, MetaSubCluster %in% c("inhibitory", "undetermined"))
  res = wilcox.test(as.formula("Values~MetaSubCluster"), test_data, paired = F)
  statTestResults [marker,4] = res$p.value
}
write.csv(statTestResults, paste0(outputFolder,'/NonParametricTestInjury.csv'))
###_____________________________________________________________________________________



##### MORE PLOTTING OPTIONS###########

###_______________________________________________________________________________________________________________________________________#
# ## Violin plots Injury Markers
# plot_list = list()
# violinMarkers = c("8OHGuano","pTDP43","ApoE4","MFN2","PolyubiK48","PolyubiK63")
# #violinMarkers = markers #  = markers ## for all markers
# for (meta in levels(all_dataSubCluster$MetaSubCluster)) {
#   #meta = c("excitatory") #run if single population # inhibitory, glia, excitatory, undetermined 
#   plot_data = subset(all_data_sampled, MetaSubCluster == meta)
#   plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
#   plot_data = subset(plot_data, Channel %in% violinMarkers)
#   g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.05, position = position_dodge(width = 0.9), outlier.size = -10) +
#     #scale_fill_manual(values = colorRampPalette(c("#000000"))(10)[1:10]) + # colorRampPalette(c("deepskyblue4", "darkorange")), colorRampPalette(brewer.pal(9,"Blue"))
#     xlab("") + ylab("") +
#     ylim(0,1.0) +
#     ggtitle(paste0("Cluster: ", meta)) +
#     theme_classic(12) +
#     theme(legend.position = "none") +
#     theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.x = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.y = ggplot2::element_text(face = "bold", size = 10),
#           axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
#   
#   # for single plot
#   #ggsave(paste0(outputFolder, "/Violins_", meta,".pdf"),g, width = 6, height = 6)
#   
#   plot_list[[meta]] <- g
# }
# picName = paste0(outputFolder, "/Violins_InjuryMarks.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 15, height = 5, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 2))
# dev.off()
# 
# 
# ## Plot scatter + boxplots for each channel separately. Show all clusters on each plot
# plot_list = list()
# violinMarkers = c("8OHGuano","pTDP43","ApoE4","MFN2","PolyubiK48","PolyubiK63")
# for (meta in levels(all_dataSubCluster$MetaSubCluster)) {
#   #meta = c("excitatory") #run if single population # inhibitory, glia, excitatory, undetermined 
#   plot_data = subset(all_data_sampled, MetaSubCluster == meta)
#   plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
#   plot_data = subset(plot_data, Channel %in% violinMarkers)
#   g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
#     geom_point(aes(color = Channel), cex = 0.5, alpha = 0.8, position = position_jitter()) +
#     geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
#     scale_color_manual(values = palette) +
#     xlab("") + ylab("") +
#     ggtitle(paste0("Cluster: ", meta)) +
#     ylim(c(0, 1)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic(12) +
#     theme(legend.position = "none")
#   plot_list[[meta]] <- g
# }
# picName = paste0(outputFolder, "/Scatter_InjuryMarks.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 15, height = 5, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 2))
# dev.off()
# 
# ###_______________________________________________________________________________________________________________________________________#
# ## Violin plots Amyloid
# plot_list = list()
# violinMarkers = c("PanAmyloidbeta1724","Amyloidbeta140", "Amyloidbeta142")
# 
# ## Plot violin + boxplots for each channel separately. Show all clusters on each plot
# palette = c("#A0A0A0","#A0A0A0","#A0A0A0")
# for (meta in levels(all_dataSubCluster$MetaSubCluster)) {
#   meta = c("undetermined") #run if single population # inhibitory, glia, excitatory, undetermined 
#   plot_data = subset(all_data_sampled, MetaSubCluster == meta)
#   plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
#   plot_data = subset(plot_data, Channel %in% violinMarkers)
#   g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.05, position = position_dodge(width = 0.9), outlier.size = 0) +
#     scale_fill_manual(values = colorRampPalette(c("#000000"))(10)[1:10]) + # colorRampPalette(c("deepskyblue4", "darkorange")), colorRampPalette(brewer.pal(9,"Blue"))
#     xlab("") + ylab("") +
#     ylim(0.0000001,0.2) +
#     ggtitle(paste0("Cluster: ", meta)) +
#     theme_classic(12) +
#     theme(legend.position = "none") +
#     theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.x = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.y = ggplot2::element_text(face = "bold", size = 10),
#           axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
#   
#   # for single plot
#   #ggsave(paste0(outputFolder, "/Violins_", meta,".pdf"),g, width = 6, height = 6)
#   
#   plot_list[[meta]] <- g
# }
# picName = paste0(outputFolder, "/Violins_DiseaseMarksUnderm.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 5, height = 5, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 1))
# dev.off()
# 
# 
# ## Plot scatter + boxplots for each channel separately. Show all clusters on each plot
# palette = c("#000000","#000000","#000000")
# for (meta in levels(all_dataSubCluster$MetaSubCluster)) {
#   meta = c("excitatory") #run if single population # inhibitory, glia, excitatory, undetermined 
#   plot_data = subset(all_data_sampled, MetaSubCluster == meta)
#   plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
#   plot_data = subset(plot_data, Channel %in% violinMarkers)
#   g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
#     geom_point(aes(color = Channel), cex = 3, alpha = 0.3, position = position_jitter()) +
#     geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
#     scale_color_manual(values = palette) +
#     xlab("") + ylab("") +
#     ggtitle(paste0("Cluster: ", meta)) +
#     ylim(c(-0.001,1)) + # watch for warnings to see how many outliers were removed in each plot.
#     theme_classic(12) +
#     theme(legend.position = "none")
#   plot_list[[meta]] <- g
# }
# picName = paste0(outputFolder, "/Scatter_DiseaseMarksExcitInhib.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 15, height = 5, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 2))
# dev.off()
# 
# 
# 
# ###_______________________________________________________________________________________________________________________________________#
# ## Violin plots Tau
# plot_list = list()
# violinMarkers = c("PHF1Tau")
# #violinMarkers = markers #  = markers ## for all markers
# 
# for (meta in levels(all_dataSubCluster$MetaSubCluster)) {
#   meta = c("excitatory") #run if single population # inhibitory, glia, excitatory, undetermined 
#   plot_data = subset(all_data_sampled, MetaSubCluster == meta)
#   plot_data$Channel = ordered(plot_data$Channel, levels = violinMarkers)
#   plot_data = subset(plot_data, Channel %in% violinMarkers)
#   g = ggplot(plot_data, aes(x = Channel, y = as.numeric(Values))) +
#     geom_violin(aes(fill = Channel), scale = "width", width = 0.6) +
#     geom_boxplot(width = 0.05, position = position_dodge(width = 0.9), outlier.size = 0) +
#     #scale_fill_manual(values = colorRampPalette(c("#000000"))(10)[1:10]) + # colorRampPalette(c("deepskyblue4", "darkorange")), colorRampPalette(brewer.pal(9,"Blue"))
#     xlab("") + ylab("") +
#     ylim(0.0000001,0.05) +
#     ggtitle(paste0("Cluster: ", meta)) +
#     theme_classic(12) +
#     theme(legend.position = "none") +
#     theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.x = ggplot2::element_text(face = "bold", size = 10),
#           axis.title.y = ggplot2::element_text(face = "bold", size = 10),
#           axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
#   
#   # for single plot
#   #ggsave(paste0(outputFolder, "/Violins_", meta,".pdf"),g, width = 6, height = 6)
#   
#   plot_list[[meta]] <- g
# }
# ## Plot violin + boxplots for each channel separately. Show all clusters on each plot
# picName = paste0(outputFolder, "/Violins_DiseaseMarksTauExcitInhib.pdf")
# cat(paste0("Plotting: ", picName, "\n"))
# pdf(picName, width = 5, height = 5, onefile = T)
# do.call(grid.arrange, list(grobs = plot_list, ncol = 3))
# dev.off()








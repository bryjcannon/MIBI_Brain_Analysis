## FOLD CHANGE NEIGHBOURHOOD ANALYSIS ON GRADED (HI MED LO) OF TANGLES AND THREADS############

### LIBRARIES ###
libs <- c('ggplot2','RColorBrewer','reshape2','devtools','FlowSOM','Rtsne',
          'pheatmap','flowCore','plyr','scales','coin','grid',
          'gridExtra','tidyr','randomcoloR','viridis','tidyverse','ggpubr')
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
setwd("/Volumes/KausaliaHD 1/R4R")
source("/Volumes/KausaliaHD 1/R4R/General/readWriteCytofFCS.R")

## loads csv file of Real data that is generated from MIBIcreateObjectNeighborMatrixEzSeg.m
dataPath = "/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/Pooled/" 

## load csv file of graded/gated PHF1Tau objects
tauGrades <- read_csv("/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/CellEngineCurated/EzDeepcellPooled/alldataDiseasePooled_TotalTauvsPHF1_GatedPopulations/PHF1Tau_Grades/TangleThread_TotalTauvsPHF1Tau/TotalTauvsPHF1TauGrades.csv")

## loads Tangle Neighbourhoods
plot_data <- read_csv(paste0(dataPath,"RealRuns/TangleThreads_MatrixByIndividual_Objects/tauTangles_TangleThread_counts_50px.csv" )) # _freqs_ or _counts_
outputFolder = paste0(dataPath,"DiseaseNeighbourhoodEZDeepcell/FoldChangeObjects/TangleNeighbourhood_50pixeldistance(11May21)/")
dir.create(outputFolder, showWarnings= F, recursive = T)

#plot_data = ceiling(plot_data)
plot_data_realunmelted = plot_data

plot_data_real = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_real) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_real$MantisPopulation)

## loads csv file of Random data that is generated from MIBIcreateObjectNeighborMatrixEzSegRandom.m
nrun = 1
fulldata = read_csv(paste0(dataPath,"RandomRuns/TangleThread_MatrixByIndividual_Objects(11May21)/tauTangles_TangleThread_counts_50px_1.csv")) # _freqs_ or _counts_  
fulldata = as.matrix(fulldata)
# for (runnum in 2:nrun){
#   ## loads Amyloid Neighbourhoods
#   plot_data = read_csv(paste0(dataPath,"RandomRuns/TangleThread_MatrixByIndividual_Objects/tauTangles_TangleThread_counts_50px_",runnum,".csv"))  
#   fulldata[,3:6] = fulldata[,3:6]+as.matrix(plot_data[,3:6])
# }
# fulldata[,3:6] = fulldata[,3:6]/nrun
fulldata[fulldata==0] = 1
#fulldata = ceiling(fulldata)
plot_data = as.data.frame(fulldata)

plot_data_randunmelted = plot_data
plot_data_rand = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_rand) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_rand$MantisPopulation)

#make violin for fold change
plot_data = plot_data_real
plot_data$Value = plot_data$Value/plot_data_rand$Value
table(plot_data$MantisPopulation)
#plot_data$MantisPopulation = ordered(plot_data$MantisPopulation, levels = c("neurons", "astrocyte_process",  "vessel_CD31_CD105", "microglia_process" ))  # in console to check if levels(plot_data$Meta)
plot_data$MantisPopulation = ordered(plot_data$MantisPopulation, levels = c("neurons", "astrocyte_process", "microglia_process" ))  # in console to check if levels(plot_data$Meta)
levels(plot_data$MantisPopulation)

#plot_data = subset(plot_data, Object_id<3) 
## astrocytes Green "#106900", endothelial Blue"#0028D4",  microglia mango "#FFBF00",neurons DarkRed "#980A0A"...Color order follows order BEFORE custome ordering
#palette = c("#980A0A","#106900","#0028D4","#FFBF00")
palette = c("#980A0A","#106900","#FFBF00")

#Group ViolinPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = log2(Value))) +   # y =  MaxDist, MeanDist, MedianDist
  geom_violin(aes(fill = MantisPopulation), scale = "width", width = 0.6, adjust = 5) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(-7,7) +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Violin_Tangles_FC_50Pixeldistance.pdf"), plot = p, width = 5, height = 5)

#Group ScatterPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = log2(Value))) +   # y =  MaxDist, MeanDist, MedianDist
  geom_point(aes(color = MantisPopulation), cex = 0.5, alpha = 0.3, position = position_jitter(height = 0)) + # for freq no vertical jitter set to 0, for counts set 0.2
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(c(-7,7)) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(15) +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Scatter_Tangles_FC_50Pixeldistance.pdf"), plot = p, width = 7, height = 5)

#Boxplot with anova
#myComparison = list(c("astrocyte_process", "microglia_process"), c("astrocyte_process", "neurons"), c("neurons", "microglia_process"), c("neurons", "vessel_CD31_CD105"), c("microglia_process", "vessel_CD31_CD105"), c("astrocyte_process", "vessel_CD31_CD105"))
myComparison = list(c("astrocyte_process", "microglia_process"), c("astrocyte_process", "neurons"), c("neurons", "microglia_process"))
plot_data$log2Value = log2(plot_data$Value)
plot_data = subset(plot_data, Value > 0)
res = aov(Value ~ MantisPopulation, plot_data)
TukeyHSD(res)
p = ggboxplot(plot_data,x = "MantisPopulation",y = "Value") +
  stat_compare_means(comparisons=myComparison, method = "anova", aes(label=..p.ad..)) 
ggsave(paste0(outputFolder, "/Tangle_Boxplot_50pxDist.pdf"), plot = p, width = 7, height = 5)

#tally of object 
table(plot_data$MantisPopulation)

for (myComp in myComparison) {
  print(myComp)
  testdata = subset(plot_data, MantisPopulation %in% myComp)
  testdata$MantisPopulation = as.numeric(testdata$MantisPopulation)
  results = aov(MantisPopulation~log2Value, testdata)
  print(summary(results))
}

#Table with T-test and means
compare_means(Value~MantisPopulation,plot_data, method = 'anova')
library(dplyr)
plot_data %>%
  group_by(MantisPopulation) %>%
  summarise(mean = mean(Value), sd = sd(Value))

##Heatmap to plot Tangle Neighbourhood Flavours
set.seed(240214)
callback=function(x,data){
  hclust_dist = dist(t(as.matrix(plot_data)), method = "euclidean")
  hclust(hclust_dist, "ward.D2")
}

plot_data =  as.matrix(plot_data_realunmelted)
colnames_plotdata = paste0("p",plot_data[,1], "_objID", plot_data[,2])
rownames_plotdata = c(colnames(plot_data_realunmelted)[3:6],"CL")

## integrate PHFtau grade/gated cell engine population
tempdata = as.data.frame(plot_data)
tempdata$tauGrade = NA
  for (i in 1:nrow(tempdata)){
   tauObject =plot_data[i,]
   tauGrade = subset(tauGrades, point_id==tauObject[["Point_num"]] & obj_id==tauObject[["Object_num"]]) 
  if (nrow(tauGrade)>0){
   tempdata[i,"tauGrade"]=tauGrade[1,"file"]  
  }
}

table(tempdata$tauGrade)

tempdata$tauGrade = as.numeric(as.factor(tempdata$tauGrade))

plot_data[,3:6] = log2(plot_data[,3:6]/as.matrix(plot_data_randunmelted[,3:6]))

#To replace neg infinity or na with minimum value for each row
for (i in 3:ncol(plot_data)){
  plot_data[which(is.infinite(plot_data[,i])),i] = min(plot_data[which(is.finite(plot_data[,i])),i])
}

#To cluster for the overlay baseed on heatmap (n=10) by FLOWSOM
## if running FlowSOM
fcsExprs <- flowCore::flowFrame(plot_data)
plot_data = as.data.frame(plot_data)

seed = 140214

numClusters = 10 # dont change cos original flowSom cluster

clusteringLineage = colnames(plot_data[,3:6])

fSOM <- FlowSOM(fcsExprs,
                compensate = F, transform = F,scale = T,
                colsToUse = clusteringLineage, nClus = numClusters, seed = seed)

#flowSOM MetaCluster size
table(fSOM$metaclustering)
metaClustering = fSOM$metaclustering

fSOM_clustering = data.frame(fSOM$FlowSOM$map$mapping)
colnames(fSOM_clustering) = c("Cluster", "Value")
fSOM_clustering$Meta = as.numeric(fSOM$metaclustering[fSOM_clustering[,1]])

plot_data$ClusteringLineage = fSOM_clustering$Meta

write.csv(plot_data, paste0('/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/EzSegData_uci2712J/ezSegDiseaseGrades/ClusteringLineageTangles.csv'))
#########################################

plot_data =  t(plot_data)
plot_data = plot_data[3:7,]
plot_data = as.data.frame(plot_data)
#plot_data = subset(plot_data, Object_num<3) 
colnames(plot_data) =  1:ncol(plot_data)
rownames(plot_data) = rownames_plotdata

plot_data[nrow(plot_data) + 1,] = tempdata$tauGrade
rownames(plot_data) = c(rownames_plotdata, "tauGrade")

plot_data = plot_data[,!is.na(plot_data[nrow(plot_data),])]

anno_data = plot_data["tauGrade", ]
anno_data = as.data.frame(t(anno_data))
anno_data$tauGrade = as.factor(anno_data$tauGrade)
levels(anno_data$tauGrade) = c("hi","lo","med")
anno_data$tauGrade = ordered(anno_data$tauGrade, levels = c("hi","med","lo"))
levels(anno_data$tauGrade)
anno_color = list(tauGrade= c(hi = "#a80077", med = "#f4dd66", lo ="#67872f"))

#tally of hi med lo
table(anno_data$tauGrade)

# to add CL annotation
cl_data = plot_data["CL", ]
cl_data = as.data.frame(t(cl_data))
cl_data$CL = as.factor(cl_data$CL)
levels(cl_data$CL) = c("a","b","c","d","e","f","g","h","i","j")
#palette <- c("#ff0000", "#F58DA2", "#675545", "#000426", "#00AED9", "#AD4B56", "#6F2789", "#3D5567", "#53DA3F","#c7b299") #1:10
cl_color = list(CL= c(a="#ff0000", b="#F58DA2", c="#675545", d="#000426", e="#00AED9", f="#AD4B56", g="#6F2789", h="#3D5567", i="#53DA3F",j="#c7b299"))

#to combine annotations
anno_data = cbind(anno_data, cl_data)
anno_color = list(tauGrade= c(hi = "#a80077", med = "#f4dd66", lo ="#67872f"),CL= c(a="#ff0000",  b="#F58DA2", c="#675545", d="#000426", e="#00AED9", f="#AD4B56", g="#6F2789", h="#3D5567", i="#53DA3F",j="#c7b299") )

plot_data =  plot_data[, order(plot_data["CL",])]
plot_data = plot_data[c(1,2,4,3),]

#Barplot for tally numbers of gated objects
plot_datacounts = plot_data[1:nrow(plot_data),0]
plot_datacounts$counts = 0
for (i in 1:nrow(plot_data)){
  plot_datacounts[i,1] = ncol(plot_data)-length(which(plot_data[i,]==min(plot_data[i,])))
}
plot_datacounts=melt(plot_datacounts)
plot_datacounts$variable = factor(rownames(plot_data), levels = rownames(plot_data))
ggplot(plot_datacounts, aes(x=variable, y=value)) +
  geom_bar(stat = "identity")+
  theme_classic(20) +
  #ylim(0,4000)
ggsave(paste0(outputFolder,"/Bargraph_tanglethreadnumberspercelltype.pdf"),  width = 12, height = 10)

##NO SCALING##
# this plots Rocket color palette heatmaps is with tauGrade Clustering
pdf(paste0(outputFolder,"/Heatmap_TangleGrade_FC_50Pixeldistance1.pdf"), width = 20, height = 15, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = F, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "none", main = "Log2 Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nPixelDistance = 50",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# these heatmaps is without tauGrade Clustering
#plot_data = plot_data[1:(nrow(plot_data) - 1),] 

# this plots heatmaps Rocket color palette with lineage Clustering
pdf(paste0(outputFolder,"/Heatmap2_Tangles_FC_50Pixeldistance1.pdf"), width = 20, height = 15, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "none", main = "Log2 Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nPixelDistance = 50",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# this plots heatmaps Viridis color palette with lineage Clustering
pdf(paste0(outputFolder,"/Heatmap1_Tangles_FC_50Pixeldistance1.pdf"), width = 20, height = 15, onefile = FALSE)
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "none", main = "Log2 Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nPixelDistance = 50",
            color = viridis(100, alpha=1, begin=0, end=1, direction=1, option = "D"), cellheight = 60)
dev.off()
dev.off()

##ROW SCALING##
# this plots Rocket color palette heatmaps is with tauGrade Clustering
pdf(paste0(outputFolder,"/Heatmap_TangleGrade_FC_50PixeldistanceScaled2.pdf"), width = 20, height = 15, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = F,  annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "row", main = "Log2 Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nRowScaled nPixelDistance = 50",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# these heatmaps is without tauGrade Clustering
#plot_data = plot_data[1:(nrow(plot_data) - 1),] 

# this plots heatmaps Rocket color palette with lineage Clustering
pdf(paste0(outputFolder,"/Heatmap2_Tangles_FC_50PixeldistanceScaled2.pdf"), width = 20, height = 15, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = F,  annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "row", main = "Log2 Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nRowScaled nPixelDistance = 50",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# this plots heatmaps Viridis color palette with lineage Clustering
pdf(paste0(outputFolder,"/Heatmap1_Tangles_FC_50PixeldistanceScale2.pdf"), width = 20, height = 15, onefile = FALSE)
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "row", main = "Log2 Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nRowScaled nPixelDistance = 50",
            color = viridis(100, alpha=1, begin=0, end=1, direction=1, option = "D"), cellheight = 60)
dev.off()
dev.off()


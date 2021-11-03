## FOLD CHANGE NEIGHBOURHOOD ANALYSIS ON GRADED (HI MED LO) PLAQUES############

### LIBRARIES ###
libs <- c('ggplot2','RColorBrewer','reshape2','devtools',
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
amyloidGrades <- read_csv("/Volumes/KausaliaHD 1/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/CellEngineCurated/EzDeepcellPooled/alldataDiseasePooled_CD56vsAbeta_GatedPopulations/Abeta_Grades/AmyloidGates_CD56vsAb1724/CD56vsAb1724Grades.csv")

## loads Amyloid Neighbourhoods
plot_data <- read_csv(paste0(dataPath,"RealRuns/AmyloidPlaques_MatrixByIndividual_Objects/amyloidPlaques_AmyloidPlaque_counts_100px.csv" )) # _freqs_ or _counts_
outputFolder = paste0(dataPath,"DiseaseNeighbourhoodEZDeepcell/FoldChangeObjects/AmyloidNeighbourhood_100pixeldistance(11May21)/")
dir.create(outputFolder, showWarnings= F, recursive = T)

#plot_data = ceiling(plot_data)
plot_data_realunmelted = plot_data

plot_data_real = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_real) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_real$MantisPopulation)

## loads csv file of Random data that is generated from MIBIcreateObjectNeighborMatrixEzSegRandom.m
nrun = 1
fulldata = read_csv(paste0(dataPath,"RandomRuns/AmyloidP_MatrixByIndividual_Objects(11May21)/amyloidPlaques_AmyloidP_counts_100px_1.csv")) # _freqs_ or _counts_  
fulldata = as.matrix(fulldata)
# for (runnum in 2:nrun){
#   ## loads Amyloid Neighbourhoods
#   plot_data = read_csv(paste0(dataPath,"RandomRuns/AmyloidP_MatrixByIndividual_Objects/amyloidPlaques_AmyloidP_counts_100px_",runnum,".csv"))  
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
plot_data$MantisPopulation = ordered(plot_data$MantisPopulation, levels = c("neurons", "astrocyte_process",  "vessel_CD31_CD105", "microglia_process" ))  # in console to check if levels(plot_data$Meta)
levels(plot_data$MantisPopulation)

#plot_data = subset(plot_data, Object_id<3) 
## astrocytes Green "#106900", endothelial Blue"#0028D4",  microglia mango "#FFBF00",neurons DarkRed "#980A0A"...Color order follows order BEFORE custome ordering
palette = c("#980A0A","#106900","#0028D4","#FFBF00")

#Group ViolinPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = log2(Value))) +  # y =  MaxDist, MeanDist, MedianDist
  geom_violin(aes(fill = MantisPopulation), scale = "width", width = 0.6, adjust = 5) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(-6,6) +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Violin_AmyloidPlaques_FC_100Pixeldistance.pdf"), plot = p, width = 5, height = 5)

#Group ScatterPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = log2(Value))) + # y =  MaxDist, MeanDist, MedianDist
  geom_point(aes(color = MantisPopulation), cex = 0.5, alpha = 0.3, position = position_jitter(height = 0)) + # for freq no vertical jitter set to 0, for counts set 0.2
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(c(-5,5)) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(15) +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/Scatter_AmyloidPlaques_FC_100Pixeldistance.pdf"), plot = p, width = 7, height = 5)

#Boxplot with anova
myComparison = list(c("astrocyte_process", "microglia_process"), c("astrocyte_process", "neurons"), c("neurons", "microglia_process"), c("neurons", "vessel_CD31_CD105"), c("microglia_process", "vessel_CD31_CD105"), c("astrocyte_process", "vessel_CD31_CD105"))
plot_data$log2Value = log2(plot_data$Value)
plot_data = subset(plot_data, Value > 0)
res = aov(Value ~ MantisPopulation, plot_data)
TukeyHSD(res)
p = ggboxplot(plot_data,x = "MantisPopulation",y = "Value") +
  stat_compare_means(comparisons=myComparison, method = "anova", aes(label=..p.ad..)) 
ggsave(paste0(outputFolder, "/Amyloid_Boxplot_100pxDist.pdf"), plot = p, width = 7, height = 5)

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

##Heatmap to plot Amyloid Neighbourhood Flavours
set.seed(240214)
callback=function(x,data){ 
  hclust_dist = dist(t(as.matrix(plot_data)), method = "euclidean")
  hclust(hclust_dist, "ward.D2")
}

plot_data =  as.matrix(plot_data_realunmelted)
colnames_plotdata = paste0("p",plot_data[,1], "_objID", plot_data[,2])
rownames_plotdata = colnames(plot_data_realunmelted)[3:6]

## integrate Amyloid grade/gated cell engine population
tempdata = as.data.frame(plot_data)
tempdata$amyloidGrade = NA
  for (i in 1:nrow(tempdata)){
    amyloidObject =plot_data[i,]
    amyloidGrade = subset(amyloidGrades, point_id==amyloidObject[["Point_num"]] & obj_id==amyloidObject[["Object_num"]]) 
    if (nrow(amyloidGrade)>0){
      tempdata[i,"amyloidGrade"]=amyloidGrade[1,"file"]  
    }
}

table(tempdata$amyloidGrade)

tempdata$amyloidGrade = as.numeric(as.factor(tempdata$amyloidGrade))

plot_data = plot_data[,3:6]
plot_data = log2(plot_data/as.matrix(plot_data_randunmelted[,3:6]))
plot_data =  t(plot_data)

#To replace neg infinity or na with minimum value for each row
for (i in 1:nrow(plot_data)){
  plot_data[i,which(is.infinite(plot_data[i,]))] = min(plot_data[i,which(is.finite(plot_data[i,]))])
}

plot_data = as.data.frame(plot_data)
#plot_data = subset(plot_data, Object_num<3) 
colnames(plot_data) =  1:ncol(plot_data)
rownames(plot_data) = rownames_plotdata

plot_data[nrow(plot_data) + 1,] = tempdata$amyloidGrade
rownames(plot_data) = c(rownames_plotdata, "amyloidGrade")
plot_data = as.data.frame(plot_data)
plot_data = plot_data[,!is.na(plot_data[nrow(plot_data),])]

anno_data = plot_data["amyloidGrade", ]
anno_data = as.data.frame(t(anno_data))
anno_data$amyloidGrade = as.factor(anno_data$amyloidGrade)
levels(anno_data$amyloidGrade) =  c("hi","med","lo")
anno_color = list(amyloidGrade= c(hi = "#a80077", med = "#e1e73d", lo = "#56ab2f"))

#tally of hi med lo
table(anno_data$amyloidGrade)

plot_data =  plot_data[, order(plot_data["amyloidGrade",])]
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
  ylim(0,500)
ggsave(paste0(outputFolder,"/Bargraph_Amyloidnumberspercelltype.pdf"))

##NO SCALING##
# this Rocket color palette heatmaps is with amyloidGrade Clustering
pdf(paste0(outputFolder,"/Heatmap_AmyloidGrade_FC_100Pixeldistance.pdf"), width = 30, height = 10, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = F, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "none", main = "Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nPixelDistance = 50",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# these heatmaps is without amyloidGrade Clustering
#plot_data = plot_data[1:(nrow(plot_data) - 1),] 

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap1_AmyloidPlaques_FC_100Pixeldistances.pdf"), width = 30, height = 10, onefile = FALSE)
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "none", main = "Fold Change From \nRandomised microglia/astrocyte/vessel objects around amyloid plaque objects \nPixelDistance = 100",
            color = viridis(100, alpha=1, begin=0, end=1, direction=1, option = "D"), cellheight = 60)
dev.off()
dev.off()

##ROW SCALING##
# this plots Rocket color palette heatmaps with amyloidGrade Clustering
pdf(paste0(outputFolder,"/Heatmap_AmyloidGrade_FC_100Pixeldistancescaled.pdf"), width = 30, height = 10, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = F, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "row", main = "Fold Change From \nRandomised neuronal/microglia/astrocyte/vessel objects around tangle objects \nPixelDistance = 50",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# these heatmaps is without amyloidGrade Clustering
#plot_data = plot_data[1:(nrow(plot_data) - 1),] 

# this plots Rocket color palette heatmaps with lineage Clustering
pdf(paste0(outputFolder,"/Heatmap2_AmyloidPlaques_FC_100Pixeldistancerowscaled.pdf"), width = 30, height = 10, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "row", main = "Fold Change From \nRandomised microglia/astrocyte/neuronal objects around amyloid plaque objects \nPixelDistance = 100",
            col = colfunc(100), direction = 1, cellheight = 60)
dev.off()
dev.off()

# this plots Viridis color palette heatmaps with lineage Clustering
pdf(paste0(outputFolder,"/Heatmap1_AmyloidPlaques_FC_100Pixeldistancescaled.pdf"), width = 30, height = 10, onefile = FALSE)
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, annotation_col = anno_data , annotation_colors = anno_color, border_color = FALSE, scale = "row", main = "Fold Change From \nRandomised microglia/astrocyte/vessel objects around amyloid plaque objects \nPixelDistance = 100",
            color = viridis(100, alpha=1, begin=0, end=1, direction=1, option = "D"), cellheight = 60)
dev.off()
dev.off()



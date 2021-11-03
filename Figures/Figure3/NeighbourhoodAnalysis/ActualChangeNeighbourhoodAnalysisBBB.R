### LIBRARIES ###
libs <- c('ggplot2','RColorBrewer','reshape2','devtools',
          'pheatmap','flowCore','plyr','scales','coin','grid',
          'gridExtra','tidyr','randomcoloR','viridis','tidyverse', 'ggpubr')
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

## loads Vessel Neighbourhoods
plot_data <- read_csv(paste0(dataPath,"RealRuns/BBB_MatrixByIndividual_Objects/vessel_CD31_CD105_BB1_counts_25px.csv" )) # counts or vessel_MCT1_BBB2neighbor_counts_25px.csv, vessel_CD31_CD105_BBB1neighbor_counts_25px.csv
outputFolder = paste0(dataPath,"BBBNeighbourhoodEZDeepcell/RealCountsObject2/")
dir.create(outputFolder, showWarnings= F, recursive = T)

plot_data = ceiling(plot_data)
plot_data_realunmelted = plot_data

plot_data_real = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_real) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_real$MantisPopulation)


## loads csv file of Random data that is generated from MIBIcreateObjectNeighborMatrixEzSegRandom.m
nrun = 10
fulldata = read_csv(paste0(dataPath,"RandomRuns/BBB_MatrixByIndividual_Objects/vessel_CD31_CD105_BBB1_counts_25px_1.csv")) # counts of vessel_MCT1_BBB2_counts_25px_1.csv, vessel_CD31_CD105_BBB1_counts_25px.csv
fulldata = as.matrix(fulldata)
for (runnum in 2:nrun){
  ## loads Vessel Neighbourhoods
  plot_data = read_csv(paste0(dataPath,"RandomRuns/BBB_MatrixByIndividual_Objects/vessel_CD31_CD105_BBB1_counts_25px_",runnum,".csv"))  # vessel_CD31_CD105_BBB1neighbor_counts_25px.csv
  fulldata[,3:5] = fulldata[,3:5]+as.matrix(plot_data[,3:5])
}
fulldata[,3:5] = fulldata[,3:5]/nrun
fulldata[fulldata==0] = 1
fulldata = ceiling(fulldata)
plot_data = as.data.frame(fulldata)

plot_data_randunmelted = plot_data
plot_data_rand = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_rand) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_rand$MantisPopulation)
sum(plot_data_realunmelted$Object_num)

#make violin for Actual Mean values
plot_data = plot_data_real

##to calculate Fold change
#plot_data$Value = plot_data$Value/plot_data_rand$Value

#plot_data = subset(plot_data, Object_id<3) 
## astrocytes Green "#106900", endothelial Blue"#0028D4",  microglia mango "#FFBF00",neurons DarkRed "#980A0A"...Color order follows order BEFORE custome ordering
palette = c("#980A0A","#106900","#FFBF00")

#Group ViolinPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = Value)) +  # y =  MaxDist, MeanDist, MedianDist
  geom_violin(aes(fill = MantisPopulation), scale = "width", width = 0.6, adjust = 5) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(c(0,30))  +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/BBB1_Violin_25pxDist.pdf"), plot = p, width = 7, height = 5)

#Group ScatterPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = Value))+  # y =  MaxDist, MeanDist, MedianDist
  geom_point(aes(color = MantisPopulation), cex = 0.5, alpha = 0.3, position = position_jitter(height = 0.2)) + # for freq no vertical jitter set to 0, for counts set 0.2
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(c(0,30)) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(15) +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/BBB1_Scatter_25pxDist.pdf"), plot = p, width = 7, height = 5)

#Boxplot with T-test
myComparison = list(c("astrocyte_process", "microglia_process"),c("astrocyte_process", "neurons"), c("neurons", "microglia_process"))
p = ggboxplot(plot_data,x = "MantisPopulation",y = "Value") +
stat_compare_means(comparisons=myComparison, method = "t.test", aes(label=..p.ad..)) 
ggsave(paste0(outputFolder, "/BBB1_Boxplot_25pxDist.pdf"), plot = p, width = 7, height = 5)

#tally of object 
table(plot_data$MantisPopulation)

#Table with T-test and means
compare_means(Value~MantisPopulation,plot_data, method = 't.test')
library(dplyr)
plot_data %>%
  group_by(MantisPopulation) %>%
  summarise(mean = mean(Value), sd = sd(Value))

##Heatmap to plot Vessel Neighbourhood Flavours
set.seed(240214)
callback=function(x,data){
  hclust_dist = dist(t(as.matrix(plot_data)), method = "euclidean")
  hclust(hclust_dist, "ward.D2")
}

plot_data =  as.matrix(plot_data_realunmelted)
colnames_plotdata = paste0("p",plot_data[,1], "_objID", plot_data[,2])
rownames_plotdata = colnames(plot_data_realunmelted)[3:5]

plot_data = plot_data[,3:5]
#plot_data = plot_data/as.matrix(plot_data_randunmelted[,3:5]) ## To calculate fold Change
plot_data =  t(plot_data)

#plot_data = apply(plot_data, 2, scale)#marker scaled
plot_data[is.na(plot_data)] <- 0

plot_data = as.data.frame(plot_data)
#plot_data = subset(plot_data, Object_num<3) 
colnames(plot_data) =  1:ncol(plot_data)
rownames(plot_data) = rownames_plotdata

# this maps  Viridis color palette
breakList = c(seq(0,0.2,0.02),0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0)*max(plot_data) #c(0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.52,0.54,0.58,0.6,0.62,0.64,0.68,0.7,0.75,0.8,0.85,0.9,0.95,1.0)
pdf(paste0(outputFolder,"/Heatmap1_BBB1_25pxDist1.pdf"), width = 200, height = 5, onefile = FALSE)
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "Number of microglia/astrocyte/neuronal objects around vessel objects \nPixelDistance = 25",
            color = viridis(length(breakList), alpha=1, begin=0, end=1, direction=1, option = "D"), breaks = breakList, cellheight = 40)
dev.off()
dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap2_BBB1_25pxDist1.pdf"), width = 200, height = 5, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "Number of microglia/astrocyte/neuronal objects around vessel objects \nPixelDistance = 25",
            col = colfunc(100), direction = 1, cellheight = 40)
dev.off()
dev.off()


tabledata = plot_data
table(plot_data_real$Value) # num
table(tabledata[1,]>tabledata[2,]) # TRUE shows number of astrocytes greater than neurons, FALSE is shows number of astrocytes lessthan and equals to than neurons


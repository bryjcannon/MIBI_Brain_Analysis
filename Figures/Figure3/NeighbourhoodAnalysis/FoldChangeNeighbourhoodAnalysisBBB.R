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
setwd("/Volumes/KausaliaHD/R4R")
source("/Volumes/KausaliaHD/R4R/General/readWriteCytofFCS.R")

## loads csv file of Real data that is generated from MIBIcreateObjectNeighborMatrixEzSeg.m
dataPath = "/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/Pooled/" 

## loads Vessel Neighbourhoods
plot_data <- read_csv(paste0(dataPath,"RealRuns/BBB_MatrixByIndividual_Objects/vessel_CD31_CD105_BB1_counts_25px.csv" )) #
outputFolder = paste0(dataPath,"BBBNeighbourhoodEZDeepcell/FoldChangeObject2/")
dir.create(outputFolder, showWarnings= F, recursive = T)

plot_data = ceiling(plot_data)
plot_data_realunmelted = plot_data

## to calculate ratios of cell types aroung objetcs rather than absolute counts (for real data)
for(n in 1:nrow(plot_data)){
  cellcount = sum(plot_data[n,3:5])
  plot_data[n,3:5] = plot_data[n,3:5]/cellcount
}
##

plot_data_real = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_real) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_real$MantisPopulation)

## loads csv file of Random data that is generated from MIBIcreateObjectNeighborMatrixEzSegRandom.m
nrun = 10
fulldata = read_csv(paste0(dataPath,"RandomRuns/BBB_MatrixByIndividual_Objects/vessel_CD31_CD105_BBB1_counts_25px_1.csv")) # 
fulldata = as.matrix(fulldata)
for (runnum in 2:nrun){
## loads Vessel Neighbourhoods
plot_data = read_csv(paste0(dataPath,"RandomRuns/BBB_MatrixByIndividual_Objects/vessel_CD31_CD105_BBB1_counts_25px_",runnum,".csv"))  
fulldata[,3:5] = fulldata[,3:5]+as.matrix(plot_data[,3:5])
}
fulldata[,3:5] = fulldata[,3:5]/nrun
fulldata[fulldata==0] = 1
fulldata = ceiling(fulldata)
plot_data = as.data.frame(fulldata)

## to calculate ratios of cell types aroung objetcs rather than absolute counts (for random data)
for(n in 1:nrow(plot_data)){
  cellcount = sum(plot_data[n,3:5]) # total cells count around every endothelial cells
  plot_data[n,3:5] = plot_data[n,3:5]/cellcount 
}
##
plot_data_randunmelted = plot_data
plot_data_rand = melt(plot_data, id.vars = c("Point_num", "Object_num"))
colnames(plot_data_rand) = c("Point", "Object_id", "MantisPopulation", "Value")
table(plot_data_rand$MantisPopulation)

#make violin for fold change
plot_data = plot_data_real
plot_data$Value = plot_data$Value/plot_data_rand$Value

#plot_data = subset(plot_data, Object_id<3) 
## astrocytes Green "#106900", endothelial Blue"#0028D4",  microglia mango "#FFBF00",neurons DarkRed "#980A0A"...Color order follows order BEFORE custome ordering
palette = c("#980A0A","#106900","#FFBF00")

#Group ViolinPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = log2(Value))) +  # y =  MaxDist, MeanDist, MedianDist
  geom_violin(aes(fill = MantisPopulation), scale = "width", width = 0.6, adjust = 5) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_fill_manual(values = palette) +
  xlab("") + ylab("") +
  ylim(-3,4) +
  ggtitle("") +
  theme_classic(12) +
  #theme(legend.position = "none") +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/BBB1_Violin_25pxDist_FCRatios.pdf"), plot = p, width = 7, height = 5)

#Group ScatterPlot (Distances per cell type per object)
p = ggplot(plot_data,aes(x = MantisPopulation, y = log2(Value) )) +   # y =  MaxDist, MeanDist, MedianDist log2(Value)
  geom_point(aes(color = MantisPopulation), cex = 0.5, alpha = 1, position = position_jitter(height = 0.2)) + # for freq no vertical jitter set to 0, for counts set 0.2
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.size = 0) +
  scale_color_manual(values = palette) +
  xlab("") + ylab("") +
  ggtitle(paste0("")) +
  ylim(-3,4) + # watch for warnings to see how many outliers were removed in each plot.
  theme_classic(15) +
  theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10),
        axis.title.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(angle = 270, hjust = 0, vjust = 0.5))
ggsave(paste0(outputFolder, "/BBB1_Scatter_25pxDist_FCRatios.pdf"), plot = p, width = 7, height = 5)

#Boxplot with T-test
myComparison = list(c("astrocyte_process", "microglia_process"),c("astrocyte_process", "neurons"), c("neurons", "microglia_process"))
p = ggboxplot(plot_data,x = "MantisPopulation",y = "Value") +
  stat_compare_means(comparisons=myComparison, method = "t.test", aes(label=..p.ad..)) 
ggsave(paste0(outputFolder, "/BBB1_Boxplot_25pxDist.pdf"), plot = p, width = 7, height = 5)

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
plot_data = plot_data/as.matrix(plot_data_randunmelted[,3:5])
plot_data =  t(plot_data)

plot_data = apply(plot_data, 2, scale) # marker scaled
plot_data[is.na(plot_data)] <- 0

plot_data = as.data.frame(plot_data)
#plot_data = subset(plot_data, Object_num<3) 
colnames(plot_data) =  1:ncol(plot_data)
rownames(plot_data) = rownames_plotdata

# this maps  Viridis color palette
pdf(paste0(outputFolder,"/Heatmap1_BBB1_25pxDist_FCRatios.pdf"), width = 20, height = 5, onefile = FALSE)
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "Fold Change From \nRandomised microglia/astrocyte/neuronal objects around vessel objects \nPixelDistance = 25",
            color = viridis(100, alpha=1, begin=0, end=1, direction=1, option = "D"), cellheight = 40)
dev.off()
dev.off()

# this maps Rocket color palette
pdf(paste0(outputFolder,"/Heatmap2_BBB1_25pxDist_FCRatios.pdf"), width = 20, height = 5, onefile = FALSE)
colfunc <- colorRampPalette(c("#000099","#E0E0E0","#990000")) #"#1565C0" Blue,"#E0E0E0"Gray,"#C62828" Red
p=pheatmap (plot_data, cluster_rows = F, cluster_cols = T, clustering_callback = callback, border_color = FALSE, scale = "none", main = "Fold Change From \nRandomised microglia/astrocyte/neuronal objects around vessel objects \nPixelDistance = 25",
            col = colfunc(100), direction = 1, cellheight = 40)
dev.off()
dev.off()

tabledata = plot_data
table(plot_data_real$Value) # num
table(tabledata[1,]>tabledata[2,]) # TRUE shows number of astrocytes greater than neurons, FALSE is shows number of astrocytes lessthan and equals to than neurons
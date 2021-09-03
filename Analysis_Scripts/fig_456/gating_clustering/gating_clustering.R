#Gating & Clustering

# START OF GATING
# all gated data set together
master_data_all$cluster_id = 0

gateNames = c("Ungated","oligodendrocytes","astrocytes","microglia","endothelial","Dopaminergic-TH","Inhibitory-VGAT")

## To Gate Population. ## if not running FlowSOM
gatedIndexes = c()

## visualization
master_data_all_cuberoot %>%
  filter(MAP2 >= 0 & MAP2 <= 2) %>%
  filter(GFAP >= 0 & GFAP <= 5) %>%
  filter(obj_type_id %in% "cells") %>%
  ggplot() +
  aes(x = MAP2, y = GFAP) +
  geom_point(size = 1L, colour = "#0c4c8a") +
  geom_density2d(colour = "#FF7F00") +
  theme_minimal()

## viz2
master_data_all_cuberoot %>%
  filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter(Iba1 > 0 & MAP2 > 0)
  ggplot() +
  aes(x = Iba1, y = MAP2) +
  geom_point(size = 0.25L, colour = "#0c4c8a") +
  geom_density2d(colour = "#FF7F00") +
  theme_minimal()

# data transform
master_data_all_cuberoot <- master_data_all
master_data_all_cuberoot[panel] <- (master_data_all_cuberoot[panel]*10000)^(1/3)
filter(master_data_all, obj_type_id == 'cells') %>% select(cluster_id) %>% table()

# gating cells based on biaxial plots (versions of above)

  # oligodendrocytes
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG > 10 & MBP > 12)
  master_data_all_cuberoot[spatial_id == gate]$cluster_id = 1
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & MAG > 10 & MBP > 12))
  #data_indexes = which(with(master_data_all_cuberoot, gate)) #from dmitry
  #data_indexes = setdiff(data_indexes, gatedIndexes) #from dmitry
  master_data_all[data_indexes, "cluster_id"] = 1
  #gatedIndexes = c(gatedIndexes, data_indexes) #from dmitry
  
  # astrocytes
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) %>%
    filter(GFAP > 8 & MAP2 < 10) %>% filter(GFAP / MAP2 > 1.2) %>% dim()
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & MAG <= 10 & MBP <= 12 & ((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) & ((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) &
                            (GFAP > 8 & MAP2 < 10)  & (GFAP / MAP2 > 1.2)))
  master_data_all[data_indexes, "cluster_id"] = 2
  
  # microglia
  gate <-  master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 >= 18) | (MCT1 <= 7 | MAP2 >= 18)) %>% filter((CD45 > 7 & MAP2 < 15) | (Iba1 > 10))
  data_indexes = which(with(gate))
  master_data_all[data_indexes, "cluster_id"] = 3
  
  # endothelial
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 > 6 & MAP2 < 15) | (CD105 > 7 & MAP2 < 18) | (MCT1 > 7 | MAP2 < 18))
  data_indexes = which(with(gate)) 
  master_data_all[data_indexes, "cluster_id"] = 4
  
  # neurons
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) %>%
    filter(GFAP <= 8 | MAP2 >= 10) %>% filter(GFAP / MAP2 <= 1.2)
  data_indexes = which(with(gate)) 
  master_data_all[data_indexes, "cluster_id"] = 5








#####################################################################







data_indexes = which(with(master_data_all, TH > 0.2)) # TH
data_indexes = setdiff(data_indexes, gatedIndexes)
master_data_all[data_indexes, "cluster_id"] = 5 # to assign cluster ID for dopaminergic neurons
#plot_histograms(gateNames[4])
gatedIndexes = c(gatedIndexes, data_indexes)

data_indexes = which(with(master_data_all, VGAT > 0.2)) # VGAT
data_indexes = setdiff(data_indexes, gatedIndexes)
master_data_all[data_indexes, "Meta"] = 6 # to assign cluster ID for VGAT neurons
#plot_histograms(gateNames[4])
gatedIndexes = c(gatedIndexes, data_indexes)
##
##
##
table(master_data_all$cluster_id)

## if running FlowSOM
fcsExprs <- flowCore::flowFrame(as.matrix(master_data_all[panel]))
suppressWarnings(flowCore::write.FCS(fcsExprs, 'allDataFcs.fcs', what = "numeric"))

library(FlowSOM)
library(Rtsne)
seed = 123

numClusters = 15 # dont change cos original flowSom cluster
numGates = 6
palette <- distinctColorPalette(numClusters)
exclude = c("C12", "Na23", "Si28", "Ca40", "Background", "Ta181", "Au197", "empty113" )
clusteringMarkers = c("VGLUT1", "VGLUT2")#,"Iba1","CD45","MCT1","CD31","MBP","MAG","MAP2","TH", "VGLUT1", "VGLUT2", "VGAT","GFAP", "PolyubiK63", "8OHGuano","ApoE","beta142","pTDP43", "PHF", "PanAbeta"
clusteringMarkers = c("Iba1","CD45","MCT1","CD31","MBP","MAG","MAP2","Reelin","GFAP")
#markers = setdiff(TMAPanel$Label, exclude)

fSOM <- FlowSOM(fcsExprs,
                compensate = F, transform = F,scale = T,
                colsToUse = clusteringMarkers, nClus = numClusters, seed = seed)

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
##
##
##
###  
# initial cluster number  
num_clusters <- 5
# color palette of choice, here using c25 from brain_data_pkg
palette <- c25

clustering_markers <- c('HistoneH3Lyo', 'MAP2', 'VGAT', 'VGLUT1', 'VGLUT2', 'X8OHGuano', 'GFAP', 'Iba1', 'CD45', 'MCT1', 'CD31')
exclude <- c("C12", "Na23", "Si28", "Ca40", "Background", "Ta181", "Au197", "empty113" )

markers <- setdiff(panel, exclude)

# create a flowFrame from previously imported, transformed data to use in FlowSOM calculation. Currently set to only cluster on 'cell' objects from deepCell data.
flowFrame_cluster <-  new("flowFrame", exprs = as.matrix(subset(master_data, obj_type_id == 'cell', select = panel)))

# run FlowSOM on the above flowFrame which has already had transformation performed
fSOM <- FlowSOM(flowFrame_cluster,
                compensate = F, scale = F, colsToUse = clustering_markers, nClus = num_clusters, seed = 140214)

table(fSOM$metaclustering) 
metaClustering <- fSOM$metaclustering

fSOM_clustering <- data.frame(fSOM$FlowSOM$map$mapping)
colnames(fSOM_clustering) <- c("Cluster", "Value")
fSOM_clustering$Meta <- as.numeric(fSOM$metaclustering[fSOM_clustering[,1]])

fSOMcodes = fSOM$FlowSOM$map$codes
write.csv (fSOMcodes, paste0(outputFolder,'/fSOM_FlowSOM_map_codes.csv'))
fSOMcodesMean = apply(fSOMcodes,2,mean)
write.csv(fSOMcodesMean, paste0(outputFolder,'/mean_of_fSOM_FlowSOM_map_codes.csv'))

# VIZ ABOVE
plot_data = data.frame(table(fSOM_clustering$Meta))
plot_ly(x = plot_data$Var1, y = plot_data$Freq, type = 'bar')
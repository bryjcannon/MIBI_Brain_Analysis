#Gating & Clustering

# START OF GATING
# all gated data set together
master_data_all_cuberoot$cluster_id = 0
# add gate ids for ez_segmented objects
data_indexes = which(with(master_data_all_cuberoot, obj_type_id == "tau-pathies"))
master_data_all_cuberoot[data_indexes, "cluster_id"] = 1
data_indexes = which(with(master_data_all_cuberoot, obj_type_id == "amyloid-pathies"))
master_data_all_cuberoot[data_indexes, "cluster_id"] = 2
data_indexes = which(with(master_data_all_cuberoot, obj_type_id == "microglia-processes"))
master_data_all_cuberoot[data_indexes, "cluster_id"] = 3
data_indexes = which(with(master_data_all_cuberoot, obj_type_id == "vessels"))
master_data_all_cuberoot[data_indexes, "cluster_id"] = 4

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
  filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter(Iba1 > 0 & MAP2 > 0) %>%
  ggplot() +
  aes(x = Iba1, y = MAP2) +
  geom_point(size = 0.25L, colour = "#0c4c8a") +
  geom_density2d(colour = "#FF7F00") +
  theme_minimal()

# size viz
master_data_all_cuberoot %>%
  filter(obj_type_id %in% "cells") %>%
  ggplot() +
  aes(MajorAxisLength / MinorAxisLength) +
  geom_histogram(size = 1L, colour = "#0c4c8a", binwidth = 1) +
  xlim(0,10) +
  theme_minimal()

# gating cells based on biaxial plots (versions of above)

  # oligodendrocytes
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG > 12 & MBP > 12)
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (MAG > 12 & MBP > 12)))
  #data_indexes = which(with(master_data_all_cuberoot, gate)) #from dmitry
  #data_indexes = setdiff(data_indexes, gatedIndexes) #from dmitry
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 5
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "oligodendrocytes" # currently doesn't work due to obj_type_id as a factor variable
  #gatedIndexes = c(gatedIndexes, data_indexes) #from dmitry
  
  # astrocytes
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 13) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) %>%
    filter(GFAP > 8 & MAP2 < 10) %>% filter(GFAP / MAP2 > 1.2) %>% dim()
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (MAG <= 10 & MBP <= 13) & ((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) & ((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) &
                            (GFAP > 8 & MAP2 < 10)  & (GFAP / MAP2 > 1.2)))
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 6
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "astrocytes"
  
  # combine astrocytes and oligodendrocytes into non-immune neuroglia bucket
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (cluster_id == 5 | cluster_id == 6)))
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 13
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "non-immune neuroglia"
  
  # microglia
  gate <-  master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 13) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 >= 18) | (MCT1 <= 7 | MAP2 >= 18)) %>% filter((CD45 > 7 & MAP2 < 15) | (Iba1 > 10))
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (MAG <= 10 & MBP <= 13) & ((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 >= 18) | (MCT1 <= 7 | MAP2 >= 18)) & ((CD45 > 7 & MAP2 < 15) | (Iba1 > 10))))
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 7
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "microglia"
  
  # endothelial
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 13) %>% filter((CD31 > 6 & MAP2 < 15) | (CD105 > 7 & MAP2 < 18) | (MCT1 > 7 | MAP2 < 18))
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (MAG <= 10 & MBP <= 13) & ((CD31 > 6 & MAP2 < 15) | (CD105 > 7 & MAP2 < 18) | (MCT1 > 7 | MAP2 < 18)))) 
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 8
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "endothelial"
  
  # neurons
  gate <- master_data_all_cuberoot %>% filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 13) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) %>%
    filter(GFAP <= 8 | MAP2 >= 10) %>% filter(GFAP / MAP2 <= 1.2)
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (MAG <= 10 & MBP <= 13) & ((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) & ((CD45 <= 7 | MAP2 >= 15) | (Iba1 <= 10)) &
                              (GFAP <= 8 | MAP2 >= 10) & (GFAP / MAP2 <= 1.2))) 
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 9
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "neurons"

  # gate out anything with low histone
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & (HistoneH3Lyo <= 5 | obj_size > 750 | MajorAxisLength / MinorAxisLength > 5)))
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 10
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "lowH3"
  
  # gate out endothelial cells at the edge of tissue (near gold)
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & cluster_id == '8' & run_type_id == 'Ctrl1' & point_id %in% c(22,24,25,27,26,55,56,85,86))) 
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & cluster_id == '8' & run_type_id == 'Ctrl2' & point_id %in% c(15,16,45,46,75,76,105)))
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & cluster_id == '8' & run_type_id == 'Ctrl3' & point_id %in% c(1,50)))
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & cluster_id == '8' & run_type_id == 'HiAD' & point_id %in% c(1,28,29,160,169,196,195)))
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & cluster_id == '8' & run_type_id == 'MedAD' & point_id %in% c(156-1,157-1,182-1,183-1,184-1,208-1,207-1,209-1,210-1,211-1,234-1,233-1,232-1,235-1,236-1,237-1,260-1,259-1,258-1,257-1,1-1,2-1,26-1,27-1,25-1)))
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 11
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "edge_endo"
  
  # gate out unwanted border regions in MedAD
  data_indexes = which(with(master_data_all_cuberoot, obj_type_id == 'cells' & run_type_id == 'MedAD' & point_id %in% c(156-1,157-1,182-1,183-1,184-1,208-1,207-1,209-1,210-1,211-1,234-1,233-1,232-1,239-1,235-1,236-1,237-1,238-1,239-1,240-1,241-1,260-1,259-1,258-1,257-1,256-1,255-1,1-1,2-1,26-1,27-1,25-1,254-1,253-1)))
  master_data_all_cuberoot[data_indexes, "cluster_id"] = 12
  master_data_all_cuberoot[data_indexes, "obj_type_id"] = "discard"
  
  # check cluster_ids
  filter(master_data_all_cuberoot, obj_type_id == 'cells') %>% select(cluster_id) %>% table()
  select(master_data_all_cuberoot, cluster_id) %>% table()




#############################________################################
# output df to csvs for matlab use - CytoMAP
  master_data_all_cuberoot_c <- master_data_all_cuberoot
  #converting levels to num for csv
  levels(master_data_all_cuberoot_c$obj_type_id) <- c('cells', 'amyloid-pathies', 'tau-pathies', 'microglia-processes', 'vessels')
  
  levels(master_data_all_cuberoot_c$Subregions) <- c("Alveus", "Au", "CA1", "CA2", "CA3", "CA4", "DG", "EC", "Fimbria", "ML", "PreSubiculum", "SRL", "Subiculum", "TH-LV")
  
  levels(master_data_all_cuberoot_c$run_type_id) <- c('HiAD_modded', 'MedAD_modded', 'Ctrl_modded', 'Ctrl_modded', 'Ctrl_modded')
  run_levels <- levels(master_data_all_cuberoot_c$run_type_id)
  
  master_data_all_cuberoot_cc <- mutate_if(master_data_all_cuberoot_c, is.factor, ~ as.numeric(.x))
  
  cluster_names <- c('tau-pathies', 'amyloid-pathies', 'microglia-processes', 'vessels', 'oligodendrocytes', 'astrocytes', 'microglia', 'endothelial', 'neurons')
  
  for (level in 1:length(run_levels)) {
    path_name <- paste0("/Volumes/BryJC_Stanford/paper1_analysis/Fig6/single_object_pipeline/cytomap_files/", run_levels[level])
    dir.create(path_name, showWarnings = FALSE)
    for (cluster_num in 1:length(cluster_names)) {
      data_to_write <- master_data_all_cuberoot_cc %>% filter(run_type_id == level & cluster_id == cluster_num) %>% na.omit()
      write_csv(data_to_write, paste0(path_name,'/',cluster_names[cluster_num],'.csv'))
    }
  }

# output df to txts for python use - SVCA
  cluster_names <- c('tau-pathies', 'amyloid-pathies', 'microglia-processes', 'vessels', 'oligodendrocytes', 'astrocytes', 'microglia', 'endothelial', 'neurons')
  use_clusters <- c(1, 2, 5, 6, 7, 8, 9)
  master_data_all_cuberoot_c <- filter(master_data_all_cuberoot, cluster_id == 1 | cluster_id == 2 | cluster_id == 5 | cluster_id == 6 | cluster_id == 7 | cluster_id == 8 | cluster_id == 9)
  levels(master_data_all_cuberoot_c$run_type_id) <- c('HiAD', 'MedAD', 'Ctrl', 'Ctrl', 'Ctrl')
  run_levels <- levels(master_data_all_cuberoot_c$run_type_id)
  
  for (level in run_levels) {
    path_name <- paste0("/Volumes/BryJC_Stanford/paper1_analysis/Fig6/single_object_pipeline/svca_files/", level)
    run_subregions <- levels(filter(master_data_all_cuberoot_c, run_type_id == 'HiAD')$Subregions)
    for (subr in run_subregions) {
      path_name_subr <- paste0(path_name, '_', subr)
      dir.create(path_name_subr, showWarnings = FALSE)
      # make expression table
      data_to_write <- master_data_all_cuberoot_c %>% filter(run_type_id == level)
      panel_data <- data_to_write[panel][5:43]
      write.table(panel_data, paste0(path_name_subr,'/expression.txt'), sep = " ", row.names = FALSE, col.names = TRUE)
      # make position table
      data_to_write <- master_data_all_cuberoot_c %>% filter(run_type_id == level)
      x_y <- data_to_write[,c(64,65)]
      write.table(x_y, paste0(path_name_subr,'/positions.txt'), sep = ",", row.names = FALSE, col.names = FALSE)
    }
  }

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
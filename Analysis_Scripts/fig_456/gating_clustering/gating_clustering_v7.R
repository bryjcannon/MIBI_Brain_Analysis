#### Manual Gating & Clustering ####
  gated_data <- data_normalized
  
  ## visualization
  gated_data %>%
    filter(MAG > 0 & MAG <= 2) %>%
    filter(GFAP > 0 & GFAP <= 5) %>%
    filter(obj_type_id %in% "cells") %>%
    ggplot() +
    aes(x = MAG, y = GFAP, colour = MBP) +
    geom_point(size = 1L ) +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal()
  
  ## viz2
  gated_data %>%
    filter(obj_type_id == 'cells') %>% filter(MAG <= 10 & MBP <= 12) %>% filter((CD31 <= 6 | MAP2 >= 15) | (CD105 <= 7 | MAP2 > 18) | (MCT1 <= 7 | MAP2 > 18)) %>% filter(Iba1 > 0 & MAP2 > 0) %>%
    ggplot() +
    aes(x = Iba1, y = MAP2) +
    geom_point(size = 0.25L, colour = "#0c4c8a") +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal()
  
  # size viz
  gated_data %>%
    filter(obj_type_id %in% "amyloid-pathies") %>%
    ggplot() +
    aes(MajorAxisLength / MinorAxisLength) +
    geom_histogram(size = 1L, colour = "#0c4c8a", binwidth = 1) +
    xlim(0,25) +
    theme_minimal()
  
  gated_data %>%
    filter(obj_type_id %in% "tau-pathies") %>%
    ggplot() +
    aes(obj_size) +
    geom_histogram(size = 1L, colour = "#0c4c8a", binwidth = 1) +
    xlim(0,10) +
    theme_minimal()
  
  # neuron gating system
  gated_data %>%
    filter(SERT > 0 & SERT < 5) %>%
    filter(GFAP > 0 & GFAP < 5) %>%
    filter(obj_type_id %in% "cells") %>%
    filter((MAG < 0.4 & MBP < 0.45) & ((CD31 < 0.25 | MAP2 > 0.5) | (CD105 < 0.3 | MAP2 > 0.6) | (MCT1 < 0.275 | MAP2 > 0.6)) & ((CD45 < 0.275 | MAP2 > 0.6) | (Iba1 < 0.5))) %>%
    ggplot() +
    aes(x = SERT, y = GFAP) +
    geom_point(size = 1L, colour = "#0c4c8a") +
    ylim(0, 0.4) +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal()
  
  
  # START OF GATING
  
  # all gated data set together
  gated_data$cluster_id = 0
  # add gate ids for ez_segmented plaques, tangles
  data_indexes = which(with(gated_data, obj_type_id == "tau-pathies" & obj_size > 4))
  gated_data[data_indexes, "cluster_id"] = 1
  data_indexes = which(with(gated_data, obj_type_id == "amyloid-pathies" & obj_size > 9))
  gated_data[data_indexes, "cluster_id"] = 2
  data_indexes = which(with(gated_data, (obj_type_id == "tau-pathies" | obj_type_id == "amyloid-pathies") & ((Orientation < 5 & Orientation > -5) & Eccentricity > 0.8)))
  gated_data[data_indexes, "cluster_id"] = 16
  # select mglia-processes with size above mean 
  gated_data %>% filter(obj_type_id == "microglia-processes") %>% pull(obj_size) %>% median()
  data_indexes = which(with(gated_data, obj_type_id == "microglia-processes" & obj_size >= 14))
  gated_data[data_indexes, "cluster_id"] = 3
  # select mglia-processes with size above mean 
  gated_data %>% filter(obj_type_id == "vessels") %>% pull(obj_size) %>% median()
  data_indexes = which(with(gated_data, obj_type_id == "vessels" & obj_size >= 35))
  gated_data[data_indexes, "cluster_id"] = 4
  
  # gating cells based on biaxial plots (versions of above)
  
    # oligodendrocytes
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (MAG > 0.4 | MBP > 0.45)))
    gated_data[data_indexes, "cluster_id"] = 5
    
    # endothelial
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (MAG <= 0.4 & MBP <= 0.45) & ((CD31 > 0.25 & MAP2 < 0.5) | (CD105 > 0.3 & MAP2 < 0.6) | (MCT1 > 0.275 & MAP2 < 0.6)))) 
    gated_data[data_indexes, "cluster_id"] = 8
    
    # microglia
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (MAG <= 0.4 & MBP <= 0.45) & ((CD31 <= 0.25 | MAP2 >= 0.5) | (CD105 <= 0.3 | MAP2 >= 0.6) | (MCT1 <= 0.275 | MAP2 >= 0.6)) & ((CD45 > 0.275 & MAP2 < 0.6) | (Iba1 > 0.5))))
    gated_data[data_indexes, "cluster_id"] = 7
    
    # astrocytes
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (MAG <= 0.4 & MBP <= 0.45) & ((CD31 <= 0.25 | MAP2 >= 0.5) | (CD105 <= 0.3 | MAP2 >= 0.6) | (MCT1 <= 0.275 | MAP2 >= 0.6)) & ((CD45 <= 0.275 | MAP2 >= 0.6) | (Iba1 <= 0.5)) &
                              (GFAP > 0.2 & MAP2 < 0.5)))
    gated_data[data_indexes, "cluster_id"] = 6
    
    # neurons
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (MAG < 0.4 & MBP < 0.45) & ((CD31 < 0.25 | MAP2 > 0.5) | (CD105 < 0.3 | MAP2 > 0.6) | (MCT1 < 0.275 | MAP2 > 0.6)) & ((CD45 < 0.275 | MAP2 > 0.6) | (Iba1 < 0.5)) &
                                (GFAP <= 0.2 | (MAP2 >= 0.5 | MFN2 >= 0.3 | TH >= 0.37 | SERT >= 0.14))))
    gated_data[data_indexes, "cluster_id"] = 9
    
    # combine astrocytes and oligodendrocytes into non-immune neuroglia bucket
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (cluster_id == 5 | cluster_id == 6)))
    gated_data[data_indexes, "cluster_id"] = 10
  
    # gate out anything with low histone
    data_indexes = which(with(gated_data, obj_type_id == 'cells' & (HistoneH3Lyo <= 0.125 | obj_size > 800 | MajorAxisLength / MinorAxisLength > 5)))
    gated_data[data_indexes, "cluster_id"] = 11
    
    # gate out endothelial cells, disease objects at the edge of tissue (near gold)
    data_indexes = which(with(gated_data, (cluster_id == 8 | cluster_id == 1 | cluster_id == 2 | cluster_id == 3 | cluster_id == 4) & run_type_id == 'Ctrl1' & point_id %in% c(22,23,24,25,26,27,55,56,85,86))) 
    gated_data[data_indexes, "cluster_id"] = 12
    data_indexes = which(with(gated_data, (cluster_id == 8 | cluster_id == 1 | cluster_id == 2 | cluster_id == 3 | cluster_id == 4) & run_type_id == 'Ctrl2' & point_id %in% c(15,16,45,46,75,76,105)))
    gated_data[data_indexes, "cluster_id"] = 12
    data_indexes = which(with(gated_data, (cluster_id == 8 | cluster_id == 1 | cluster_id == 2 | cluster_id == 3 | cluster_id == 4) & run_type_id == 'Ctrl3' & point_id %in% c(1,50)))
    gated_data[data_indexes, "cluster_id"] = 12
    data_indexes = which(with(gated_data, (cluster_id == 8 | cluster_id == 1 | cluster_id == 2 | cluster_id == 3 | cluster_id == 4) & run_type_id == 'HiAD' & point_id %in% c(1,28,29,168,169,196,195)))
    gated_data[data_indexes, "cluster_id"] = 12
    data_indexes = which(with(gated_data, (cluster_id == 8 | cluster_id == 1 | cluster_id == 2 | cluster_id == 3 | cluster_id == 4) & run_type_id == 'MedAD' & point_id %in% c(156-1,157-1,182-1,183-1,184-1,208-1,207-1,209-1,210-1,211-1,234-1,233-1,232-1,235-1,236-1,237-1,260-1,259-1,258-1,257-1,1-1,2-1,26-1,27-1,25-1)))
    gated_data[data_indexes, "cluster_id"] = 12
    
    # gate out unwanted border regions in MedAD
    data_indexes = which(with(gated_data, run_type_id == 'MedAD' & point_id %in% c(156-1,157-1,182-1,183-1,184-1,208-1,207-1,209-1,210-1,211-1,234-1,233-1,232-1,239-1,235-1,236-1,237-1,238-1,239-1,240-1,241-1,260-1,259-1,258-1,257-1,256-1,255-1,1-1,2-1,26-1,27-1,25-1,254-1,253-1)))
    gated_data[data_indexes, "cluster_id"] = 13
    
    # assign any mglia or endothelial cells near projections (7) / vessels (8) respectively to null
    for (cell in seq_along(gated_data)) {
      if (gated_data[cell,]$cluster_id == 8) {
        x = gated_data[cell,]$x_centroid
        y = gated_data[cell,]$y_centroid
        #filter projections by -10 & + 10 in x and y direction, cluster_id's 3 and 4 to match 7 and 8 above, respectively.
        data_indexes = which(with(gated_data, (cluster_id == 4 & (x_centroid < abs(x-MajorAxisLength/2) | x_centroid < abs(x+MajorAxisLength/2) & (y_centroid < abs(x-MajorAxisLength/2) | y_centroid < abs(x+MajorAxisLength/2))))))
        gated_data[data_indexes, "cluster_id"] = 14
      }
    }
                         
    # set ungated objects and cells to non-zero cluster_id
    data_indexes = which(with(gated_data, cluster_id == 0))
    gated_data[data_indexes, "cluster_id"] = 15
    
    # check cluster_ids
    filter(gated_data, obj_type_id == 'cells') %>% select(cluster_id) %>% table()
    select(gated_data, cluster_id) %>% table()




#############################________################################

#############################________################################
#### Clustering Methods ####
  
  #### Dmitry code ####
    ## if running FlowSOM
    fcsExprs <- flowCore::flowFrame(as.matrix(data_normalized[panel]))
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

    
    
  #### Bryan Code ####
    
    # initial cluster number  
    num_clusters <- 30
    # color palette of choice, here using c25 from brain_data_pkg
    palette <- c25
    
    #clustering_markers <- c('TH', 'Calretinin', 'Calbindin', 'Parvalbumin', 'PanGAD6567', 'VGAT', 'VGLUT2', 'MAP2', 'MBP', 'MAG', 'GFAP', 'Iba1', 'CD45', 'MCT1', 'CD31', 'CD105')
    clustering_markers <- c('MAP2', 'MBP', 'MAG', 'GFAP', 'Iba1', 'CD45', 'MCT1', 'CD31', 'CD105')
    exclude <- c("C12", "Na23", "Si28", "Ca40", "Background", "Ta181", "Au197", "empty113" )
    
    markers <- setdiff(panel, exclude)
    
    # create a flowFrame from previously imported, transformed data to use in FlowSOM calculation. Currently set to only cluster on 'cell' objects from deepCell data.
    fsom_clustering_data <- filter(gated_data, obj_type_id == "cells" & cluster_id %in% c(7:10))
      #used to compare manual clustering on un-normalized data
      master_data_all_cuberoot_orig_cid <- master_data_all_cuberoot %>% rename(og_cid = cluster_id)
      cluster_comparisons <- inner_join(fsom_clustering_data, master_data_all_cuberoot_orig_cid, by = "spatial_id", copy = TRUE) 
      #convert all factors to numeric                                  
      fsom_clustering_data %>% select(cluster_id) %>% table()
      fsom_clustering_data <- mutate_if(fsom_clustering_data, is.factor, ~ as.numeric(.x))
    
    flowFrame_cluster <- flowCore::flowFrame(as.matrix(fsom_clustering_data))
    
    
    # run FlowSOM on the above flowFrame which has already had transformation performed
    fSOM2 <- FlowSOM(flowFrame_cluster,
                    compensate = F, scale = F, colsToUse = clustering_markers, nClus = num_clusters, seed = 123, xdim=20, ydim=20)
    
    table(fSOM2$metaclustering) 
    metaClustering_consensus <- metaClustering_consensus(fSOM2$FlowSOM$map$codes,k=5)
    metaClustering2 <- fSOM2$metaclustering
    
    # isolate cluster and metaclustering data
    fSOM2_clustering <- data.frame(fSOM2$FlowSOM$map$mapping)
    colnames(fSOM2_clustering) <- c("Cluster", "Value")
    fSOM2_clustering$Meta <- as.numeric(fSOM2$metaclustering[fSOM2_clustering[,1]])
    # attach above data to manual gated data to compare
    fsom_clustering_data$FlowSOM_ids <- fSOM2_clustering$Cluster
    
    # visualize manual gating over metaclusters
    PlotStars(fSOM2[[1]], backgroundValues = as.factor(fSOM2[[2]]))
    PlotPies(fSOM2$FlowSOM, cellTypes=cluster_comparisons$og_cid, backgroundValues = as.factor(metaClustering2))
    PlotNumbers(fSOM2$FlowSOM, backgroundValues = as.factor(metaClustering2))
    PlotMarker(fSOM2[[1]], "TotalTau")
    
    # calculate heatmap from all samples
    fsom_clustering_data %>%
      dplyr::group_by(labels) %>%
      dplyr::summarize_if(is.numeric, funs(median)) %>%
      ungroup() ->
      heat_f
    rownames(heat_f) <- heat_f$labels
      #heat[heat>1] <- 1
    heatmaply(heat_f[,panel])
    heatmaply(heat_f[,panel][-c(1:4,44:45)])
    heatmaply(heat_f[,clustering_markers])

      
      ### reassign metaclusters based upon plot inspection ###
        # neurons
        data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(11,1,2,12,22,32,31,21,71,72,62,61,51,52,42,41,6,16,5,4,93,83,30,25,3,23,15,17,29,26,24,14,13))) 
        fsom_clustering_data[data_indexes, "FlowSOM_ids"] = 109
        # non-immune glia
        data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(28,27,37,38,48,47,46,56,57,36,55,45,35,44,54,43,53,65,75,76,66,67,77,78,68,58,40,39,50,49,94,84,74,69,74,69,73,63,85,60,59,64,34,33)))
        fsom_clustering_data[data_indexes, "FlowSOM_ids"] = 110 
        # microglia
        data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(100,90,99,89,80,79,70,69,98,88,97,87,96,86,95)))
        fsom_clustering_data[data_indexes, "FlowSOM_ids"] = 107 
        # endothelial
        data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(19,20,10,9,8,7,91,92,81,82,18)))
        fsom_clustering_data[data_indexes, "FlowSOM_ids"] = 108 
        # unknown
        data_indexes = which(with(fsom_clustering_data, !FlowSOM_ids %in% c(109,110,107,108)))
        check <- fsom_clustering_data[data_indexes, "FlowSOM_ids"]
        fsom_clustering_data[data_indexes, "FlowSOM_ids"] = 200
        
        
      ### create new data.frame 'fully_FlowSOM' with combined newly labeled cells and ez_objects ###
        non_flowSOM_data <- anti_join(gated_data, fsom_clustering_data, by="spatial_id")
        non_flowSOM_data$FlowSOM_ids <- 200
        fully_FlowSOM <- rbind(fsom_clustering_data, non_flowSOM_data)
        
        # fix ez_labeling for FlowSOM_ids (not technically a FlowSOM_id)
        data_indexes = which(with(fully_FlowSOM, cluster_id == 1))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 1 # tau tangles
        data_indexes = which(with(fully_FlowSOM, cluster_id == 2))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 2 # amyloid plaques
        data_indexes = which(with(fully_FlowSOM, cluster_id == 3))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 3 # microglia-processes
        data_indexes = which(with(fully_FlowSOM, cluster_id == 4))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 4 # vessels
        
        # gate out FlowSOM_ids of endothelial cells at the edge of tissue (near gold)
        data_indexes = which(with(fully_FlowSOM, (FlowSOM_ids == 108) & run_type_id == 'Ctrl1' & point_id %in% c(22,23,24,25,26,27,55,56,85,86))) 
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 200
        data_indexes = which(with(fully_FlowSOM, (FlowSOM_ids == 108) & run_type_id == 'Ctrl2' & point_id %in% c(15,16,45,46,75,76,105)))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 200
        data_indexes = which(with(fully_FlowSOM, (FlowSOM_ids == 108) & run_type_id == 'Ctrl3' & point_id %in% c(1,50)))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 200
        data_indexes = which(with(fully_FlowSOM, (FlowSOM_ids == 108) & run_type_id == 'HiAD' & point_id %in% c(1,28,29,168,169,196,195)))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 200
        data_indexes = which(with(fully_FlowSOM, (FlowSOM_ids == 108) & run_type_id == 'MedAD' & point_id %in% c(156-1,157-1,182-1,183-1,184-1,208-1,207-1,209-1,210-1,211-1,234-1,233-1,232-1,235-1,236-1,237-1,260-1,259-1,258-1,257-1,1-1,2-1,26-1,27-1,25-1)))
        fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 200
        
        # assign any mglia or endothelial cells near projections (107) / vessels (108) respectively to null
        for (cell in seq_along(fully_FlowSOM)) {
          if (fully_FlowSOM[cell,]$FlowSOM_ids == 107) {
            x = fully_FlowSOM[cell,]$x_centroid
            y = fully_FlowSOM[cell,]$y_centroid
            #filter projections by -10 & + 10 in x and y direction (3 for mglia-processes, 4 for vessels)
            data_indexes = which(with(fully_FlowSOM, (FlowSOM_ids == 3 & (x_centroid < abs(x-MajorAxisLength/2) | x_centroid < abs(x+MajorAxisLength/2) & (y_centroid < abs(x-MajorAxisLength/2) | y_centroid < abs(x+MajorAxisLength/2))))))
            fully_FlowSOM[data_indexes, "FlowSOM_ids"] = 200
          }
        }
        
        filter(fully_FlowSOM, obj_type_id == 'cells') %>% select(FlowSOM_ids) %>% table() # check FlowSOM_ids
    
    # fSOMcodes = fSOM$FlowSOM$map$codes
    # write.csv (fSOMcodes, paste0(outputFolder,'/fSOM_FlowSOM_map_codes.csv'))
    # fSOMcodesMean = apply(fSOMcodes,2,mean)
    # write.csv(fSOMcodesMean, paste0(outputFolder,'/mean_of_fSOM_FlowSOM_map_codes.csv'))
    
    # VIZ ABOVE
    plot_data = data.frame(table(fSOM2_clustering$Meta))
    plot_ly(x = plot_data$Var1, y = plot_data$Freq, type = 'bar')
    
    # viz manual gating to manual metaclustering
    master_data_all_cuberoot %>%
           filter(cluster_id %in% c(7,8,9,10)) %>%
           ggplot() +
           aes(x=cluster_id) + 
           geom_histogram() +
           theme_minimal()
  
  #### Felix's code ####
    # make a flowframe for FlowSOM
    ff_1 <- flowFrame(exprs = data.matrix(fsom_clustering_data[,panel]), desc = list(FIL = 1))
    ff_1
    # run FlowSOM (with set.seed for reproducibility)
    set.seed(123)
    out_fSOM <- FlowSOM::ReadInput(ff_1, transform = FALSE, scale = FALSE, compensate = FALSE)
    out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = panel[-c(1:4,44:45)], xdim=10, ydim=10)
    out_fSOM_MST <- FlowSOM::BuildMST(out_fSOM)
    labels <- out_fSOM$map$mapping[,1]
    # save labels
    #write_rds(labels, “labels_sub.rds”)
    #labels <- read_rds(“labels_sub.rds”)
    fsom_clustering_data$labels <- as.factor(labels)
    # calculate heatmap from all samples
    fsom_clustering_data %>%
      dplyr::group_by(labels) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) %>%
      ungroup() ->
      heat
    rownames(heat) <- heat$labels
    heat[heat>1] <- 1
    heatmaply(heat[,panel])
    heatmaply(heat[,panel][-c(1:4,44:45)])
    heatmaply(heat[,clustering_markers])
    
  #### testing with umap/tsne ####
      obj_umap <- fsom_clustering_data[panel][-c(1:4,44:45)]
      obj_umap <- as.matrix(obj_umap)
      out_obj_umap <- umap(obj_umap)
      
      obj_umap_plot <- as.data.frame(out_obj_umap$layout)
      colnames(obj_umap_plot) <- c("UMAP1", "UMAP2")
      obj_umap_plot_all_data <- obj_umap_plot
      obj_umap_plot_all_data <- cbind(obj_umap_plot, obj_normalized)
      
      p1_U <- ggplot(obj_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = ApoE4)) +
        geom_point(size = 1) + 
        coord_fixed(ratio = 1)
      p1_U
      
      # UMAP analysis - SUBSETTED DATA
      obj_sub_umap <- subsetted_norm_data[, c(panel)]
      obj_sub_umap <- as.matrix(obj_sub_umap)
      out_obj_sub_umap <- umap(obj_sub_umap)
      
      obj_sub_umap_plot <- as.data.frame(out_obj_sub_umap$layout)
      colnames(obj_sub_umap_plot) <- c("UMAP1", "UMAP2")
      obj_sub_umap_plot_all_data <- obj_sub_umap_plot
      obj_sub_umap_plot_all_data <- cbind(obj_sub_umap_plot, subsetted_norm_data)
      
      p1_sub_U <- ggplot(obj_sub_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = CD47)) +
        geom_point(size = 1) + 
        coord_fixed(ratio = 1)
      
      p1_sub_U
      
      # viz umap
      obj_umap_plot_all_data %>%
        +     filter(UMAP1 >= -21L & UMAP1 <= 23L) %>%
        +     filter(UMAP2 >= -20L & UMAP2 <= 25L) %>%
        +     ggplot() +
        +     aes(x = UMAP1, y = UMAP2, colour = labels) +
        +     geom_point(size = 1L, show.legend=FALSE) +
        +     scale_color_hue() +
        +     theme_minimal()
      
      #tSNE
        # tSNE analysis
        # prepare object data - pick a transformed data.frame and revert back to matrix for RtSNE
        obj_rtsne <- fsom_clustering_data[panel][-c(1:4,44:45)]
        obj_rtsne <- as.matrix(obj_rtsne)
        head(obj_rtsne)
        colnames(obj_rtsne)
        dim(obj_rtsne)
        # run RtSNE
        set.seed(123)
        out_obj_rtsne <- Rtsne(obj_rtsne, dims = 2, perplexity = 50, theta = 0.5, #Run Rtnse.multicore if using Ubuntu/S3IT
                               max_iter = 1000, verbose = T, pca = F, check_duplicates=F)
        
        # prepare for plotting (double check assignments are correct between categories and objects)
        obj_tsne_plot <- as.data.frame(out_obj_rtsne$Y)
        colnames(obj_tsne_plot) <- c("tSNE1", "tSNE2")
        obj_tsne_plot_all_data <- obj_tsne_plot
        obj_tsne_plot_all_data <- cbind(obj_tsne_plot, fsom_clustering_data)
        
        # plot tSNE
        p1 <- ggplot(obj_tsne_plot_all_data, aes(x = tSNE1, y = tSNE2, color = obj_type_id)) +
          geom_point(size = 1) + 
          coord_fixed(ratio = 1)
        p1

        
#############################________################################

#############################________################################
#### Write Files for Spatial Analysis Pipelines ####
  # output df to csvs for matlab use - CytoMAP
  cytomap_data <- fully_FlowSOM
        
  #converting levels to num for csv
  levels(cytomap_data$obj_type_id) <- c('cells', 'plaques', 'tangles', 'microglia-processes', 'vessels')
  
  levels(cytomap_data$Subregions) <- c("Alveus", "Au", "CA1", "CA2", "CA3", "CA4", "DG", "EC", "Fimbria", "ML", "PreSubiculum", "SRL", "Subiculum", "TH-LV")
  subr <- levels(cytomap_data$Subregions) 
  
  levels(cytomap_data$run_type_id) <- c('HiAD', 'MedAD', 'Ctrl', 'Ctrl', 'Ctrl')
  run_levels <- levels(cytomap_data$run_type_id)
  
  cytomap_data_num <- mutate_if(cytomap_data, is.factor, ~ as.numeric(.x))
  
  cluster_names <- c('tangles', 'plaques', 'microglia-processes', 'vessels', 'oligodendrocytes', 'astrocytes', 'microglia', 'endothelial', 'neurons', 'non_immune_neuroglia')
  
  by_subregion <- FALSE
  
  for (level in 1:length(run_levels)) {
    path_name <- paste0("/Volumes/BryJC_Stanford/paper1_analysis/Fig6/single_object_pipeline/cytomap_files/", run_levels[level], "_FSOM3")
    dir.create(path_name, showWarnings = FALSE)
    for (cluster_num in c(1:4,107:110)) {
      if (by_subregion == TRUE) {
        for (sr in 1:14) {
          data_to_write <- cytomap_data_num %>% filter(run_type_id == level & FlowSOM_ids == cluster_num & Subregions == sr) %>% na.omit()
          if (cluster_num > 100) {cluster_num <- cluster_num-100}
          if (!empty(data_to_write)) {
            readr::write_csv(data_to_write, paste0(path_name,'/', 'final2subr_', subr[sr], '_', cluster_names[cluster_num], '.csv'))
          }
        }
      }
      else {
        data_to_write <- cytomap_data_num %>% filter(run_type_id == level & FlowSOM_ids == cluster_num) %>% na.omit()
        if (cluster_num > 100) {cluster_num <- cluster_num-100}
        readr::write_csv(data_to_write, paste0(path_name,'/', cluster_names[cluster_num], '.csv'))
      }
    }
  }
  
  # output df to txts for python use - SVCA
  cluster_names <- c('tau-pathies', 'amyloid-pathies', 'microglia-processes', 'vessels', 'oligodendrocytes', 'astrocytes', 'microglia', 'endothelial', 'neurons')
  use_clusters <- c(1, 2, 3, 5, 6, 7, 8, 9)
  SVCA_data <- filter(gated_data, cluster_id == 1 | cluster_id == 2 | cluster_id == 3 | cluster_id == 5 | cluster_id == 6 | cluster_id == 7 | cluster_id == 8 | cluster_id == 9)
  levels(SVCA_data$run_type_id) <- c('HiAD', 'MedAD', 'Ctrl', 'Ctrl', 'Ctrl')
  run_levels <- levels(SVCA_data$run_type_id)
  
  for (level in run_levels) {
    path_name <- paste0("/Volumes/BryJC_Stanford/paper1_analysis/Fig6/single_object_pipeline/svca_files/", level)
    run_subregions <- levels(filter(SVCA_data, run_type_id == 'HiAD')$Subregions)
    for (subr in run_subregions) {
      path_name_subr <- paste0(path_name, '_', subr)
      dir.create(path_name_subr, showWarnings = FALSE)
      # make expression table
      data_to_write <- SVCA_data %>% filter(run_type_id == level)
      panel_data <- data_to_write[panel][5:43]
      write.table(panel_data, paste0(path_name_subr,'/expression.txt'), sep = " ", row.names = FALSE, col.names = TRUE)
      # make position table
      data_to_write <- SVCA_data %>% filter(run_type_id == level)
      x_y <- data_to_write[,c(64,65)]
      write.table(x_y, paste0(path_name_subr,'/positions.txt'), sep = ",", row.names = FALSE, col.names = FALSE)
    }
  }

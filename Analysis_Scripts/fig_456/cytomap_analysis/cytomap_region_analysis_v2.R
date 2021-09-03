#### Region Analysis ####
  #marker subsets
  just_markers <- panel[-c(1:4,44:45)]
  disease_struct <- c('Amyloidbeta140','Amyloidbeta142','PHF1Tau', 'X8OHGuano', 'pTDP43', 'PolyubiK48', 'PolyubiK63')
  glia_struct <- c('Iba1', 'CD45', 'ApoE4','PanApoE2E3E4','pTDP43','CD33Lyo', 'PHF1Tau','Amyloidbeta140','Amyloidbeta142', 'PanAmyloidbeta1724','Synaptophysin','PSD95','CD47', 'VGLUT1','VGLUT2')
  glia_struct2 <- c('Iba1', 'CD45', 'ApoE4', 'CD33Lyo', 'CD47', 'EEA1', 'pTDP43', 'PHF1Tau', 'Amyloidbeta140','Amyloidbeta142','Synaptophysin','PSD95')
  neuron_struct <- c('Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','TH','SERT','PSD95','MAP2','MFN2')
  neuron_struct2 <- c('Parvalbumin','PanGAD6567','VGAT','VGLUT1','VGLUT2','PSD95','Synaptophysin','MFN2','CD47','EEA1')
  other_struct <- c('Amyloidbeta140','Amyloidbeta142','PHF1Tau','X8OHGuano','CD47','Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','TH','SERT','VGLUT1','VGLUT2','PSD95','Synaptophysin','MAP2','MFN2')
  morpho_struct <- c("Area","Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Perimeter")
  
#Data carpentry
    # prepare quantile normalized neighborhood data (note: neighborhood counts are summed from each object within each neighborhood)
    master_regions_c <- master_regions
    
    master_regions_c$run_type_id <- as.factor(master_regions_c$run_type_id)
    levels(master_regions_c$run_type_id) <- c('HiAD', 'MedAD', 'Ctrl')
    #levels(master_regions_c$run_type_id) <- relevel(master_regions_c$run_type_id, 'Ctrl','MedAD','HiAD') #doesn't work :(
    
    master_regions_c$FlowSOM_ids <- as.factor(master_regions_c$FlowSOM_ids)
    levels(master_regions_c$FlowSOM_ids) <- c('tangles','plaques','microglia','endothelial','neurons','non-immune glia')
    
    master_regions_c$OtherModel4 <- as.factor(master_regions_c$OtherModel4)
    nu_panel <- c(panel,"Area","Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Perimeter")
    
    #convert gate columns to single column
    just_gates <- master_regions[c('spatial_id',"DG","CA4","CA3","CA2","CA1")]
    just_gates %>%
      gather(gate_tag, val, -spatial_id) %>%   # Gather color cols from wide to long format
      filter(!val==0) %>%    # Drop rows with 0 values
      select(-val) -> j_gates             # Remove the unnecessary `val` column, assign to j_gates
    # attach gate_info to master_regions, renaming columns
    master_gated_regions <- inner_join(master_regions_c, j_gates, by='spatial_id')
    
#Heatmaps
    #calculate heatmap
    master_regions_c %>%
      dplyr::group_by(OtherModel4) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) %>%
      ungroup() ->
      heat_region
    
    #organize naming, ordering schema  
      rownames(heat_region) <- paste0(heat_region$OtherModel4)
      #heat[heat>1] <- 1
      heatmaply(heatmaply::normalize(heat_region[c(1:6),morpho_struct]), grid_gap = 1, Rowv = FALSE, labRow = rownames(heat_region)[c(1:6)]) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
      
    #heatmaply(normalize(heat[1:6,nu_panel]))
    heatmaply(normalize(heat_region[c(2:7),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # Ctrl
    heatmaply(normalize(heat_region[c(9:14),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # MedAD
    heatmaply(normalize(heat_region[c(16:21),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # HiAD
    #heatmaply(heat[,clustering_markers])
  
  # Full distributions --> expr data distributions for object types, by region
    panel_neighbors <- names(master_regions)
    #filtered_data <- filter(master_regions_cNU, FlowSOM_ids == 'microglia' & OtherModel4 %in% c(1,2,3,4,5,6) & CD47 <1)
    filtered_data <- filter(master_gated_regions, FlowSOM_ids == 'neurons' & OtherModel4 %in% c(1,4,5,6))
    gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(master_gated_regions), nu_panel)))
    #gathered_data$expr <- gathered_data$expr+1 # for log+1 transformations
    gathered_data_disease_struct <- filter(gathered_data, markers %in% disease_struct)
    gathered_data_glia_struct <- filter(gathered_data, markers %in% glia_struct2)
    gathered_data_neuron_struct <- filter(gathered_data, markers %in% neuron_struct2)
    gathered_data_other_struct <- filter(gathered_data, markers %in% other_struct)
    
    #violin distribtiuon -> by marker, region
    ggplot(gathered_data_glia_struct) +
      aes(x = run_type_id, y = expr, color = run_type_id) +
      geom_violin(adjust = 1L, scale = "area") +
      #scale_y_continuous(trans = "log") +
      scale_fill_gradient() +
      theme_minimal() +
      facet_wrap(vars(markers, OtherModel4))
    
    #violin distro - region grouping
    ggplot(gathered_data_glia_struct) +
      aes(x = OtherModel4, y = expr, fill = run_type_id) +
      geom_violin(adjust = 1L, scale = "area") +
      #scale_y_continuous(trans = "log") +
      scale_fill_hue() +
      labs(title = "", subtitle = "") +
      theme_minimal() +
      labs(fill = "Sample Type") +
      labs(x = "Region", y = "Expression", title = "") +
      theme(plot.title = element_text(size=22),
            axis.text = element_text(size=14),
            axis.title = element_text(size=14,face="bold"),
            strip.text = element_text(size=14),
            axis.text.x = element_text(face="bold", size=14), 
            axis.text.y = element_text(face="bold", size=14),
            legend.title = element_text(size=12, face="bold"),
            legend.text = element_text(size=12, face="bold")) +
      facet_wrap(vars(markers), scales = "free_y")
    
    #biaxial plot for deep dive into marker / object combos
    gathered_data_neuron_struct %>%
      filter(markers %in% c('PHF1Tau','CD47','VGLUT1')) %>%
      
      ggplot() +
      aes(x = markers, y = markers, fill = run_type_id) +
      geom_point(size = 1L ) +
      geom_density2d(colour = "#FF7F00") +
      scale_colour_gradientn(colours=rainbow(4)) +
      theme_minimal() +
      facet_wrap(vars(run_type_id))
    
    #plot - zScore of mean signal / region/sample // mean signal/region --> for select markers in object subset
        #calculate heatmap
        filtered_data %>%
          dplyr::group_by(run_type_id, OtherModel4) %>%
          dplyr::summarize_if(is.numeric, funs(mean)) %>%
          ungroup() ->
          filtered_region
        #calculate means of regions
        filtered_region %>% 
          dplyr::group_by(OtherModel4) %>%
          dplyr::summarize_if(is.numeric, funs(mean)) -> filtered_mean
        #calculate mean difference
        scored_means <- c()
        for (i in c(1,2,4,5,6)) {
          filt_region_select <- filter(filtered_region, OtherModel4 == i)
          filt_mean_select <- apply(filt_region_select[nu_panel], 2, function(x) mean(x)) #grab mean of means for region
          filt_region_select[nu_panel] <- t(apply(filt_region_select[nu_panel], 1, function(x) x - filt_mean_select[nu_panel])) #normalize each samples mean over region mean for each marker
          scored_means <-rbind(scored_means, filt_region_select)
        }
        gathered_scored_data <- gather(scored_means, key = "markers", value = "expr", -c('run_type_id', 'OtherModel4'))
        #gathered_data$expr <- gathered_data$expr+1 # for log+1 transformations
        gathered_data_disease_scored <- filter(gathered_scored_data, markers %in% disease_struct)
        gathered_data_glia_scored <- filter(gathered_scored_data, markers %in% glia_struct)
        #plot differences from mean
        ggplot(gathered_data_glia_scored) +
          aes(x = OtherModel4, fill = run_type_id, weight = expr) +
          geom_bar(position = "dodge") +
          scale_fill_hue() +
          theme_minimal() +
          facet_wrap(vars(markers), scales = "free_y")
        
    #plot - zScore of mean signal / region/sample // mean signal/region --> for select markers in object subset
    
    #corr plot - looking at specific marker correlations in subsetted markers, objects, regions
    for (i in c('Ctrl','MedAD','HiAD')) {
      filtered_data_corr <- filter(master_gated_regions, run_type_id == i & OtherModel4 == 5  & FlowSOM_ids == 'neurons')
      m <- cor(filtered_data_corr[neuron_struct2])
      corrplot(m, type="upper", tl.col="black", tl.srt=45, tl.cex = 0.5, insig = "blank", title = paste0(i,5))
    }
        #full corr differences
        for (i in c(1,5)) {
          filtered_data_corrA <- filter(master_regions_c, run_type_id == 'Ctrl' & OtherModel4 == i  & FlowSOM_ids == 'neurons')
          mA <- cor(filtered_data_corrA[glia_struct])
          filtered_data_corrB <- filter(master_regions_c, run_type_id == 'MedAD' & OtherModel4 == i  & FlowSOM_ids == 'neurons')
          mB <- cor(filtered_data_corrB[glia_struct])
          filtered_data_corrC <- filter(master_regions_c, run_type_id == 'HiAD' & OtherModel4 == i  & FlowSOM_ids == 'neurons')
          mC <- cor(filtered_data_corrC[glia_struct])
          mD <- mC - mA
          corrplot(mC, type="upper", tl.col="black", tl.srt=45, tl.cex = 0.5, insig = "blank")
        }
    
    #boxplot
      #cell/obj distribution of each region as percentage of total cells in region
      ggplot(master_regions_c) +
        aes(x = run_type_id, fill = FlowSOM_ids) +
        geom_bar(position = "fill") +
        scale_fill_hue() +
        theme_minimal() +
        labs(fill = "Cell / Object Type") +
        scale_x_discrete(limits=c("HiAD", "MedAD", "Ctrl")) +
        labs(x = "Sample Type", y = "Percent of total cells + objects", title = "") +
        theme(plot.title = element_text(size=22),
              axis.text = element_text(size=14),
              axis.title = element_text(size=14,face="bold"),
              strip.text = element_text(size=14),
              axis.text.x = element_text(face="bold", size=14), 
              axis.text.y = element_text(face="bold", size=14),
              legend.title = element_text(size=12, face="bold"),
              legend.text = element_text(size=12, face="bold")) +
        facet_wrap(vars(OtherModel4), scales = "free")
      #region distribution of each sample as percentage of total cells/obj
      ggplot(master_regions_c) +
        aes(x = run_type_id, fill = OtherModel4) +
        geom_bar(position = "fill") +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        labs(fill = "Region") +
        scale_x_discrete(limits=c("HiAD", "MedAD", "Ctrl")) +
        labs(x = "Sample Type", y = "Percent of total cells + objects", title = "") +
        theme(plot.title = element_text(size=22),
              axis.text = element_text(size=14),
              axis.title = element_text(size=14,face="bold"),
              strip.text = element_text(size=14),
              legend.title = element_text(size=12, face="bold"),
              legend.text = element_text(size=12, face="bold"))

      
# For Microglia - specific analysis
  filt_data_img <- dplyr::semi_join(fully_FlowSOM, filtered_data, by='spatial_id') #used to make masks in object_image_writer script
  filt_data_img <- dplyr::left_join(filt_data_img, filtered_data[c('spatial_id', 'OtherModel4')], by='spatial_id') #include region_ids  
  filt_data_img$region_cluster_id <- filt_data_img$OtherModel4 #rename region_ids
  
  # after making masks of data in matlab using getMaskedData.m, pull files into colocr for pixel co-localization / correlation
  # need to run these lines of code for each sample, each region, all channels you intend to look at it.
  tif_paths <- c()
  tif_samples <- c('HiAD','MedAD','Ctrl')
  tif_regions <- c(1,2,4,5,6)
  for (i in tif_samples) {
    for (j in tif_region) {
      #read in images
      curr_path <- tif_paths[j]
      img_names <- list.files(path = curr_path, recursive = T, full.names = T, pattern = paste0("_", j, ".tif"))
      img <- colocr::image_load(img_names)
      # calculate co-localization statistics
      img %>%
        roi_select(threshold = 0) %>%
        roi_show() #%>%
        #roi_test(type = 'both')
    }
  }
  obj_data_raw <- lapply(csv_names, read.csv) %>% bind_rows() # grab data from csv's then convert data to data.frame
  
  
# For DG->CA1 Neurons - specific analysis
  hi_tau_neurons <- filter(master_gated_regions, PHF1Tau > 0)
  lo_tau_neurons <- filter(master_gated_regions, PHF1Tau == 0)
  
  #calculate heatmap
    master_gated_regions %>%
      filter(FlowSOM_ids %in% c('neurons')) %>%
      dplyr::group_by(run_type_id, gate_tag) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) %>%
      ungroup() ->
      heat_gated_regions
    
    #organize naming, ordering schema  
    rownames(heat_gated_regions) <- paste0(heat_gated_regions$run_type_id, ' ', heat_gated_regions$gate_tag)
    #normalize to DG of each sample
    DG_vector <- apply(heat_gated_regions[c(5,10,15),other_struct], 1, function(x) x)
    heat_gr_HiAD <- heat_gated_regions[c(1:5),other_struct]
    heat_gr_MedAD <- heat_gated_regions[c(6:10),other_struct]
    heat_gr_Ctrl <- heat_gated_regions[c(11:15),other_struct]
    
    #divergent color palette
    div_color = ggplot2::scale_fill_gradient2(
      low = "blue", 
      high = "red", 
      midpoint = 0, 
      limits = c(-2, 2)
    )
    
    #full heatmap
    heatmaply(heat_gated_regions[c(1:15),other_struct], scale='column', grid_gap = 1, Rowv = TRUE, Colv = TRUE, labRow = rownames(heat_gated_regions)[c(1:15)], row_side_colors = c(rep('HiAD',5), rep('MedAD', 5), rep('Ctrl',5)))
    
    #by sample heatmap
    h1 <- heatmaply(heat_gr_HiAD, limit = c(-2,2), scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_gated_regions)[c(1:5)], scale_fill_gradient_fun = div_color) # Ctrl
    h2 <- heatmaply(heat_gr_MedAD, limit = c(-2,2), scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_gated_regions)[c(6:10)], scale_fill_gradient_fun = div_color) # MedAD
    h3 <- heatmaply(heat_gr_Ctrl, limit = c(-2,2), scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_gated_regions)[c(11:15)], scale_fill_gradient_fun = div_color) # HiAD
    #plot together
    plotly::subplot(h1,h2,h3)
  
  #biaxial plot
  h4 <- hi_tau_neurons %>%
    ggplot() +
    aes(x = CD47, y = VGLUT1, color = run_type_id) +
    geom_point(size = 1L ) +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal() +
    facet_wrap(vars(run_type_id))
  
  h5 <- lo_tau_neurons %>%
    ggplot() +
    aes(x = CD47, y = VGLUT1, color = run_type_id) +
    geom_point(size = 1L ) +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal() +
    facet_wrap(vars(run_type_id))
  
  plotly::subplot(h4,h5)
  
# For Regional Cell types (Neurons, Microglia) - specific analysis
  hi_tau_neurons <- filter(master_regions_c, FlowSOM_ids == 'neurons' & PHF1Tau > 0)
  lo_tau_neurons <- filter(master_regions_c, FlowSOM_ids == 'neurons' & PHF1Tau == 0)
  hi_tau_microglia <- filter(master_regions_c, FlowSOM_ids == 'microglia' & PHF1Tau > 0)
  lo_tau_microglia <- filter(master_regions_c, FlowSOM_ids == 'microglia' & PHF1Tau == 0)
  
  #calculate heatmap
  master_gated_regions %>%
    filter(FlowSOM_ids %in% c('microglia')) %>%
    dplyr::group_by(run_type_id, OtherModel4) %>%
    dplyr::summarize_if(is.numeric, funs(mean)) %>%
    ungroup() ->
    heat_regions
  
  #organize naming, ordering schema  
  rownames(heat_regions) <- paste0(heat_regions$run_type_id, ' ', heat_regions$OtherModel4)
  #normalize to DG of each sample
  heat_r_HiAD <- heat_regions[c(1:6),morpho_struct]
  heat_r_MedAD <- heat_regions[c(7:12),morpho_struct]
  heat_r_Ctrl <- heat_regions[c(13:18),morpho_struct]
  
  #divergent color palette
  div_color = ggplot2::scale_fill_gradient2(
    low = "blue", 
    high = "red", 
    midpoint = 0, 
    limits = c(-2.5, 2.5)
  )
  
  #full heatmap
  heatmaply(heat_regions[c(1:15),other_struct], scale='column', grid_gap = 1, Rowv = TRUE, Colv = TRUE, labRow = rownames(heat_regions)[c(1:15)], row_side_colors = c(rep('HiAD',5), rep('MedAD', 5), rep('Ctrl',5)))
  
  #by sample heatmap
  h1r <- heatmaply(heat_r_HiAD, limit = c(-2.5,2.5), scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_regions)[c(1:6)], scale_fill_gradient_fun = div_color) # Ctrl
  h2r <- heatmaply(heat_r_MedAD, limit = c(-2.5,2.5), scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_regions)[c(7:12)], scale_fill_gradient_fun = div_color) # MedAD
  h3r <- heatmaply(heat_r_Ctrl, limit = c(-2.5,2.5), scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_regions)[c(13:18)], scale_fill_gradient_fun = div_color) # HiAD
  #plot together
  plotly::subplot(h1r,h2r,h3r)
  
  #biaxial plot
  h4 <- hi_tau_neurons %>%
    ggplot() +
    aes(x = CD47, y = VGLUT1, color = run_type_id) +
    geom_point(size = 1L ) +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal() +
    facet_wrap(vars(run_type_id))
  
  h5 <- lo_tau_neurons %>%
    ggplot() +
    aes(x = CD47, y = VGLUT1, color = run_type_id) +
    geom_point(size = 1L ) +
    geom_density2d(colour = "#FF7F00") +
    theme_minimal() +
    facet_wrap(vars(run_type_id))
  
  plotly::subplot(h4,h5)

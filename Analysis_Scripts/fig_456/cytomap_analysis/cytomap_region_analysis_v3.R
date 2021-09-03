#### Region Analysis of CytoMAMP Per-Cell Data ####
# Includes anatomical region labeling and de novo region clustering labeling

  #### Marker subsets ####
    just_markers <- panel[-c(1:4,44:45)]
    disease_struct <- c('Amyloidbeta140','Amyloidbeta142',"PanAmyloidbeta1724",'PHF1Tau', 'X8OHGuano', 'pTDP43', 'PolyubiK48', 'PolyubiK63')
    glia_struct <- c('Iba1', 'CD45', 'ApoE4','PanApoE2E3E4','pTDP43','CD33Lyo', 'PHF1Tau','Amyloidbeta140','Amyloidbeta142', 'PanAmyloidbeta1724','Synaptophysin','PSD95','CD47','VGLUT1','VGLUT2')
    glia_struct2 <- c('Iba1', 'CD45', 'ApoE4', 'CD33Lyo', 'CD47', 'EEA1', 'MFN2', 'pTDP43', 'PHF1Tau', 'Amyloidbeta140','Amyloidbeta142','Synaptophysin','PSD95')
    glia_struct3 <- c('Iba1', 'CD45', 'ApoE4', 'CD33Lyo', 'CD47')
    neuron_struct <- c('Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','TH','SERT','PSD95','MAP2','MFN2')
    neuron_struct2 <- c('Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','VGLUT1','VGLUT2','PSD95','Synaptophysin','MAP2','MFN2','CD47','EEA1','PolyubiK48', 'PolyubiK63')
    neuron_struct3 <- c('Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','VGLUT1','VGLUT2','CD47')
    neuron_struct4 <- c('Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','VGLUT1','VGLUT2','PSD95','Synaptophysin','MFN2','CD47','EEA1','TH','SERT','PolyubiK48','PolyubiK63','X8OHGuano','MAP2','CD56Lyo','pTDP43','PHF1Tau')
    other_struct <- c('Amyloidbeta140','Amyloidbeta142','PHF1Tau','X8OHGuano','CD47','Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','TH','SERT','VGLUT1','VGLUT2','MAP2','MFN2')
    other_struct2 <- c('PHF1Tau','CD47','Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','VGLUT1','VGLUT2','MAP2','Iba1', 'CD45', 'ApoE4', 'CD33Lyo', 'MBP', 'MAG', 'GFAP')
    morpho_struct <- c("Area","Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Perimeter")
    
  #### Initial Data Carpentry ####
      # prepare quantile normalized neighborhood data (note: neighborhood counts are summed from each object within each neighborhood)
      master_regions_c <- master_regions
      names(master_regions_c)[77] <- c('de_novo_regions')
      
      master_regions_c$run_type_id <- factor(master_regions_c$run_type_id, levels = c(3,2,1)) #judt use factor in the futre, filling in variables
      levels(master_regions_c$run_type_id) <- c('Ctrl', 'MedAD', 'HiAD')
      
      master_regions_c$FlowSOM_ids <- as.factor(master_regions_c$FlowSOM_ids)
      levels(master_regions_c$FlowSOM_ids) <- c('tangles','plaques','microglia','endothelial','neurons','non-immune glia')
      
      master_regions_c$obj_type_id <- as.factor(master_regions_c$obj_type_id)
      levels(master_regions_c$obj_type_id) <- c('cells', 'plaques', 'tangles')
      
      master_regions_c$de_novo_regions <- as.factor(master_regions_c$de_novo_regions)
      
      nu_panel <- c(panel,"Area","Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Perimeter")
      
      
      #convert gate columns to single column
      just_gates <- master_regions[c('spatial_id',"DG","CA4","CA3","CA2","CA1")]
      just_gates %>%
        gather(anat_id, val, -spatial_id) %>% # Gather color cols from wide to long format
        filter(!val==0) %>% # Drop rows with 0 values
        select(-val) -> j_gates # Remove the unnecessary `val` column, assign to j_gates
      
      # attach gate_info to master_regions, renaming columns
      j_gates$anat_id <- factor(j_gates$anat_id, levels = c('DG','CA4','CA3','CA2','CA1'))
      master_regions_img <- left_join(master_regions_c, j_gates, by='spatial_id')
      master_gated_anat <- inner_join(master_regions_c, j_gates, by='spatial_id')
      
      #changes for image generation
        master_regions_img$anat_id <- as.numeric(as.factor((master_regions_img$anat_id))) #convert anat_id to factor, then numeric dtype for image generation
        master_regions_img$FlowSOM_ids <- as.numeric(master_regions_img$FlowSOM_ids) #convert anat_id to factor dtype for image generation
        levels(master_regions_img$obj_type_id) <- c('cells', 'amyloid-pathies', 'tau-pathies') #rename levels for reading mapping csv file paths
        #rename Ctrl runtypes to Ctrl1, Ctrl2, Ctrl3 for reading mapping csv file paths
        img_run <- fully_FlowSOM[,c("spatial_id", "run_type_id")]
        master_regions_img$run_type_id = NULL
        master_regions_img <- inner_join(master_regions_img, img_run, by='spatial_id')
      
      
  #### Heatmaps - general template ####
      #calculate heatmap
      master_regions_c %>%
        dplyr::group_by(de_novo_regions) %>%
        dplyr::summarize_if(is.numeric, funs(mean)) %>%
        ungroup() ->
        heat_region
      
      #organize naming, ordering schema  
        rownames(heat_region) <- paste0(heat_region$de_novo_regions)
        #heat[heat>1] <- 1
        heatmaply(heatmaply::normalize(heat_region[markers]), grid_gap = 1, Rowv = FALSE, labRow = rownames(heat_region)[c(1:5)]) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
        
      #heatmaply(normalize(heat[1:6,nu_panel]))
      heatmaply(normalize(heat_region[c(2:7),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # Ctrl
      heatmaply(normalize(heat_region[c(9:14),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # MedAD
      heatmaply(normalize(heat_region[c(16:21),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # HiAD
      #heatmaply(heat[,clustering_markers])
    
  #### Targeted marker processing --> expr data distributions for object types, by region ####
      panel_neighbors <- names(master_regions)
      #filtered_data <- filter(master_regions_cNU, FlowSOM_ids == 'microglia' & OtherModel4 %in% c(1,2,3,4,5,6) & CD47 <1)
      filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia' & de_novo_regions %in% c(1,2,3,4,5))
      gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(disease_cell_overlap), nu_panel)))
      #gathered_data$expr <- gathered_data$expr+1 # for log+1 transformations
      gathered_data_disease_struct <- filter(gathered_data, markers %in% disease_struct)
      gathered_data_glia_struct <- filter(gathered_data, markers %in% glia_struct3)
      gathered_data_neuron_struct <- filter(gathered_data, markers %in% neuron_struct2)
      gathered_data_other_struct <- filter(gathered_data, markers %in% other_struct)
      gathered_data_morpho_struct <- filter(gathered_data, markers %in% morpho_struct)
      
  #### Marker distributions #### 
      
      #violin distribtiuon -> by marker, region
      ggplot(gathered_data_glia_struct) +
        aes(x = run_type_id, y = expr, color = run_type_id) +
        geom_violin(adjust = 1L, scale = "area") +
        #scale_y_continuous(trans = "log") +
        scale_fill_gradient() +
        theme_minimal() +
        facet_wrap(vars(markers, de_novo_regions))
      
      #violin distro - region grouping
      ggplot(gathered_data_glia_struct) +
        aes(x = de_novo_regions, y = expr, fill = run_type_id) +
        geom_violin(adjust = 1L, scale = "area") +
        #scale_y_continuous(trans = "log") +
        scale_fill_hue() +
        labs(title = "", subtitle = "") +
        theme_minimal() +
        labs(fill = "Sample Type") +
        labs(x = "Region", y = "Expression", title = "") +
        theme1 +
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
            dplyr::group_by(run_type_id, de_novo_regions) %>%
            dplyr::summarize_if(is.numeric, funs(mean)) %>%
            ungroup() ->
            filtered_region
          #calculate means of regions
          filtered_region %>% 
            dplyr::group_by(de_novo_regions) %>%
            dplyr::summarize_if(is.numeric, funs(mean)) -> filtered_mean
          #calculate mean difference
          scored_means <- c()
          for (i in c(1,2,4,5,6)) {
            filt_region_select <- filter(filtered_region, de_novo_regions == i)
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
            aes(x = de_novo_regions, fill = run_type_id, weight = expr) +
            geom_bar(position = "dodge") +
            scale_fill_hue() +
            theme_minimal() +
            facet_wrap(vars(markers), scales = "free_y")
          
      #plot - zScore of mean signal / region/sample // mean signal/region --> for select markers in object subset
      
      #corr plot - looking at specific marker correlations in subsetted markers, objects, regions
      for (i in c('Ctrl','MedAD','HiAD')) {
        filtered_data_corr <- filter(master_gated_anat, run_type_id == i & de_novo_regions == 5  & FlowSOM_ids == 'neurons')
        m <- cor(filtered_data_corr[neuron_struct2])
        corrplot(m, type="upper", tl.col="black", tl.srt=45, tl.cex = 0.5, insig = "blank", title = paste0(i,5))
      }
          #full corr differences
          for (i in c(1,5)) {
            filtered_data_corrA <- filter(master_regions_c, run_type_id == 'Ctrl' & de_novo_regions == i  & FlowSOM_ids == 'neurons')
            mA <- cor(filtered_data_corrA[glia_struct])
            filtered_data_corrB <- filter(master_regions_c, run_type_id == 'MedAD' & de_novo_regions == i  & FlowSOM_ids == 'neurons')
            mB <- cor(filtered_data_corrB[glia_struct])
            filtered_data_corrC <- filter(master_regions_c, run_type_id == 'HiAD' & de_novo_regions == i  & FlowSOM_ids == 'neurons')
            mC <- cor(filtered_data_corrC[glia_struct])
            mD <- mC - mA
            corrplot(mC, type="upper", tl.col="black", tl.srt=45, tl.cex = 0.5, insig = "blank")
          }
    
  #### Marking general distributions ####   
          
      #boxplot
        #cell/obj distribution of each region as percentage of total cells in region
        ggplot(master_regions_c) +
          aes(x = run_type_id, fill = FlowSOM_ids) +
          geom_bar(position = "fill") +
          scale_fill_manual(values = color1) +
          theme_minimal() +
          labs(fill = "Cell / Object Type") +
          scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
          labs(x = "Sample Type", y = "Percent of total cells and objects", title = "Hippocampus relative cell or object composition") +
          theme(plot.title = element_text(size=22, hjust = 0.5),
                axis.text = element_text(size=14),
                axis.title = element_text(size=14,face="bold"),
                strip.text = element_text(size=14),
                axis.text.x = element_text(face="bold", size=14), 
                axis.text.y = element_text(face="bold", size=14),
                legend.title = element_text(size=12, face="bold"),
                legend.text = element_text(size=12, face="bold"))
          #facet_wrap(vars(run_type_id), scales = "free")
          
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
  
        
  #### For Microglia - specific analysis ####
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
    
    
  #### For DG->CA1 Neurons - specific analysis ####
    hi_tau_neurons <- filter(master_gated_anat, PHF1Tau > 0)
    lo_tau_neurons <- filter(master_gated_anat, PHF1Tau == 0)
    
    #calculate heatmap
      master_gated_anat %>%
        filter(FlowSOM_ids %in% c('neurons')) %>%
        dplyr::group_by(run_type_id, anat_id) %>%
        dplyr::summarize_if(is.numeric, funs(mean)) %>%
        ungroup() ->
        heat_gated_regions
      
      #organize naming, ordering schema  
      rownames(heat_gated_regions) <- paste0(heat_gated_regions$run_type_id, ' ', heat_gated_regions$anat_id)
      #normalize to DG of each sample
      DG_vector <- apply(heat_gated_regions[c(5,10,15),neuron_struct2], 1, function(x) x)
      heat_gr_HiAD <- heat_gated_regions[c(1:5),neuron_struct2]
      heat_gr_MedAD <- heat_gated_regions[c(6:10),neuron_struct2]
      heat_gr_Ctrl <- heat_gated_regions[c(11:15),neuron_struct2]
      
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
      h1 <- heatmaply(heat_gr_HiAD, limit = c(-2,2), scale = 'column', Rowv = FALSE, Colv = TRUE, labRow = rownames(heat_gated_regions)[c(1:5)], scale_fill_gradient_fun = div_color) # Ctrl
      h2 <- heatmaply(heat_gr_MedAD, limit = c(-2,2), scale = 'column', Rowv = FALSE, Colv = TRUE, labRow = rownames(heat_gated_regions)[c(6:10)], scale_fill_gradient_fun = div_color) # MedAD
      h3 <- heatmaply(heat_gr_Ctrl, limit = c(-2,2), scale = 'column', Rowv = FALSE, Colv = TRUE, labRow = rownames(heat_gated_regions)[c(11:15)], scale_fill_gradient_fun = div_color) # HiAD
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
    
  #### For Regional Cell types (Neurons, Microglia) - specific analysis ####
    hi_tau_neurons <- filter(master_regions_c, FlowSOM_ids == 'neurons' & PHF1Tau > 0)
    lo_tau_neurons <- filter(master_regions_c, FlowSOM_ids == 'neurons' & PHF1Tau == 0)
    hi_tau_microglia <- filter(master_regions_c, FlowSOM_ids == 'microglia' & PHF1Tau > 0)
    lo_tau_microglia <- filter(master_regions_c, FlowSOM_ids == 'microglia' & PHF1Tau == 0)
    
    #calculate heatmap
    disease_cell_overlap %>%
      filter(FlowSOM_ids %in% c('neurons')) %>%
      dplyr::group_by(run_type_id, de_novo_regions) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) %>%
      ungroup() ->
      heat_regions
    
    #organize naming, ordering schema  
    rownames(heat_regions) <- paste0(heat_regions$run_type_id, ' ', heat_regions$de_novo_regions)
    #normalize to DG of each sample
    heat_r_HiAD <- heat_regions[c(1:5),neuron_struct2]
    heat_r_MedAD <- heat_regions[c(6:10),neuron_struct2]
    heat_r_Ctrl <- heat_regions[c(11:15),neuron_struct2]
    
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
    h1r <- heatmaply(heat_r_HiAD, limit = c(-2.5,2.5), scale = 'column', Rowv = FALSE, Colv = TRUE, labRow = rownames(heat_regions)[c(1:5)], scale_fill_gradient_fun = div_color) # Ctrl
    h2r <- heatmaply(heat_r_MedAD, limit = c(-2.5,2.5), scale = 'column', Rowv = FALSE, Colv = TRUE, labRow = rownames(heat_regions)[c(6:10)], scale_fill_gradient_fun = div_color) # MedAD
    h3r <- heatmaply(heat_r_Ctrl, limit = c(-2.5,2.5), scale = 'column', Rowv = FALSE, Colv = TRUE, labRow = rownames(heat_regions)[c(11:15)], scale_fill_gradient_fun = div_color) # HiAD
    #plot together
    plotly::subplot(h1r,h2r,h3r, nrows = 1)
    
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
  
  #### For tau_object and amyloid_object positive cells ####
    
    # Summary stats (object distribution)
    ggplot(disease_cell_overlap) +
      aes(x = run_type_id, fill = tau_obj_POS) +
      geom_bar(position = "fill") +
      scale_fill_hue() +
      theme_minimal() +
      facet_grid(rows = vars(FlowSOM_ids), cols = vars(de_novo_regions), scales = "free_y") +
      theme1
    
    filter(master_gated_anat2, FlowSOM_ids %in% c('neurons', 'microglia')) %>% 
      ggplot() +
      aes(x = run_type_id, fill = tau_obj_POS) +
      geom_bar(position = "fill") +
      scale_fill_hue() +
      theme_minimal() +
      facet_grid(rows = vars(FlowSOM_ids), cols = vars(anat_id), scales = "free_y") +
      theme1
    
    # Marker distributions
    ggplot(gathered_data_glia_struct) +
      aes(x = de_novo_regions, y = expr, fill = tau_obj_POS) +
      geom_boxplot(adjust = 1L, scale = "area") +
      #scale_y_continuous(trans = "log") +
      scale_fill_hue() +
      labs(title = "", subtitle = "") +
      theme_minimal() +
      labs(fill = "Sample Type") +
      labs(x = "Region", y = "Expression", title = "") +
      theme1 +
      facet_wrap(vars(markers), scales = "free_y")
    
    # more specific distro of markers - neurons
    filtered_data <- filter(master_gated_anat2, anat_id %in% c('CA1','CA2') & FlowSOM_ids == 'neurons' & tau_obj_POS == 0)
    gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(master_gated_anat2), nu_panel)))
    ggplot(filter(gathered_data, markers %in% setdiff(c(disease_struct, neuron_struct2), c('X8OHGuano', 'pTDP43', 'EEA1')))) +
      aes(x = run_type_id, y = expr, fill = anat_id) +
      geom_boxplot(adjust = 1L, scale = "area") +
      #scale_y_continuous(trans = "log") +
      scale_fill_hue() +
      labs(title = "", subtitle = "") +
      theme_minimal() +
      labs(fill = "Sample Type") +
      labs(x = "Region", y = "Expression", title = "") +
      theme1 +
      facet_wrap(vars(markers), scales = "free_y")
    
    #DG->CA1
    master_gated_anat2 <- inner_join(disease_cell_overlap, j_gates, by='spatial_id')
    
    # Cell type object prediction (logistic regression)
    # Check percentage of regions and tau objects - healthy region and where the objects fall vs unhealthy regions
    
  #### For cell positive tau_objects and amyloid_objects ####
    
    # Summary stats (object distribution) - tau
    filter(cell_disease_overlap, FlowSOM_ids %in% c('tangles')) -> p1t
    p1t$de_novo_regions <- factor(p1t$de_novo_regions, levels=c(1,2,5,3,4))
    filter(p1t, de_novo_regions %in% c(1,5,4)) %>%
      ggplot() +
      aes(x = de_novo_regions, fill = cells_tau_obj_POS) +
      geom_bar(position = "dodge") +
      scale_fill_hue() +
      theme_minimal(20) +
      labs(title='tau')
      #facet_grid(rows = vars(FlowSOM_ids), cols = vars(de_novo_regions), scales = "free_y")
    
    # Summary stats (object distribution) - amyloid
    filter(cell_disease_overlap, FlowSOM_ids %in% c('plaques')) -> p1a
    p1a$de_novo_regions <- factor(p1a$de_novo_regions, levels=c(1,2,5,3,4))
    filter(p1a, de_novo_regions %in% c(1,5,3)) %>%
      ggplot() +
      aes(x = de_novo_regions, fill = cells_amyloid_obj_POS) +
      geom_bar(position = "dodge") +
      scale_fill_hue() +
      theme_minimal(20) +
      labs(title='amyloid')
      #facet_grid(rows = vars(FlowSOM_ids), cols = vars(de_novo_regions), scales = "free_y")
    
    #tangles, plaques
    filter(master_gated_anat3, FlowSOM_ids %in% c('tangles', 'plaques')) %>% 
      ggplot() +
      aes(x = run_type_id, fill = cells_tau_obj_POS) +
      geom_bar(position = "fill") +
      scale_fill_hue() +
      theme_minimal(10) +
      facet_grid(rows = vars(FlowSOM_ids), cols = vars(anat_id), scales = "free_y") +
    
    # Marker distributions
    ggplot(gathered_data_glia_struct) +
      aes(x = de_novo_regions, y = expr, fill = tau_obj_POS) +
      geom_boxplot(adjust = 1L, scale = "area") +
      #scale_y_continuous(trans = "log") +
      scale_fill_hue() +
      labs(title = "", subtitle = "") +
      theme_minimal() +
      labs(fill = "Sample Type") +
      labs(x = "Region", y = "Expression", title = "") +
      theme1 +
      facet_wrap(vars(markers), scales = "free_y")
    
    #DG->CA1
    master_gated_anat3 <- inner_join(cell_disease_overlap, j_gates, by='spatial_id')
    
  #### Cells in Objects vs Objects in cells ####
    disease_cell_overlap %>%
      dplyr::group_by(run_type_id, de_novo_regions) %>%
      dplyr::summarise(ratio_t = sum(tau_obj_POS==1)/sum(tau_obj_POS==0), ratio_a = sum(amyloid_obj_POS==1)/sum(amyloid_obj_POS==0), n = n()) -> dco3
    cell_disease_overlap %>%
      dplyr::group_by(run_type_id, de_novo_regions, FlowSOM_ids) %>%
      dplyr::summarise(ratio_ct = sum(cells_tau_obj_POS==1)/sum(cells_tau_obj_POS==0), ratio_ca = sum(cells_amyloid_obj_POS==1)/sum(cells_amyloid_obj_POS==0), n_c = n()) -> cdo3
    dco3 <- left_join(dco3, filter(cdo3, FlowSOM_ids == 'tangles')[c('run_type_id', 'de_novo_regions', 'ratio_ct', 'n_c')], by=c('run_type_id', 'de_novo_regions'))
    dco3 <- left_join(dco3, filter(cdo3, FlowSOM_ids == 'plaques')[c('run_type_id', 'de_novo_regions', 'ratio_ca', 'n_c')], by=c('run_type_id', 'de_novo_regions'))
    
    ggplot(dco3, aes(x = ratio_a, y = ratio_ca)) +
      geom_point(aes(color = run_type_id), size=5, alpha=0.75) +
      scale_color_manual(values = color4) +
      #scale_x_continuous(limits = c(0,1.5), breaks = seq(0, 1.5, 0.5)) +
      #scale_y_continuous(limits = c(0,2.5), breaks = seq(0, 2.5, 0.5)) + 
      theme_bw(20) +
      labs(x="A object Overlap Ratio", y="Cell overlap A object Ratio", title="De novo regions: Ratio of cells containing A objects", color="Severity", shape="Cell type") +
      #theme1 +
      facet_wrap(vars(de_novo_regions), nrow = 1, scales = "fixed")
      
  #### Heatmaps 2.0 --> by regions, overlap, cell types ####
    disease_cell_overlap$de_novo_regions <- factor(disease_cell_overlap$de_novo_regions, levels=c(1,2,5,3,4))
    #TAU
      disease_cell_overlap %>%
        dplyr::group_by(de_novo_regions, FlowSOM_ids, tau_obj_POS) %>%
        dplyr::summarize_if(is.numeric, funs(median)) %>%
        ungroup() ->
        heat_regions_2t
      
      #organize naming, ordering schema  
      rownames(heat_regions_2t) <- paste0(heat_regions_2t$de_novo_regions, ' ', heat_regions_2t$FlowSOM_ids, ' ', heat_regions_2t$tau_obj_POS)
      #heat[heat>1] <- 1
      choose_rows <- which(heat_regions_2t$FlowSOM_ids=='neurons' & heat_regions_2t$de_novo_regions %in% c(1,5,4) & heat_regions_2t$tau_obj_POS %in% c(0,1))
      heatmaply(heat_regions_2t[choose_rows,markers], grid_gap = 1, Rowv = TRUE, labRow = rownames(heat_regions_2t)[choose_rows]) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
      
    #AMYLOID
      disease_cell_overlap %>%
        dplyr::group_by(de_novo_regions, FlowSOM_ids, amyloid_obj_POS) %>%
        dplyr::summarize_if(is.numeric, funs(median)) %>%
        ungroup() ->
        heat_regions_2a
      
      #organize naming, ordering schema  
      rownames(heat_regions_2a) <- paste0(heat_regions_2a$de_novo_regions, ' ', heat_regions_2a$FlowSOM_ids, ' ', heat_regions_2a$amyloid_obj_POS)
      #heat[heat>1] <- 1
      choose_rows <- which(heat_regions_2a$FlowSOM_ids=='neurons' & heat_regions_2a$de_novo_regions %in% c(1,5,3) & heat_regions_2a$amyloid_obj_POS %in% c(0,1))
      heatmaply(heat_regions_2a[choose_rows,neuron_struct3], grid_gap = 1, Rowv = TRUE, labRow = rownames(heat_regions_2a)[choose_rows]) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
      
  #### Heatmaps 3.0 --> by digging down into object positivity marker differences ####
    #TAU
      disease_cell_overlap %>%
      dplyr::group_by(de_novo_regions, FlowSOM_ids) %>%
      dplyr::summarise(is.numeric, funs(ratio_t = sum(tau_obj_POS==1)/sum(tau_obj_POS==0), ratio_a = sum(amyloid_obj_POS==1)/sum(amyloid_obj_POS==0))) -> ratios_hm1
      
      disease_cell_overlap %>%
      dplyr::group_by(de_novo_regions, FlowSOM_ids, tau_obj_POS) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) -> ratio_hm_t
      ratio_hm_t1 <- filter(ratio_hm_t, tau_obj_POS == 1)
      ratio_hm_t0 <- filter(ratio_hm_t, tau_obj_POS == 0)
      ratio_hm_tFinal <- ratio_hm_t1[markers] / ratio_hm_t0[markers]
      ratio_hm_tFinal <- log(ratio_hm_tFinal, base=2)
      ratio_hm_tFinal <- dplyr::bind_cols(ratio_hm_t1[c('de_novo_regions', 'FlowSOM_ids')], ratio_hm_tFinal)
    
      #organize naming, ordering schema  
      rownames(ratio_hm_tFinal) <- paste0(ratio_hm_tFinal$de_novo_regions, ' ', ratio_hm_tFinal$FlowSOM_ids, ' ', ratio_hm_tFinal$tau_obj_POS)
      #heat[heat>1] <- 1
      heatmaply(ratio_hm_tFinal[markers[-39]], grid_gap = 1, Rowv = FALSE, labRow = rownames(ratio_hm_tFinal)) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
      
      #heatmaply(normalize(heat[1:6,nu_panel]))
      heatmaply(ratio_hm_tFinal[which(heat_regions_2t$FlowSOM_ids=='neurons'),neuron_struct2], grid_gap = 1, Rowv = FALSE, labRow = rownames(ratio_hm_tFinal)[c(3,19,15)]) # Neurons
      heatmaply(ratio_hm_tFinal[which(heat_regions_2t$FlowSOM_ids=='microglia'),glia_struct3], grid_gap = 1, Rowv = FALSE, labRow = rownames(ratio_hm_tFinal)[c(1,13,17)]) # Microglia
      #heatmaply(heat[,clustering_markers])
      
    #AMYLOID
    disease_cell_overlap %>%
      dplyr::group_by(de_novo_regions, FlowSOM_ids, amyloid_obj_POS) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) %>%
      ungroup() ->
      heat_regions_2a
    
    #organize naming, ordering schema  
    rownames(heat_regions_2a) <- paste0(heat_regions_2a$de_novo_regions, ' ', heat_regions_2a$FlowSOM_ids, ' ', heat_regions_2a$amyloid_obj_POS)
    #heat[heat>1] <- 1
    heatmaply(heat_regions_2a[markers], grid_gap = 1, Rowv = TRUE, labRow = rownames(heat_regions_2a)) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
    
  #### Post-UMAP Differences ####
    ggplot(filter(disease_cell_overlap, FlowSOM_ids == 'neurons'), aes(x = PolyubiK48, y = PHF1Tau)) +
      geom_point(aes(color = tau_obj_POS), size=2, alpha=0.75) +
      #scale_color_gradient(low='blue', high='red', limits=c(0,1.5)) +
      #scale_shape_manual(values=c(15,17,19,18)) +
      scale_color_manual(values = color4) +
      #xlim(0,2.25) +
      #ylim(0,0.8) +
      theme_bw(20) +
      #labs(x="Tau object Overlap Ratio", y="Amyloid object Overlap Ratio", title="DG to CA1: Ratio of cells containing disease objects", color="Anantomical Region", shape="Cell type") +
      #theme1 +
      facet_wrap(vars(run_type_id), scales = "fixed")
    
    # Tangles -> Polyubik48 vs PHF1Tau inverse relationship
    tanglesPk48 <- filter(cell_disease_overlap, FlowSOM_ids == 'tangles' & PolyubiK48 >= 0.5)
    
    # Iba+ CD45+ ApoE4+ Microglia (also looking at CD33 and pTDP43)
    
  
  ############ Neuron only dive ############  
    #### SCORPIUS Trajectory Inference ####
    install.packages('SCORPIUS')
    library('SCORPIUS')
    scorp_data <- list()
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', de_novo_regions %in% c(5))
    set.seed(321)
    
    scorp_n_markers <- c('CD47','TotalTau','Calretinin','CD56Lyo','Parvalbumin','MAP2','PolyubiK48','MFN2',
                          'PanGAD6567','VGAT','PolyubiK63','Synaptophysin','VGLUT1','VGLUT2','Calbindin','PSD95',
                          'Presenilin1NTF',
                          'PHF1Tau','Amyloidbeta142','Amyloidbeta140','PanAmyloidbeta1724')
    
    scorp_data$expression <- as.matrix(filtered_data[,scorp_n_markers])
    scorp_data$sample_info <- as.data.frame(filtered_data[,c('run_type_id', 'de_novo_regions', 'tau_obj_POS', 'amyloid_obj_POS')])
    names(scorp_data$sample_info) <- c('run_type_id', 'de_novo_regions', 'tau_obj_POS', 'amyloid_obj_POS')
    
    scorp_expression <- scorp_data$expression
    scorp_group_name <- scorp_data$sample_info$de_novo_regions
    scorp_space <- reduce_dimensionality(scorp_expression, "euclidean", ndim = 3)
    draw_trajectory_plot(scorp_space, progression_group = scorp_group_name, contour = TRUE, progression_group_palette = color3, point_size = 0.2)
    
    scorp_traj <- infer_trajectory(scorp_space)
    # save as ggplot2 object
    test1 <- draw_trajectory_plot(
      scorp_space, 
      progression_group = scorp_group_name,
      progression_group_palette = color3,
      path = scorp_traj$path,
      contour = TRUE,
      point_size = 0.2
    )    
    ggsave(filename = 'draw_traj_path_denovoregions.pdf', plot = test1, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v2/'
    )
    
    gimp <- gene_importances(scorp_expression, scorp_traj$time, num_permutations = 0, num_threads = 8)
    gene_sel <- gimp[1:nrow(gimp),]
    expr_sel <- scorp_expression[,gene_sel$gene]
    
    gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
    gene_sel <- gimp$gene[gimp$qvalue < 1.1]
    expr_sel <- scale_quantile(scorp_expression[,gene_sel])
    
    # save as heatmap object (exclude everything after first line if unable to find good proportions)
    draw_trajectory_heatmap(expr_sel, scorp_traj$time, scorp_group_name, show_labels_row = TRUE, progression_group_palette = color3,
                            cellwidth = 0.01, cellheight = 10, fontsize = 10, width = 11, height = 4, 
                            filename = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v2/severity_SCORPIUS_v2_denovoregions.pdf')

    modules <- extract_modules(scale_quantile(expr_sel), scorp_traj$time, verbose = FALSE)
    draw_trajectory_heatmap(expr_sel, scorp_traj$time, scorp_group_name, modules, show_labels_row = TRUE, progression_group_palette = color4)
    
    #notes: redo for tangles and plaques, check groupings for neurons (severity and diseasePOS, denovoregions & diseasePOS)
      #for John
        write_csv(filtered_data, path = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v2/bjc_neurons_cytoskel.csv')
      #for John - change factors to numeric (changes diseasePOS to 1,2 instead of 0,1)
        filtered_data <- (mutate_if(filtered_data, is.factor, ~ as.numeric(as.character(.x))))
        write_csv(filtered_data, path = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v2/bjc_neurons_cytoskel_num.csv')
    
    
    #### Hier-Clustering ####
    neuroPlus_markers <- c(scorp_n_markers, 'HistoneH3Lyo', 'C12')
        
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', de_novo_regions %in% c(1))
        filtered_data[filtered_data > 0.5] <- 0.5
        pheatmap(filtered_data[, c(neuroPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', de_novo_regions %in% c(4))
        filtered_data[filtered_data > 0.5] <- 0.5
        pheatmap(filtered_data[, c(neuroPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', de_novo_regions %in% c(5))
        filtered_data[filtered_data > 0.5] <- 0.5
        pheatmap(filtered_data[, c(neuroPlus_markers)], cluster_cols = F)
        
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', run_type_id %in% c('Ctrl'))
        filtered_data[filtered_data > 0.5] <- 0.5
        pheatmap(filtered_data[, c(neuroPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', run_type_id %in% c('MedAD'))
        filtered_data[filtered_data > 0.5] <- 0.5
        pheatmap(filtered_data[, c(neuroPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'neurons', run_type_id %in% c('HiAD'))
        filtered_data[filtered_data > 0.5] <- 0.5
        pheatmap(filtered_data[, c(neuroPlus_markers)], cluster_cols = F)
        
    #### Neuron subclustering proportions ####
    ggplot() +
          aes(x=)
        
  ############ Microglia only dive ############  
    #### Hier-Clustering ####
    mgliaPlus_markers <- c(hgr_panel_m, 'HistoneH3Lyo', 'C12')
    
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia', de_novo_regions %in% c(1))
    filtered_data[filtered_data > 0.5] <- 0.5
    pheatmap(filtered_data[, c(mgliaPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia', de_novo_regions %in% c(4))
    filtered_data[filtered_data > 0.5] <- 0.5
    pheatmap(filtered_data[, c(mgliaPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia', de_novo_regions %in% c(5))
    filtered_data[filtered_data > 0.5] <- 0.5
    pheatmap(filtered_data[, c(mgliaPlus_markers)], cluster_cols = F)
    
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia', run_type_id %in% c('Ctrl'))
    filtered_data[filtered_data > 0.5] <- 0.5
    pheatmap(filtered_data[, c(mgliaPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia', run_type_id %in% c('MedAD'))
    filtered_data[filtered_data > 0.5] <- 0.5
    pheatmap(filtered_data[, c(mgliaPlus_markers)], cluster_cols = F)
    filtered_data <- filter(disease_cell_overlap, FlowSOM_ids == 'microglia', run_type_id %in% c('HiAD'))
    filtered_data[filtered_data > 0.5] <- 0.5
    pheatmap(filtered_data[, c(mgliaPlus_markers)], cluster_cols = F)
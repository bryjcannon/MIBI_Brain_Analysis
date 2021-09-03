#### Neighborhood Analysis ####
  
  #### Initial Data Carpentry ####
    # prepare quantile normalized neighborhood data (note: neighborhood counts are summed from each object within each neighborhood)
    master_neighborhood_c <- master_neighborhood
    names(master_neighborhood_c)[80] <- c('de_novo_regions')
    master_neighborhood_c$de_novo_regions <- as.factor(master_neighborhood_c$de_novo_regions)

  #### Heatmaps ####
    # prepare quantile normalized neighborhood data (note: neighborhood counts are summed from each object within each neighborhood)
    master_neighborhood_c <- master_neighborhood
    nu_panel <- c(panel,"Area","Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Perimeter")
    # normalize to total MFI for each neighborhood
    norm_sumMFI_vector <- apply(master_neighborhood_c[panel][-c(1:4,44:45)], 1, function(x) sum(x))
    master_neighborhood_c[panel][-c(1:4,44:45)] <- (master_neighborhood_c[panel][-c(1:4,44:45)]) / as.numeric(norm_sumMFI_vector)
    master_neighborhood_c[is.na(master_neighborhood_c)] = 0
    # normalize to max MFI for each neighborhood
    norm_maxMFI_vector <- apply(master_neighborhood_c[panel][-c(1:4,44:45)], 1, function(x) max(x))
    master_neighborhood_c[panel][-c(1:4,44:45)] <- (master_neighborhood_c[panel][-c(1:4,44:45)]) / as.numeric(norm_maxMFI_vector)
    master_neighborhood_c[is.na(master_neighborhood_c)] = 0
    # normalize to num of cells in each neighborhood
    master_neighborhood_c[panel][-c(1:4,44:45)] <- (master_neighborhood_c[panel][-c(1:4,44:45)]) / master_neighborhood_c$NCells
    master_neighborhood_c[is.na(master_neighborhood_c)] = 0
    # quant normalize per sample
    for (run in c('Ctrl', 'MedAD', 'HiAD')) {
      data_indexes = which(with(master_neighborhood_c, run_type_name == run))
      norm_range_vector <- apply(master_neighborhood_c[data_indexes,nu_panel], 2, function(x) quantile(x, 0.99, names = F))
      master_neighborhood_c[data_indexes,nu_panel] <- t(t(master_neighborhood_c[data_indexes,nu_panel]) / as.numeric(norm_range_vector))
    }
    # quant normalize for entire dataset
    norm_range_vector <- apply(master_neighborhood_c[nu_panel], 2, function(x) quantile(x, 0.99, names = F)) 
    master_neighborhood_c[nu_panel] <- t(t(master_neighborhood_c[nu_panel]) / as.numeric(norm_range_vector)) 
    
    #calculate heatmap
    master_neighborhood_c %>%
      filter(NCells > 0) %>%
      #dplyr::group_by(OtherModel4, run_type_name) %>%
      dplyr::group_by(run_type_name) %>%
      dplyr::summarize_if(is.numeric, funs(median), na.rm=TRUE) %>%
      ungroup() ->
      heat_neighborhood
    
    #display heatmaps
      #organize naming, ordering schema  
        #rownames(heat_neighborhood) <- paste0(heat_neighborhood$run_type_name, heat_neighborhood$OtherModel4)
        rownames(heat_neighborhood) <- heat_neighborhood$run_type_name
        # order by sample
        heatmaply(heatmaply::normalize(heat_neighborhood[c(1:3),panel][neuron_struct2]), grid_gap = 1, Rowv = FALSE, labRow = rownames(heat_neighborhood)[c(1:3)])
        # order by region
        heatmaply(heatmaply::normalize(heat_neighborhood[c(1:6),panel[-c(1:4,44:45)]]), grid_gap = 1, Rowv = FALSE, labRow = rownames(heat_neighborhood)[c(1:6)]) #row_side_colors = c(rep('1',3), rep('2', 3), rep('3',3), rep('4',3), rep('5',3), rep('6',3)
        
        #heatmaply(normalize(heat[1:6,nu_panel]))
        heatmaply((heat_neighborhood[c(2:7),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # Ctrl
        heatmaply((heat_neighborhood[c(9:14),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # MedAD
        heatmaply((heat_neighborhood[c(16:21),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # HiAD
        #heatmaply(heat[,clustering_markers])
  
  #### Full distributions --> expr data distributions for object types, by region ####
    panel_neighbors <- names(master_neighborhood)
    gathered_data <- gather(master_neighborhood_c, key = "markers", value = "expr", -c(setdiff(names(master_neighborhood), nu_panel)))
    gathered_data$expr <- gathered_data$expr+1 # for log+1 transformations
    gathered_data_disease_struct <- filter(gathered_data, markers == c('Amyloidbeta140','Amyloidbeta142','ApoE4', 'PanAmyloidbeta1724', 'PanApoE2E3E4', 'PHF1Tau', 'X8OHGuano', 'pTDP43'))
    gathered_data_glia_struct <- filter(gathered_data, markers == c('CD47', 'ApoE4', 'GFAP', 'CD33Lyo', 'Iba1'))
    gathered_data_neuron_struct <- filter(gathered_data, markers == c('Calretinin','Calbindin','Parvalbumin','PanGAD6567','VGAT','TH','SERT','PSD95','MAP2','MFN2'))
    gathered_data_other_struct <- filter(gathered_data, markers == c('PolyubiK63', 'PolyubiK48'))
    
    #violin distriubtiuon
    ggplot(gathered_data_other_struct) +
      aes(x = run_type_name, y = expr, color = run_type_name) +
      geom_violin(adjust = 1L, scale = "area") +
      scale_y_continuous(trans = "log") +
      scale_fill_gradient() +
      theme_minimal() +
      facet_wrap(vars(markers, OtherModel4))
    
    #### UMAP ####
      neigh_umap_data <- master_neighborhood_c[c('MAP2','CD31','CD105','MCT1','Iba1','CD45','MAG','MBP','GFAP')]
      neigh_umap_data <- as.matrix(neigh_umap_data)
      neigh_umap <- umap(neigh_umap_data)
      
      neigh_umap_plot <- as.data.frame(neigh_umap$layout)
      colnames(neigh_umap_plot) <- c("UMAP1", "UMAP2")
      neigh_umap_plot_all_data <- neigh_umap_plot
      neigh_umap_plot_all_data <- cbind(neigh_umap_plot, master_neighborhood_c)
      
      plot_UMAP1 <- ggplot(obj_umap_plot_all_data, aes(x = UMAP1, y = UMAP2)) +
        geom_point(size = 1) + 
        coord_fixed(ratio = 1)
      plot_UMAP1
    
  #### Linear Axis Sorting (Pseudo-space stand in) ####
    #save template
      ggsave(filename = 'PolyubiK48_Y_PHF1Tau_X_runtypename_denovoregions145.pdf', plot = t, width = 30, height = 10, dpi = 300, device = 'pdf',
             path = '/Volumes/BryJC_Stanford/paper1_analysis/Fig6/plots/fig_plots_final/pseudo_axis/'
      )
      
    # Template - line graph with tangles on X, cell pops on Y
      #
      t <- ggplot(master_neighborhood_c, aes(x = tangles)) +
             geom_line(aes(y = neurons), color = 'darkred') +
             geom_line(aes(y = microglia), color="steelblue", linetype="twodash") +
             facet_wrap(vars(de_novo_regions))
      
    # Template - line graph with tangles on X, cell pops on Y - shifts to have jitter points as well
      # broken down by de novo region across severity
      t <- ggplot(filter(master_neighborhood_c, de_novo_regions %in% c(1,2,3,4,5)), aes(x = neurons, y = PolyubiK48, color = de_novo_regions)) +
             geom_jitter(alpha = 0.01, height = 0.05) +
             geom_smooth(method = "gam", formula = y ~ s(x)) +
             scale_color_manual(values = color3) +
             facet_wrap(vars(run_type_name), scale = 'free')
      # same as above but with all on one line, by severity
      t <- ggplot(filter(master_neighborhood_c, de_novo_regions %in% c(1,4,5)), aes(x = PHF1Tau, y = PolyubiK48, color = de_novo_regions)) +
             geom_jitter(alpha = 0.01, height = 0.05) +
             geom_smooth(method = "gam", formula = y ~ s(x)) +
             scale_color_manual(values = color3[c(1,4,5)]) +
             facet_wrap(vars(run_type_name))
    
    # Template - line graph with _ on X, _ on Y - shifts to have jitter points as well
      #
      t <- ggplot(filter(master_neighborhood_c, de_novo_regions %in% c(1,4,5)), aes(x = PHF1Tau, y = PolyubiK48, color = run_type_name)) +
             geom_jitter(alpha = 0.01, height = 0.05, aes()) +
             geom_smooth(method = "gam", formula = y ~ s(x), span = 0.1) +
             scale_color_manual(values = color4) +
             facet_wrap(vars(de_novo_regions))
      
#### Neighborhood Analysis ####
  #Heatmaps
    # prepare quantile normalized neighborhood data (note: neighborhood counts are summed from each object within each neighborhood)
    master_neighborhood_c <- master_neighborhood
    nu_panel <- c(panel,"Area","Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","Orientation","Perimeter")
    # quant normalize per sample
    for (run in c('Ctrl', 'MedAD', 'HiAD')) {
      data_indexes = which(with(master_neighborhood_c, run_type_id == run))
      norm_range_vector <- apply(master_neighborhood_c[data_indexes,nu_panel], 2, function(x) quantile(x, 0.99, names = F))
      master_neighborhood_c[data_indexes,nu_panel] <- t(t(master_neighborhood_c[data_indexes,nu_panel]) / as.numeric(norm_range_vector))
    }
    # quant normalize for entire dataset
    norm_range_vector <- apply(master_neighborhood_c[nu_panel], 2, function(x) quantile(x, 0.99, names = F)) 
    master_neighborhood_c[nu_panel] <- t(t(master_neighborhood_c[nu_panel]) / as.numeric(norm_range_vector)) 
    
    #calculate heatmap
    master_neighborhood_c %>%
      dplyr::group_by(run_type_name, OtherModel4) %>%
      dplyr::summarize_if(is.numeric, funs(median)) %>%
      ungroup() ->
      heat_neighborhood
    
    #organize naming, ordering schema  
      rownames(heat_neighborhood) <- paste0(heat_neighborhood$run_type_name, heat_neighborhood$OtherModel4)
      #heat[heat>1] <- 1
      heatmaply(heatmaply::normalize(heat_neighborhood[c(2:7,9:14,16:21),panel][-c(1:4,44:45)]), Rowv = FALSE, grid_gap = 1, row_side_colors = c(rep(1,6), rep(2, 6), rep(3,6)))
      
    #heatmaply(normalize(heat[1:6,nu_panel]))
    heatmaply(normalize(heat_neighborhood[c(2:7),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # Ctrl
    heatmaply(normalize(heat_neighborhood[c(9:14),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # MedAD
    heatmaply(normalize(heat_neighborhood[c(16:21),panel][-c(1:4,44:45)]), Rowv = FALSE, Colv = FALSE) # HiAD
    #heatmaply(heat[,clustering_markers])
  
  # Full distributions --> expr data distributions for object types, by region
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

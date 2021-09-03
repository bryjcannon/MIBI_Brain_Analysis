#### Figure SubNeurons-> Identifying cell and disease object compositional differences in data-driven segregated human brain by AD status ####
  
#### A) UMAP of subclustered neurons ####  
  # Meta  
  p1_U <- filter(fsom_neurons) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color=Meta)) +
      geom_density_2d(color="#CCCCCC") +
      geom_point_rast(size = 5) + # original size = 0.5
      scale_color_manual(values = color6) +
      coord_fixed(ratio = 1) +
      xlim(-7,7) +
      ylim(-9,9) +
      labs(x="UMAP1", y="UMAP2", title="Neuronal subclusters", color="") +
      theme_bw() + 
      theme4 +
      guides(colour = guide_legend(override.aes = list(size=4)))
    #facet_wrap(vars(markers))
    p1_U
    ggsave(filename = 'UMAP_neuron_subclusters_Meta_blank_pt5.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
  # MetaDisease  
  p1_U <- filter(fsom_neurons) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color=Meta)) +
      geom_density_2d(color="#CCCCCC") +
      geom_point_rast(size = 5) + # original size = 0.5
      scale_color_manual(values = color7) +
      coord_fixed(ratio = 1) +
      xlim(-7,7) +
      ylim(-9,9) +
      labs(x="UMAP1", y="UMAP2", title="Neuronal subclusters", color="") +
      theme_bw() + 
      theme4 +
      guides(colour = guide_legend(override.aes = list(size=4)))
    #facet_wrap(vars(markers))
    p1_U
    ggsave(filename = 'UMAP_neuron_subclusters_MetaDisease_blank_pt5.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )  
  # Severity  
  p1_U <- filter(fsom_neurons) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color=run_type_id)) +
      geom_density_2d(color="#CCCCCC") +
      geom_point_rast(size = 5) + # original size = 0.5
      scale_color_manual(values = color4) +
      coord_fixed(ratio = 1) +
      xlim(-7,7) +
      ylim(-9,9) +
      labs(x="UMAP1", y="UMAP2", title="Neuronal subclusters", color="") +
      theme_bw() + 
      theme4 +
      guides(colour = guide_legend(override.aes = list(size=4)))
    #facet_wrap(vars(markers))
    p1_U
    ggsave(filename = 'UMAP_neuron_subclusters_Severity_blank_pt5.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
  # de novo regions  
  p1_U <- filter(fsom_neurons) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color=de_novo_regions)) +
      geom_density_2d(color="#CCCCCC") +
      geom_point_rast(size = 5) + # original size = 0.5
      scale_color_manual(values = color3) +
      coord_fixed(ratio = 1) +
      xlim(-7,7) +
      ylim(-9,9) +
      labs(x="UMAP1", y="UMAP2", title="Neuronal subclusters", color="") +
      theme_bw() + 
      theme4 +
      guides(colour = guide_legend(override.aes = list(size=4)))
    #facet_wrap(vars(markers))
    p1_U
    ggsave(filename = 'UMAP_neuron_subclusters_denovoregions_blank_pt5.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
  # markers
    p1_U <- filter(fsom_neurons) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color=MFN2)) +
      geom_density_2d(color="#CCCCCC") +
      geom_point_rast(size = 5) + # original size = 0.5
      scale_color_gradientn(colours = inferno(100)) +
      coord_fixed(ratio = 1) +
      xlim(-7,7) +
      ylim(-9,9) +
      labs(x="UMAP1", y="UMAP2", title="Neuronal subclusters", color="") +
      theme_bw() + 
      theme4 +
      guides(colour = guide_legend(override.aes = list(size=4))) +
      facet_wrap(vars(de_novo_regions))
    p1_U
    
#### B) Heatmap of markers used to subcluster neurons ####
    single_mean_norm_vector <- apply(fsom_neurons[,panel], 2, mean)
    fsom_neurons_c2 <- fsom_neurons
    fsom_neurons_c2[,panel] <- fsom_neurons_c2[,panel] - single_mean_norm_vector
    
    fsom_neurons_c2 %>%
      dplyr::group_by(Meta) %>%
      dplyr::summarize_if(is.numeric, funs(mean)) %>%
      ungroup() ->
      heat_f
    
    heat_f <- as.data.frame(heat_f)
    heat_f %>% mutate(`Proteopathy-associated`= ifelse(Meta %in% c(1,2,6,8,10), "Proteopathy clusters", "Non-proteopathy clusters")) -> heat_f
    heat_f %>% arrange(`Proteopathy-associated`) -> heat_f
    rownames(heat_f) <- heat_f$Meta
    anno_data = heat_f[,c('Meta', 'Proteopathy-associated')]
    anno_data = as.data.frame(anno_data)
    names(anno_data) = c('Meta', 'Proteopathy-associated')
    mycolors <- list(
      Meta = c(`1`='dodgerblue2',`2`='gold1',`3`='skyblue2',`4`="#FB9A99",`5`="palegreen2",`6`="#CAB2D6",`7`="#FDBF6F",`8`="#002B03",`9`="maroon",`10`="orchid1",`11`="deeppink1"),
      `Proteopathy-associated` = c("Proteopathy clusters" = "red", "Non-proteopathy clusters" = "blue"))
    
    h1u <- pheatmap(heat_f[,clustering_markers], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Neurons", 
                    annotation_row = anno_data, annotation_colors = mycolors, show_rownames = T, annotation_names_row = F,
                    col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 25, width = 15, height = 15,
                    filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4/hm_meta_neurons_meannorm_mean_inferno_rowNOclust.pdf')
    h1u
#### C) Barplot: Subclustered neuron distros across severity, severity ####
    #collect neuron counts and normalize by area
    fsom_neurons %>%
      dplyr::group_by(run_type_id, anat_id, Meta) %>%
      dplyr::summarise(n = n()) %>%
      rowwise() %>%
      mutate(meta_area_norm = n/run_areas[run_type_id]) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Disease clusters", "Non-disease clusters")) -> fsom_neurons_meta_norm
    # barplot: Severity
    C1 <- ggplot(fsom_neurons_meta_norm) +
      aes(x = run_type_id, fill = Meta, weight = meta_area_norm) +
      geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
      scale_fill_manual(labels = c(1,2,3,4,5,6,7,8,9,10,11), values = color6) +
      labs(x = "Severity", y = bquote(bold("Total neurons "~(mm^2))), title = "Neuronal populations across disease status", fill = "Cluster") +
      theme_minimal() +
      theme1a +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#gets rid of lines
      facet_grid(rows = vars(MetaDisease)) +
      theme(panel.spacing = unit(2, "lines"))
    
    C1
    ggsave(filename = 'bar_neuron_cluster_composition_severity.pdf', plot = C1, width = 12, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    ) #2=20,h=20 in v1
    
    #### creating anat & severity area normalized df ####
    fsom_neurons %>%
      dplyr::group_by(run_type_id, anat_id, Meta) %>%
      dplyr::summarise(n = n()) %>%
      rowwise() %>%
      mutate(meta_area_norm = n/run_areas[run_type_id]) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Disease clusters", "Non-disease clusters")) -> fsom_neurons_meta_norm
    
    nu_meta_norm <- c()
    for (i in c("Ctrl","MedAD","HiAD")) {
      for (j in c("DG","CA4","CA3","CA2","CA1")) {
        current_df <- filter(fsom_neurons_meta_norm, run_type_id == i & anat_id == j)
        current_df$meta_area_norm <- current_df$n / filter(pixel_anat_area, Run == i & Mask == j)$area
        nu_meta_norm <- rbind(nu_meta_norm, current_df)
      }
    }
    rm(current_df)
    fsom_neurons_meta_norm %>% na.omit() -> current_df
    nu_meta_norm$MetaDisease <- current_df$MetaDisease
    nu_meta_norm$anat_id <- as.factor(as.character(nu_meta_norm$anat_id))
    nu_meta_norm <- as.data.frame(nu_meta_norm)
      #nu_meta_norm <- semi_join(nu_meta_norm, fsom_neurons_meta_norm, by = c("run_type_id","anat_id","Meta"))
    # barplot: Severity & Anat
    C1a <- filter(nu_meta_norm, anat_id %in% c("CA1","CA2")) %>%
      ggplot() +
      aes(x = anat_id, fill = Meta, weight = meta_area_norm) +
      geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
      scale_fill_manual(labels = c(1,2,3,4,5,6,7,8,9,10,11), values = color6) +
      labs(x = "Severity", y = bquote(bold("Total neurons "~(mm^2))), title = "Neuronal populations across disease status", fill = "Cluster") +
      theme_minimal() +
      theme1a +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#gets rid of lines
      facet_grid(rows = vars(MetaDisease), cols = vars(run_type_id)) +
      theme(panel.spacing = unit(2, "lines"))
    
    C1a
    ggsave(filename = 'bar_neuron_cluster_composition_severity_anat_CA1CA2.pdf', plot = C1a, width = 12, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    ) #2=20,h=20 in v1
    
#### D) Barplot: Subclustered neuron distros across severity, de novo regions ####
    #collect neuron counts and normalize by num of neurons in de novo region
    fsom_neurons_total_denovo <- t(fsom_neurons %>% dplyr::group_by(de_novo_regions) %>% summarise(n = n()) %>% select(n))
    fsom_neurons %>%
      dplyr::group_by(de_novo_regions, Meta) %>%
      dplyr::summarise(n = n()) %>%
      rowwise() %>%
      mutate(meta_n_norm = n/fsom_neurons_total_denovo[de_novo_regions]) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Proteopathy clusters", "Non-proteopathy clusters")) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm %>% mutate(MetaMFN = ifelse(Meta %in% c(1,2,3,4,10), "MFN2 +", "MFN2 -")) -> fsom_neurons_meta_norm
      
    
    # barplot: De novo regions - proportional
    fsom_neurons_meta_norm$dnr2 = factor(fsom_neurons_meta_norm$de_novo_regions, levels = c(1,2,5,3,4), labels = c('ND','GD','MD','APD','TTD'))
      D1 <- ggplot(filter(fsom_neurons_meta_norm, de_novo_regions %in% c(1,2,3,4,5))) +
        aes(x = dnr2, fill = Meta, weight = meta_n_norm) +
        geom_bar(position = "stack", width = 0.8) +
        #geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
        scale_fill_manual(labels = c(1,2,3,4,5,6,7,8,9,10,11), values = color6) +
        labs(x = "Regions", y = "Cluster proportions (relative to total region)", title = "MFN2+ Populations persist in proteopathy-laden regions", fill = "Cluster") +
        #scale_x_discrete(labels=c("Disease", "No Disease")) +
        theme_minimal() +
        theme1a +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#gets rid of lines
        facet_grid(rows = vars(MetaDisease)) +
        theme(panel.spacing = unit(2, "lines"))
      
      D1
      ggsave(filename = 'bar_neuron_cluster_composition_denovoregions_STACKED_all_proportional.pdf', plot = D1, width = 15, height = 15, dpi = 300, device = 'pdf',
             path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4'
      ) #w=20,h=20 in v1, w=12,h=10 in v2 
    
    # barplot: De novo regions - absolute
      D1 <- ggplot(filter(fsom_neurons_meta_norm, de_novo_regions %in% c(1,2,3,4,5))) +
        aes(x = dnr2, fill = MetaMFN, weight = meta_n_norm) +
        geom_bar(position = "stack", width = 0.8) +
        #geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
        #scale_fill_manual(labels = c(1,2,3,4,5,6,7,8,9,10,11), values = color6) +
        scale_fill_manual(labels = c("MFN2 lo","MFN2 hi"), values = c("purple","orange")) +
        labs(x = "Regions", y = "Cluster proportions (relative to total region)", title = "MFN2 hi Populations persist in proteopathy-laden regions", fill = "Cluster groups") +
        #scale_x_discrete(labels=c("Disease", "No Disease")) +
        theme_minimal() +
        theme1a +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#gets rid of lines
        #facet_grid(rows = vars(MetaDisease)) +
        theme(panel.spacing = unit(2, "lines"))
      
      D1
      ggsave(filename = 'bar_neuron_cluster_composition_denovoregions_STACKED_all_proportional_MetaMFN.pdf', plot = D1, width = 15, height = 15, dpi = 300, device = 'pdf',
             path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4'
      ) #w=20,h=20 in v1, w=12,h=10 in v2 
      
  #severity & de novo -> selected cluster breakdown
    #collect neuron counts and normalize by num of neurons in de novo region
    fsom_neurons_total_denovo <- filter(fsom_neurons, anat_id %in% c("DG","CA4","CA3","CA2","CA1")) %>% dplyr::group_by(de_novo_regions, run_type_id) %>% summarise(n_total = n())
    filter(fsom_neurons, anat_id %in% c("DG","CA4","CA3","CA2","CA1")) %>%
      dplyr::group_by(de_novo_regions, run_type_id, Meta) %>%
      dplyr::summarise(n = n()) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm <- full_join(fsom_neurons_meta_norm, fsom_neurons_total_denovo, by = c("de_novo_regions","run_type_id"))
    fsom_neurons_meta_norm %>% 
      rowwise() %>%
      mutate(meta_n_norm = n/n_total) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Disease clusters", "Non-disease clusters")) -> fsom_neurons_meta_norm
    # barplot: De novo regions by severity
    fsom_neurons_meta_norm$dnr2 = factor(fsom_neurons_meta_norm$de_novo_regions, levels = c(1,2,5,3,4), labels = c('ND','GD','MD','APD','TTD'))
    D1 <- ggplot(filter(fsom_neurons_meta_norm, de_novo_regions %in% c(1,4,5) & Meta %in% c(1,2,3,4))) +
      aes(x = dnr2, fill = Meta, weight = meta_n_norm) +
      geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
      scale_fill_manual(labels = c(1,2,3,4,5,6,7,8,9,10,11), values = color6) +
      labs(x = "Cluster types", y = "Cluster proportions (relative to total region)", title = "MFN2+ neuronal populations survive as disease pathology worsens", fill = "Cluster") +
      #scale_x_discrete(labels=c("Disease", "No Disease")) +
      theme_minimal() +
      theme1a +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#gets rid of lines
      facet_grid(rows = vars(MetaDisease)) +
      theme(panel.spacing = unit(2, "lines"))
    
    D1
    ggsave(filename = 'bar_neuron_cluster_composition_denovoregions_145ONLY.pdf', plot = D1, width = 12, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    ) #2=20,h=20 in v1
    
  #just de novo -> selected cluster breakdown
    #collect neuron counts and normalize by num of neurons in de novo region
    fsom_neurons_total_denovo <- filter(fsom_neurons) %>% dplyr::group_by(de_novo_regions) %>% summarise(n_total = n())
    filter(fsom_neurons) %>%
      dplyr::group_by(de_novo_regions, Meta) %>%
      dplyr::summarise(n = n()) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm <- full_join(fsom_neurons_meta_norm, fsom_neurons_total_denovo, by = c("de_novo_regions"))
    fsom_neurons_meta_norm %>% 
      rowwise() %>%
      mutate(meta_n_norm = n/n_total) -> fsom_neurons_meta_norm
    fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Disease clusters", "Non-disease clusters")) -> fsom_neurons_meta_norm
    # barplot: De novo regions by USE THIS LATE JUNE
    fsom_neurons_meta_norm$dnr2 = factor(fsom_neurons_meta_norm$de_novo_regions, levels = c(1,2,5,3,4), labels = c('ND','GD','MD','APD','TTD'))
    D1 <- ggplot(filter(fsom_neurons_meta_norm, dnr2 %in% c('ND','GD','MD','APD','TTD') & Meta %in% c(1,2,3,4,10))) +
      aes(x = dnr2, fill = Meta, weight = meta_n_norm) +
      geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +
      scale_fill_manual(labels = c(1,2,3,4,10), values = color6) +
      labs(x = "De novo regions", y = "Cluster proportions (relative to total region)", title = "MFN2+ neuronal populations survive as disease pathology worsens", fill = "Cluster") +
      #scale_x_discrete(labels=c("Disease", "No Disease")) +
      theme_minimal() +
      theme1a +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +#gets rid of lines
      facet_grid(rows = vars(MetaDisease)) +
      theme(panel.spacing = unit(2, "lines"))
    
    D1
    ggsave(filename = 'bar_neuron_cluster_composition_denovoregions_145ONLY.pdf', plot = D1, width = 12, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    ) #2=20,h=20 in v1
#### E) Lineplot: MFN2 changes across subclusters neurons ####
    install.packages("ggnewscale")
    library("ggnewscale")
    # De novo regions
    fsom_neurons$dnr2 = factor(fsom_neurons$de_novo_regions, levels = c(1,2,5,3,4), labels = c('ND','GD','MD','APD','TTD'))
    fsom_neurons$MetaDisease_subclass <- NULL
    fsom_neurons$MetaDisease_subclass[fsom_neurons$Meta %in% c(1,2)] = "Both"
    fsom_neurons$MetaDisease_subclass[fsom_neurons$Meta %in% c(10)] = "Amyloid-disease alone"
    fsom_neurons$MetaDisease_subclass[fsom_neurons$Meta %in% c(8)] = "Tau disease alone"
    fsom_neurons$MetaDisease_subclass[fsom_neurons$Meta %in% c(3,4)] = "Disease-free"
    fsom_neurons$MetaDisease_subclass[fsom_neurons$Meta %in% c(5,6,7,9,11)] = "MFN2-lo"
    fsom_neurons$MetaDisease_subclass <- factor(fsom_neurons$MetaDisease_subclass, levels = c("Amyloid-disease alone","Tau disease alone",
                                                                                              "Both","Disease-free","MFN2-lo"), labels = c("Amyloid-disease alone","Tau disease alone",
                                                                                                                                           "Both","Disease-free","MFN2-lo"))
    # De novo regions - lineplot
    E1a <- filter(fsom_neurons, Meta %in% c(1,2,3,4,8,10)) %>%
      ggplot() +
      aes(x = dnr2, y = MFN2, group = Meta, color = Meta) +
      geom_point(stat='summary', fun=mean, size = 5) +
      stat_summary(fun=mean, geom="line", size = 0.25) + # normally size=2
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      scale_color_manual(values = color6) +
      scale_x_discrete(labels=c('ND','GD','MD','APD','TTD')) +
      labs(x = "De novo region", y = "Normalized neuronal MFN2 expression", title = "MFN2 decreases as tau disease spreads in surviving neuronal populations", fill = "Cluster") +
      theme_minimal() +
      theme1a +
      facet_grid(cols = vars(MetaDisease_subclass)) +
      theme(panel.spacing = unit(2, "lines"))
    E1a
    ggsave(filename = 'pointline_neuron_cluster_MFN_denovoregions.pdf', plot = E1a, width = 20, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4'
    )
      # De novo regions - barplot version of lineplot above
      E1a <- filter(fsom_neurons, Meta %in% c(1,2,3,4,8,10)) %>%
        ggplot() +
        aes(x = dnr2, y = MFN2, group = Meta, fill = Meta) +
        #geom_bar(stat='summary', fun.y=mean, size = 5) +
        stat_summary(fun = mean, geom="bar", position=position_dodge2(width=0.95, preserve = "single")) + # normal`ly size=2
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position=position_dodge2(width=0.95, preserve = "single")) +
        scale_fill_manual(values = color6) +
        scale_x_discrete(labels=c('ND','GD','MD','APD','TTD')) +
        labs(x = "De novo region", y = "Normalized neuronal MFN2 expression", title = "MFN2 decreases as tau disease spreads in surviving neuronal populations", fill = "Cluster") +
        theme_minimal() +
        theme1a +
        facet_grid(cols = vars(MetaDisease_subclass)) +
        theme(panel.spacing = unit(2, "lines"))
      E1a
      ggsave(filename = 'bar_neuron_cluster_MFN_denovoregions_separated.pdf', plot = E1a, width = 20, height = 10, dpi = 300, device = 'pdf',
             path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4'
      )
    
    
    # Severity
    E2 <- filter(fsom_neurons, de_novo_regions %in% c(1,2,3,4,5) )%>%#& Meta %in% c(1,2,3,4,10)) %>%
      mutate(MetaMFN = ifelse(Meta %in% c(1,2,3,4,10), "MFN2 hi", "MFN2 lo")) %>%
      ggplot() +
      aes(x = run_type_id, y = MFN2, group = Meta, color = Meta) +
      geom_point(stat='summary', fun=mean, size = 5) +
      stat_summary(fun=mean, geom="line", size = 2) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      scale_color_manual(values = color6) +
      labs(x = "Severity", y = "Mean normalized expression", title = "MFN2 higher in resistant neuronal populations", fill = "Cluster") +
      theme_minimal() +
      theme1a +
      facet_grid(cols = vars(MetaDisease)) +
      theme(panel.spacing = unit(2, "lines"))
    E2
    ggsave(filename = 'pointline_neuron_cluster_MFN_severity.pdf', plot = E2, width = 20, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    # Severity vs Anat - with Meta
    E2a <- filter(fsom_neurons, de_novo_regions %in% c(1,2,3,4,5) & anat_id %in% c("CA1","CA2","CA3"))%>%#& Meta %in% c(1,2,3,4,10)) %>%
      mutate(MetaMFN = ifelse(Meta %in% c(1,2,3,4,10), "MFN2 hi", "MFN2 lo")) %>%
      ggplot() +
      aes(x = anat_id, y = MFN2, group = Meta, color = Meta) +
      geom_point(stat='summary', fun=mean, size = 5) +
      stat_summary(fun=mean, geom="line", size = 2) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      scale_color_manual(values = color6) +
      labs(x = "Severity", y = "Mean normalized expression", title = "MFN2 higher in resistant neuronal populations", fill = "Cluster") +
      theme_minimal() +
      theme1a +
      facet_grid(rows = vars(MetaDisease), cols = vars(run_type_id)) +
      theme(panel.spacing = unit(2, "lines"))
    E2a
    ggsave(filename = 'pointline_neuron_cluster_MFN_severity_anat.pdf', plot = E2a, width = 26, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    # Anat vs Severity - no Meta
    E2b <- filter(fsom_neurons, de_novo_regions %in% c(1,2,3,4,5) )%>%#& Meta %in% c(1,2,3,4,10)) %>%
      mutate(MetaMFN = ifelse(Meta %in% c(1,2,3,4,10), "MFN2 hi", "MFN2 lo")) %>%
      ggplot() +
      aes(x = anat_id, y = MFN2) +
      geom_point(stat='summary', fun=sum, size = 5) +
      stat_summary(fun=sum, geom="line", size = 2) +
      #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      scale_color_manual(values = color6) +
      labs(x = "Severity", y = "Mean normalized expression", title = "MFN2 higher in resistant neuronal populations", fill = "Cluster") +
      theme_minimal() +
      theme1a +
      facet_grid(cols = vars(run_type_id)) +
      theme(panel.spacing = unit(2, "lines"))
    E2b
    ggsave(filename = 'pointline_neuron_cluster_MFN_severity_anat_nometa_summed.pdf', plot = E2b, width = 26, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
#### F) Correlation ####
    #multiple levels - severity
    for (i in c(1:11)) {
      cor_sample <- filter(fsom_neurons, Meta == i)
      for (j in c("Ctrl","MedAD","HiAD")) {
        cor_sample2 <- filter(cor_sample, run_type_id == j)
        neu_cor <- cor(cor_sample2[c("MFN2","PHF1Tau","Amyloidbeta142","Amyloidbeta140")])
        pdf(file = paste0("/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2/corrplots/corrplot_",i,"_",j,".pdf"))
        corrplot(neu_cor, type = "upper", 
                 tl.col = "black", tl.srt = 45)
        dev.off()
      }
    }
    #multiple levels - de novo regions
    for (i in c(1:11)) {
      cor_sample <- filter(fsom_neurons, Meta == i)
      for (j in c(1,2,3,4,5)) {
        cor_sample2 <- filter(cor_sample, de_novo_regions == j)
        if (dim(cor_sample2)[1] < 2) next
        neu_cor <- cor(cor_sample2[c("MFN2","PolyubiK48","Calbindin","PHF1Tau","HistoneH3Lyo")])
        pdf(file = paste0("/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2/corrplots/nu_corrplot_",i,"_",j,".pdf"))
        corrplot(neu_cor, type = "upper", method = "number",
                 tl.col = "black", tl.srt = 45)
        dev.off()
      }
    }
    #multiple levels - anat_id
    for (i in c("Ctrl","MedAD","HiAD")) {
      cor_sample <- filter(fsom_neurons, run_type_id == i)
      for (j in c("DG","CA4","CA3","CA2","CA1")) {
        cor_sample2 <- filter(cor_sample, anat_id == j)
        if (dim(cor_sample2)[1] < 2) next
        neu_cor <- cor(cor_sample2[c("MFN2","PHF1Tau","Amyloidbeta142","Amyloidbeta140")])
        pdf(file = paste0("/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2/corrplots/corrplot_",i,"_",j,".pdf"))
        corrplot(neu_cor, type = "upper",
                 tl.col = "black", tl.srt = 45)
        dev.off()
      }
    }
    #single level
    for (i in c(1:11)) {
      cor_sample <- filter(fsom_neurons, Meta == i)
      neu_cor <- cor(cor_sample[c("MFN2","PHF1Tau","Amyloidbeta142","Amyloidbeta140")])
      pdf(file = paste0("/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2/corrplots/corrplot_",i,".pdf"))
      corrplot(neu_cor, type = "upper", method = 'number',
               tl.col = "black", tl.srt = 45)
      dev.off()
    }
#### G1) Biaxial plot (MFN2 & PHF1Tau) ####
    anat_dive_data <- fsom_neurons %>% filter(anat_id %in% c("CA1","CA2","DG") & Meta %in% c(1,2,3,4,8,10))
    anat_dive_subset <- filter(anat_dive_data, MFN2 > 0 & PHF1Tau == 0)
    anat_dive_subset_inverse <- filter(anat_dive_data, MFN2 == 0 & PHF1Tau > 0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[c(1:3,10:15),] #grab just CA1, CA2, DG rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Percent")
    round(anat_dive_percents["Percent"], digits=2) %>% mutate(textPercent = paste0(Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    
    G1 <- filter(fsom_neurons, run_type_id %in% c("Ctrl","MedAD","HiAD") & anat_id %in% c("CA1") & Meta %in% c(1,2,3,4,8,10)) %>%
      ggplot() +
      aes(x = PHF1Tau, y = MFN2, color = Meta) +
      geom_jitter(size = 1L, width = 0.02) +
      scale_color_manual(values = color6) +
      #stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      #geom_rect(xmin=-0.04,xmax=0.04,ymin=0.02,ymax=1.1,color="black",alpha=0) +
      #geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      coord_fixed(ratio = 1) +
      theme_bw() +
      theme1a +
      labs(title="Tau-Disease Resistant Neurons predominant in DG are MFN2+", color="Meta") +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      facet_grid(rows=vars(run_type_id), cols=vars(anat_id)) +
      theme(panel.spacing = unit(2, "lines"))
      #geom_text(data=anat_dive_percents, aes(x=0.2, y=0.9, label=textPercent, size=50), show.legend = FALSE, inherit.aes = FALSE) #add percent dbl positive to each subplot
      
    G1
    ggsave(filename = 'biaxial_MFN2_PHF1Tau_SubNeurons_DGCA2CA1_MedAD.pdf', plot = G1, width = 15, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    
    # CA1 & CA2 --> MetaDisease (PHF1Tau and Histone H3Lyo)
    G1a <- filter(fsom_neurons, run_type_id %in% c("Ctrl","MedAD","HiAD") & anat_id %in% c("CA1","CA2") & Meta %in% c(1,2,3,4)) %>%
      ggplot() +
      aes(x = PHF1Tau, y = MFN2, color = Meta) +
      geom_jitter(size = 1L, width = 0.02) +
      scale_color_manual(values = color6) +
      #stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      #geom_rect(xmin=-0.04,xmax=0.04,ymin=0.02,ymax=1.1,color="black",alpha=0) +
      #geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      coord_fixed(ratio = 1) +
      xlim(-0.05,1.25) +
      ylim(-0.05,0.7) +
      theme_bw() +
      theme1a +
      labs(title="Tau-Disease Resistant Neurons predominant in DG are MFN2+", color="Meta") +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      facet_grid(cols=vars(anat_id)) +
      theme(panel.spacing = unit(2, "lines"))
    
    G1a
    ggsave(filename = 'biaxial_MFN2_PHF1Tau_SubNeurons_DGCA2CA1_MedAD.pdf', plot = G1a, width = 15, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    
#### G2) Biaxial plot (VGLUT1 & PHF1Tau) ####
    anat_dive_data <- fsom_neurons %>% filter(anat_id %in% c("CA1","CA2","DG") & Meta %in% c(1,2,3,4,8,10))
    anat_dive_subset <- filter(anat_dive_data, VGLUT1 > 0 & PHF1Tau == 0)
    anat_dive_subset_inverse <- filter(anat_dive_data, VGLUT1 == 0 & PHF1Tau > 0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[c(1:3,10:15),] #grab just CA1, CA2, DG rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Percent")
    round(anat_dive_percents["Percent"], digits=2) %>% mutate(textPercent = paste0(Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    
    G2 <- filter(fsom_neurons, anat_id %in% c("DG","CA2","CA1") & Meta %in% c(1,2,3,4,8,10)) %>%
      ggplot() +
      aes(x = PHF1Tau, y = VGLUT1, color = Meta) +
      geom_jitter(size = 1L, width = 0.02) +
      scale_color_manual(values = color6) +
      #stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      #geom_rect(xmin=-0.04,xmax=0.04,ymin=0.02,ymax=1.1,color="black",alpha=0) +
      #geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      coord_fixed(ratio = 1) +
      theme_bw() +
      theme1a +
      labs(title="Tau-Disease Resistant Neurons have lower VGLUT1", color="Meta") +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      facet_grid(rows=vars(run_type_id), cols=vars(anat_id)) +
      theme(panel.spacing = unit(2, "lines"))
      #geom_text(data=anat_dive_percents, aes(x=0.2, y=0.9, label=textPercent, size=50), show.legend = FALSE, inherit.aes = FALSE) #add percent dbl positive to each subplot
    
    G2
    ggsave(filename = 'biaxial_VGLUT1_PHF1Tau_SubNeurons_DGCA2CA1.pdf', plot = G2, width = 15, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    
#### G3) Biaxial plot (HistoneH3Lyo & PHF1Tau) ####
    anat_dive_data <- fsom_neurons %>% filter(anat_id %in% c("CA1","CA2","DG") & Meta %in% c(1,2,3,4,8,10))
    anat_dive_subset <- filter(anat_dive_data, HistoneH3Lyo > 0 & PHF1Tau == 0)
    anat_dive_subset_inverse <- filter(anat_dive_data, HistoneH3Lyo == 0 & PHF1Tau > 0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[c(1:3,10:15),] #grab just CA1, CA2, DG rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Percent")
    round(anat_dive_percents["Percent"], digits=2) %>% mutate(textPercent = paste0(Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    
    G3 <- filter(fsom_neurons, anat_id %in% c("DG","CA2","CA1") & Meta %in% c(1,2,3,4,8,10)) %>%
      ggplot() +
      aes(x = PHF1Tau, y = HistoneH3Lyo, color = Meta) +
      geom_jitter(size = 1L, width = 0.02) +
      scale_color_manual(values = color6) +
      #stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      #geom_rect(xmin=-0.04,xmax=0.04,ymin=0.02,ymax=1.1,color="black",alpha=0) +
      #geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      coord_fixed(ratio = 1) +
      theme_bw() +
      theme1a +
      labs(title="Neurons have similar DNA markers across anatomy", color="Meta") +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      facet_grid(rows=vars(run_type_id), cols=vars(anat_id)) +
      theme(panel.spacing = unit(2, "lines"))
      #geom_text(data=anat_dive_percents, aes(x=0.2, y=0.9, label=textPercent, size=50), show.legend = FALSE, inherit.aes = FALSE) #add percent dbl positive to each subplot
    
    G3
    ggsave(filename = 'biaxial_HistoneH3Lyo_PHF1Tau_SubNeurons_DGCA2CA1.pdf', plot = G3, width = 15, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    
#### I) Biaxials #####
    # basic biaxial - PHFTau & PolyubiK48,PolyubiK63
    master_gated_anat2 %>%
      filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2")) %>%
      ggplot() +
      aes(x = PHF1Tau, y = PolyubiK48, color = run_type_id) +
      geom_point(size = 1L) +
      scale_color_manual(values = color4) +
      theme_bw() +
      theme1a +
      facet_grid(vars(anat_id), vars(run_type_id))
    
    #### PHFTau & PolyubiK48 -- with regression on positive points (also include double positive event %'s) ####
    library("ggpubr")
    denovo_dive_data <- fsom_neurons
    denovo_dive_subset <- filter(denovo_dive_data, PolyubiK63>0 & PHF1Tau>0)
    denovo_dive_subset_inverse <- filter(denovo_dive_data, PolyubiK63==0 | PHF1Tau==0)
    denovo_dive_percents <- table(denovo_dive_subset$dnr2) / table(denovo_dive_data$dnr2) * 100
    denovo_dive_percents
    denovo_dive_percents <- as.data.frame(denovo_dive_percents)[1:5,] #grab just CA1, CA2 rows
    names(denovo_dive_percents) <- c("De novo regions", "Dbl_Pos_Percent")
    round(denovo_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> denovo_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- denovo_dive_data %>%
      ggplot() +
      aes(x = PHF1Tau, y = PolyubiK48) +
      geom_point(size = 1L) +
      scale_color_manual(values = color7) +
      stat_smooth(data = filter(denovo_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.42,color="black",alpha=0) +
      geom_point(data=denovo_dive_subset_inverse, size=1L, color="gray") +
      theme_bw() +
      theme1a +
      labs(title="Polyubiquitins & Disease within Neurons", color="Cluster type") +
      facet_grid(cols=vars(dnr2)) +
      stat_cor(data=anat_dive_subset, label.x=0.3, label.y=0.4) + #calculate correlation of each dbl positive subpopulation
      stat_regline_equation(data=denovo_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.3, label.y=0.38) + #caluclate and add information on corr, regression line to each subplot
      geom_text(data=denovo_dive_percents, aes(x=0.9, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      theme(legend.position="none") +
      theme(panel.spacing = unit(2, "lines"))
    
    G5
    ggsave(filename = 'biaxial_Pk48_PHF1Tau_DGtoCA1_MetaDisease_noColor.pdf', plot = G5, width = 15, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
  
    #### PHFTau & PolyubiK63 -- with regression on positive points (also include double positive event %'s) ####
    anat_dive_data <- master_gated_anat2 %>% filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2"))
    anat_dive_subset <- filter(anat_dive_data, PolyubiK63>0 & PHF1Tau>0)
    anat_dive_subset_inverse <- filter(anat_dive_data, PolyubiK63==0 | PHF1Tau==0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[10:15,] #grab just CA1, CA2 rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Dbl_Pos_Percent")
    round(anat_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- anat_dive_data %>%
      ggplot() +
      aes(x = PHF1Tau, y = PolyubiK63, color = run_type_id) +
      geom_point(size = 1L) +
      scale_color_manual(values = color4) +
      stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.365,color="black",alpha=0) +
      geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      theme_bw() +
      theme1a +
      labs(title="Polyubiquitins & Disease within Neurons", color="Severity") +
      facet_grid(rows=vars(anat_id), cols=vars(run_type_id)) +
      stat_cor(data=anat_dive_subset, label.x=0.3, label.y=0.34) + #calculate correlation of each dbl positive subpopulation
      stat_regline_equation(data=anat_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.3, label.y=0.32) + #caluclate and add information on corr, regression line to each subplot
      geom_text(data=anat_dive_percents, aes(x=0.9, y=0.29, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      theme(legend.position="none")
    
    G5
    ggsave(filename = 'biaxial_Pk63_PHF1Tau_DGtoCA1.pdf', plot = G5, width = 15, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
    )
    
    #### Amyloidbeta142 & PolyubiK48 -- with regression on positive points (also include double positive event %'s) ####
    anat_dive_data <- master_gated_anat2 %>% filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2"))
    anat_dive_subset <- filter(anat_dive_data, PolyubiK48>0 & Amyloidbeta142>0)
    anat_dive_subset_inverse <- filter(anat_dive_data, PolyubiK48==0 | Amyloidbeta142==0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[10:15,] #grab just CA1, CA2 rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Dbl_Pos_Percent")
    round(anat_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- anat_dive_data %>%
      ggplot() +
      aes(x = Amyloidbeta142, y = PolyubiK48, color = run_type_id) +
      geom_point(size = 1L) +
      scale_color_manual(values = color4) +
      stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      geom_rect(xmin=0.03,xmax=0.45,ymin=0.01,ymax=0.42,color="black",alpha=0) +
      geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      theme_bw() +
      theme1a +
      labs(title="Polyubiquitins & Disease within Neurons", color="Severity") +
      facet_grid(rows=vars(anat_id), cols=vars(run_type_id)) +
      stat_cor(data=anat_dive_subset, label.x=0.2, label.y=0.4) + #calculate correlation of each dbl positive subpopulation
      stat_regline_equation(data=anat_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.2, label.y=0.38) + #caluclate and add information on corr, regression line to each subplot
      geom_text(data=anat_dive_percents, aes(x=0.4, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      theme(legend.position="none")
    
    G5
    ggsave(filename = 'biaxial_Pk48_Abeta42_DGtoCA1.pdf', plot = G5, width = 15, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
    )
    
    #### Amyloidbeta140 & PolyubiK48 -- with regression on positive points (also include double positive event %'s) ####
    anat_dive_data <- master_gated_anat2 %>% filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2"))
    anat_dive_subset <- filter(anat_dive_data, PolyubiK48>0 & Amyloidbeta140>0)
    anat_dive_subset_inverse <- filter(anat_dive_data, PolyubiK48==0 | Amyloidbeta140==0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[10:15,] #grab just CA1, CA2 rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Dbl_Pos_Percent")
    round(anat_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- anat_dive_data %>%
      ggplot() +
      aes(x = Amyloidbeta140, y = PolyubiK48, color = run_type_id) +
      geom_point(size = 1L) +
      scale_color_manual(values = color4) +
      stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      geom_rect(xmin=0.03,xmax=0.7,ymin=0.01,ymax=0.42,color="black",alpha=0) +
      geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      theme_bw() +
      theme1a +
      labs(title="Polyubiquitins & Disease within Neurons", color="Severity") +
      facet_grid(rows=vars(anat_id), cols=vars(run_type_id)) +
      stat_cor(data=anat_dive_subset, label.x=0.2, label.y=0.4) + #calculate correlation of each dbl positive subpopulation
      stat_regline_equation(data=anat_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.2, label.y=0.38) + #caluclate and add information on corr, regression line to each subplot
      geom_text(data=anat_dive_percents, aes(x=0.55, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      theme(legend.position="none")
    
    G5
    ggsave(filename = 'biaxial_Pk48_Abeta40_DGtoCA1.pdf', plot = G5, width = 15, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
    )
    
    
    #### Amyloidbeta142 & PolyubiK63 -- with regression on positive points (also include double positive event %'s) ####
    anat_dive_data <- master_gated_anat2 %>% filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2"))
    anat_dive_subset <- filter(anat_dive_data, PolyubiK63>0 & Amyloidbeta142>0)
    anat_dive_subset_inverse <- filter(anat_dive_data, PolyubiK63==0 | Amyloidbeta142==0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[10:15,] #grab just CA1, CA2 rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Dbl_Pos_Percent")
    round(anat_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- anat_dive_data %>%
      ggplot() +
      aes(x = Amyloidbeta142, y = PolyubiK63, color = run_type_id) +
      geom_point(size = 1L) +
      scale_color_manual(values = color4) +
      stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      geom_rect(xmin=0.03,xmax=0.45,ymin=0.01,ymax=0.388,color="black",alpha=0) +
      geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      theme_bw() +
      theme1a +
      labs(title="Polyubiquitins & Disease within Neurons", color="Severity") +
      facet_grid(rows=vars(anat_id), cols=vars(run_type_id)) +
      stat_cor(data=anat_dive_subset, label.x=0.16, label.y=0.375) + #calculate correlation of each dbl positive subpopulation
      stat_regline_equation(data=anat_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.16, label.y=0.355) + #caluclate and add information on corr, regression line to each subplot
      geom_text(data=anat_dive_percents, aes(x=0.4, y=0.325, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      theme(legend.position="none")
    
    G5
    ggsave(filename = 'biaxial_Pk63_Abeta42_DGtoCA1.pdf', plot = G5, width = 15, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
    )
    
    #### Amyloidbeta140 & PolyubiK63 -- with regression on positive points (also include double positive event %'s) ####
    anat_dive_data <- master_gated_anat2 %>% filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2"))
    anat_dive_subset <- filter(anat_dive_data, PolyubiK63>0 & Amyloidbeta140>0)
    anat_dive_subset_inverse <- filter(anat_dive_data, PolyubiK63==0 | Amyloidbeta140==0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[10:15,] #grab just CA1, CA2 rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Dbl_Pos_Percent")
    round(anat_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- anat_dive_data %>%
      ggplot() +
      aes(x = Amyloidbeta140, y = PolyubiK63, color = run_type_id) +
      geom_point(size = 1L) +
      scale_color_manual(values = color4) + 
      stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      geom_rect(xmin=0.03,xmax=0.7,ymin=0.01,ymax=0.365,color="black",alpha=0) +
      geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
      theme_bw() +
      theme1a +
      labs(title="Polyubiquitins & Disease within Neurons", color="Severity") +
      facet_grid(rows=vars(anat_id), cols=vars(run_type_id)) +
      stat_cor(data=anat_dive_subset, label.x=0.2, label.y=0.3) + #calculate correlation of each dbl positive subpopulation
      stat_regline_equation(data=anat_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.2, label.y=0.28) + #caluclate and add information on corr, regression line to each subplot
      geom_text(data=anat_dive_percents, aes(x=0.55, y=0.25, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      theme(legend.position="none")
    
    G5
    ggsave(filename = 'biaxial_Pk63_Abeta40_DGtoCA1.pdf', plot = G5, width = 15, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
    )
    
    
    #### MFN2 & PHF1Tau /  Abeta42 (Final w/ Gates)####
    denovo_dive_data <- fsom_neurons
      #denovo_dive_data[panel] <- log(denovo_dive_data[panel] + 1) # for any log transformations
    denovo_dive_data_subset <- filter(denovo_dive_data, MFN2>0 & PHF1Tau>0)
    denovo_dive_data_subset_inverse <- filter(denovo_dive_data, MFN2==0 | PHF1Tau==0)
    denovo_dive_percents <- table(denovo_dive_subset$dnr2) / table(denovo_dive_data$dnr2) * 100
    denovo_dive_percents
    denovo_dive_percents <- as.data.frame(denovo_dive_percents)[1:5,] #grab just CA1, CA2 rows
    names(denovo_dive_percents) <- c("De novo regions", "Dbl_Pos_Percent")
    round(denovo_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> denovo_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- filter(fsom_neurons, Meta %in% c(1,2,3,4,10)) %>%
      ggplot() +
      aes(x = PHF1Tau, y = MFN2, color = Meta) +
      geom_point(size = 2L) +
      scale_color_manual(values = color6) +
      #stat_smooth(data = filter(denovo_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
      #geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.42,color="black",alpha=0) +
      #geom_point(data=denovo_dive_subset_inverse, size=1L, color="gray") +
      coord_fixed(ratio = 1) +
      theme_bw() +
      theme1a +
      labs(title="MFN2 & PHF1Tau within Neurons", color="Cluster type") +
      facet_grid(rows = vars(Meta), cols=vars(dnr2)) +
      #stat_cor(data=anat_dive_subset, label.x=0.3, label.y=0.4) + #calculate correlation of each dbl positive subpopulation
      #stat_regline_equation(data=denovo_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.3, label.y=0.38) + #caluclate and add information on corr, regression line to each subplot
      #geom_text(data=denovo_dive_percents, aes(x=0.9, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      #theme(legend.position="none") +
      theme(panel.spacing = unit(2, "lines"))
    
    G5
    ggsave(filename = 'biaxial_MFN2_PHF1Tau_all_denovoregions_Color.pdf', plot = G5, width = 25, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4'
    )#width = 15, height = 12 v1, 
    
    #plot2 --> USE
    G5 <- denovo_dive_data %>%
      ggplot() +
      aes(x = Amyloidbeta142, y = MFN2, color = de_novo_regions) +
      geom_jitter(size = 0.5L, width = 0.01, height = 0.01) +
      scale_color_manual(labels = c("N","G","AP","TT","M"), values = color3) +
      #stat_smooth(data = filter(denovo_dive_data_subset), method = "glm", formula = y ~ exp(-x^3), color = "black") +
      #geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.42,color="black",alpha=0) +
      #geom_point(data=denovo_dive_data_subset_inverse, size=1L, color="gray") +
      xlim(-0.02,1) +
      ylim(-0.02,1) +
      coord_fixed(ratio=1) +
      theme_bw() +
      theme1a +
      labs(title="MFN2 & Abeta42 within Neurons", color="Cluster type") +
      facet_grid(cols=vars(dnr2)) +
      
      #stat_cor(data=denovo_dive_data, label.x=0.05, label.y=1.1) + #calculate correlation of each dbl positive subpopulation
      #stat_regline_equation(data=denovo_dive_data, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y ~ log(x+1), label.x=0.05, label.y=1) + #caluclate and add information on corr, regression line to each subplot
      #geom_text(data=denovo_dive_percents, aes(x=0.9, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      #theme(legend.position="none") +
      theme(panel.spacing = unit(2, "lines"))
    
    G5
    ggsave(filename = 'biaxial_MFN2_Abeta42_all_denovoregions_Color.pdf', plot = G5, width = 25, height = 12, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4/MFN2_biaxials'
    )
    
  #plot3 --> breaking down MFN2+ clusters into relationship between MFN2 and progressive disease markers
    denovo_dive_data <- fsom_neurons
    #denovo_dive_data[panel] <- log(denovo_dive_data[panel] + 1) # for any log transformations
    #percents
    denovo_dive_data <- filter(denovo_dive_data, MFN2>0 & PHF1Tau>0)
      denovo_dive_data_topright <- filter(denovo_dive_data, MFN2>=0.25 & PHF1Tau>=0.25)
      denovo_dive_data_bottomright <- filter(denovo_dive_data, MFN2<0.25 & PHF1Tau>=0.25)
      denovo_dive_data_topleft <- filter(denovo_dive_data, MFN2>=0.25 & PHF1Tau<0.25)
      denovo_dive_data_bottomleft <- filter(denovo_dive_data, MFN2<0.25 & PHF1Tau<0.25)
      
      denovo_dive_data_percents <- c(dim(denovo_dive_data_topright)[1], dim(denovo_dive_data_bottomright)[1], dim(denovo_dive_data_topleft)[1], dim(denovo_dive_data_bottomleft)[1])
      denovo_dive_data_percents <- 100*(denovo_dive_data_percents / dim(denovo_dive_data)[1])
      
    denovo_dive_data <- filter(denovo_dive_data, MFN2>0 & Amyloidbeta142>0)
      denovo_dive_data_topright <- filter(denovo_dive_data, MFN2>=0.25 & Amyloidbeta142>=0.25)
      denovo_dive_data_bottomright <- filter(denovo_dive_data, MFN2<0.25 & Amyloidbeta142>=0.25)
      denovo_dive_data_topleft <- filter(denovo_dive_data, MFN2>=0.25 & Amyloidbeta142<0.25)
      denovo_dive_data_bottomleft <- filter(denovo_dive_data, MFN2<0.25 & Amyloidbeta142<0.25)
      
      denovo_dive_data_percents <- c(dim(denovo_dive_data_topright)[1], dim(denovo_dive_data_bottomright)[1], dim(denovo_dive_data_topleft)[1], dim(denovo_dive_data_bottomleft)[1])
      denovo_dive_data_percents <- 100*(denovo_dive_data_percents / dim(denovo_dive_data)[1])
    
    #
    denovo_dive_data_subset_inverse <- filter(denovo_dive_data, MFN2==0 | PHF1Tau==0)
    denovo_dive_percents <- table(denovo_dive_subset$dnr2) / table(denovo_dive_data$dnr2) * 100
    denovo_dive_percents
    denovo_dive_percents <- as.data.frame(denovo_dive_percents)[1:5,] #grab just CA1, CA2 rows
    names(denovo_dive_percents) <- c("De novo regions", "Dbl_Pos_Percent")
    round(denovo_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> denovo_dive_percents["textPercent"] #convert to text percents and string for figures
    
    G5 <- filter(denovo_dive_data) %>%
      ggplot() +
      aes(x = Amyloidbeta142, y = MFN2, color = FlowSOM_ids) +
      geom_jitter(size = 3L, width = 0.01, height = 0.01) +
      scale_color_manual(labels = c("N","G","AP","TT","M"), values = colorKV1[5]) +
      geom_hline(yintercept=0.25, linetype="dashed", color = "blue", size=0.8) +
      geom_vline(xintercept=0.25, linetype="dashed", color = "blue", size=0.8) +
      #stat_smooth(data = filter(denovo_dive_data, Meta %in% c(1,2,3,4,8,10)), method = "gam", formula = y~exp(-x^2), color = "black", span = 11) +
      #geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.42,color="black",alpha=0) +
      #geom_jitter(data=denovo_dive_data_subset_inverse, size = 0.5L, width = 0.01, height = 0.01, color="gray") +
      xlim(-0.02,1) +
      ylim(-0.02,1) +
      coord_fixed(ratio=1) +
      theme_bw() +
      theme1a +
      labs(title="MFN2 & Amyloidbeta142 within Neurons", color="Cluster type") +
      #facet_grid(cols = vars(dnr2)) +
      
      #stat_cor(data=denovo_dive_data, label.x=0.05, label.y=1.1) + #calculate correlation of each dbl positive subpopulation
      #stat_regline_equation(data=denovo_dive_data, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y ~ log(x+1), label.x=0.05, label.y=1) + #caluclate and add information on corr, regression line to each subplot
      #geom_text(data=denovo_dive_percents, aes(x=0.9, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
      #theme(legend.position="none") +
      theme(panel.spacing = unit(2, "lines"))
    
    G5
    ggsave(filename = 'biaxial_MFN2_Abeta42_all_3pts_quads.pdf', plot = G5, width = 15, height = 25, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4/MFN2_biaxials'
    )
    
    #plot4 --> breaking down MFN2+ clusters into relationship between MFN2 and progressive disease markers across de novo regions
      denovo_dive_data <- fsom_neurons
      #denovo_dive_data[panel] <- log(denovo_dive_data[panel] + 1) # for any log transformations
      #percents
      denovo_dive_data <- filter(denovo_dive_data)
      denovo_dive_data$gate_level <- "NA"
      
      #PHF1Tau & MFN2
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2>=0.25 & PHF1Tau>=0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "TR"
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2<0.25 & PHF1Tau>=0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "BR"
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2>=0.25 & PHF1Tau<0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "TL"
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2<0.25 & PHF1Tau<0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "BL" 
      
      denovo_dive_data_top <- table(denovo_dive_data$gate_level, denovo_dive_data$de_novo_regions)
      denovo_dive_data_bottom <- table(denovo_dive_data$de_novo_regions)
      denovo_dive_data_percents <- apply(denovo_dive_data_top, 1, function(x) (x/denovo_dive_data_bottom)*100)
      
      #PHF1Tau & Abeta42
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2>=0.25 & Amyloidbeta142>=0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "TR"
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2<0.25 & Amyloidbeta142>=0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "BR"
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2>=0.25 & Amyloidbeta142<0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "TL"
      denovo_dive_data_topright = which(with(denovo_dive_data, MFN2<0.25 & Amyloidbeta142<0.25))
      denovo_dive_data[denovo_dive_data_topright, "gate_level"] = "BL" 
      
      denovo_dive_data_top <- table(denovo_dive_data$gate_level, denovo_dive_data$de_novo_regions)
      denovo_dive_data_bottom <- table(denovo_dive_data$de_novo_regions)
      denovo_dive_data_percents <- apply(denovo_dive_data_top, 1, function(x) (x/denovo_dive_data_bottom)*100)
      
      G5 <- filter(denovo_dive_data) %>%
        ggplot() +
        aes(x = Amyloidbeta142, y = MFN2, color = de_novo_regions) +
        geom_jitter(size = 1L, width = 0.01, height = 0.01) +
        scale_color_manual(labels = c("N","G","AP","TT","M"), values = color3) +
        geom_hline(yintercept=0.25, linetype="dashed", color = "blue", size=0.8) +
        geom_vline(xintercept=0.25, linetype="dashed", color = "blue", size=0.8) +
        #stat_smooth(data = filter(denovo_dive_data, Meta %in% c(1,2,3,4,8,10)), method = "gam", formula = y~exp(-x^2), color = "black", span = 11) +
        #geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.42,color="black",alpha=0) +
        #geom_jitter(data=denovo_dive_data_subset_inverse, size = 0.5L, width = 0.01, height = 0.01, color="gray") +
        xlim(-0.02,1) +
        ylim(-0.02,1) +
        coord_fixed(ratio=1) +
        theme_bw() +
        theme1a +
        labs(title="MFN2 & Abeta42 within Neurons", color="Cluster type") +
        facet_grid(cols = vars(dnr2)) +
        
        #stat_cor(data=denovo_dive_data, label.x=0.05, label.y=1.1) + #calculate correlation of each dbl positive subpopulation
        #stat_regline_equation(data=denovo_dive_data, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y ~ log(x+1), label.x=0.05, label.y=1) + #caluclate and add information on corr, regression line to each subplot
        #geom_text(data=denovo_dive_percents, aes(x=0.9, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
        #theme(legend.position="none") +
        theme(panel.spacing = unit(2, "lines"))
      
      G5
      ggsave(filename = 'biaxial_MFN2_Abeta42_all_lines_quads.pdf', plot = G5, width = 25, height = 12, dpi = 300, device = 'pdf',
             path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v4/MFN2_biaxials'
      )
      
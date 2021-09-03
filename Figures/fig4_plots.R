#### Figure Anat -> Identifying cell and disease object compositional differences in anatomically segregated human brain by AD status ####

#### A) Segmenting workflows: ez_segmenter and deepcell processes for identifying and extracting disease objects and cells, respectively ####  
#### B) Image overlays: DG->CA1 (use object_image_writer scripts, clustering on cell/object type, anatomical regions) ####
#### C1) Barplot: cell + object distribution in DG->CA1 anatomical regions ####
  
  master_gated_anat %>%
    dplyr::group_by(run_type_id, anat_id, FlowSOM_ids) %>%
    dplyr::summarise(n = n()) %>%
    rowwise() %>%
    mutate(anat_area_norm = n / select(filter(pixel_anat_area, Run == run_type_id & Mask == anat_id), area)[[1]]) -> mga_area_norm

  #all objects, area normalized  
  C1 <- ggplot(master_gated_anat) +
    aes(x = anat_id, fill = FlowSOM_ids) +
    geom_bar(position = "fill", width = 0.9) +
    scale_fill_manual(labels = c("tangles-threads", "amyloid plaques", "microglia", "endothelial", "neurons", "non-immune glia"), values = colorKV1) +
    theme_minimal() +
    labs(fill = "Cell / Object Type") +
    #scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
    labs(x = "Anatomical Region", y = "Proportion of cells or objects / region", title = "DG to CA1 cell and object composition") +
    theme1a +
    theme(axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.5)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    facet_wrap(vars(run_type_id), nrow = 1, scales = "fixed") +
    theme(panel.spacing = unit(2, "lines"))
  
  C1
  ggsave(filename = 'bar_composition_fill_DGtoCA1.pdf', plot = C1, width = 20, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  
  #all objects, area normalized
  C2 <- ggplot(mga_area_norm) +
    aes(x = anat_id, fill = FlowSOM_ids, weight = anat_area_norm) +
    geom_bar(position = "dodge", width = 0.9) + 
    scale_fill_manual(labels = c("tangles-threads", "amyloid plaques", "microglia", "endothelial", "neurons", "non-immune glia"), values = colorKV1) +
    scale_x_discrete(labels = c("DG","CA4","CA3","CA2","CA1")) +    theme_minimal() +
    labs(fill = "Cell / Object Type") +
    #scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
    labs(x = "Anatomical Region", y = bquote(bold("Total cells and disease objects "~(mm^2))), title = "Hippocampal band disease object composition") +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    facet_wrap(vars(run_type_id), nrow = 1, scales = "fixed")
  
  C2
  ggsave(filename = 'bar_all_objects_toArea_DGtoCA1.pdf', plot = C2, width = 20, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )

#### D) Barplot: disease object amounts in DG->CA1 anatomical regions ####
  # normalize cells & disease object by area of anatomical regions
  
  D <- ggplot(filter(mga_area_norm, FlowSOM_ids %in% c('tangles', 'plaques'))) +
    aes(x = anat_id, fill = FlowSOM_ids, weight = anat_area_norm) +
    geom_bar(position = "dodge", width = 0.9) + 
    scale_fill_manual(labels = c("tangles-threads", "amyloid plaques"), values = colorKV1[c(1,2)]) +
    theme_minimal() +
    labs(fill = "Disease object type") +
    #scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
    labs(x = "Anatomical Region", y = bquote("Disease objects "~(mm^2)), title = "Hippocampal band disease object composition") +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    facet_wrap(vars(run_type_id), nrow = 1, scales = "fixed") +
    theme(panel.spacing = unit(2, "lines"))
  
  D
  ggsave(filename = 'bar_disease_objects_toArea_DGtoCA1.pdf', plot = D, width = 20, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )

#### E) Scatterplot/Barplot/Dotplot: disease object overlap ratio in different cell types across anatomy, severity ####
  master_gated_anat2 %>%
    dplyr::count(run_type_id, anat_id, FlowSOM_ids) -> mga_simp
  master_gated_anat2 %>%
    dplyr::group_by(run_type_id, anat_id, FlowSOM_ids) %>%
    dplyr::summarise(ratio_t = sum(tau_obj_POS==1)/sum(tau_obj_POS==0), ratio_a = sum(amyloid_obj_POS==1)/sum(amyloid_obj_POS==0), n = n()) -> mga_simp
  
  # scatterplot of cell types and object positivity ratios
  E1 <- ggplot(mga_simp, aes(x = ratio_t, y = ratio_a)) +
    geom_point(aes(color = run_type_id, shape = FlowSOM_ids), size=10) +
    geom_point(aes(fill = run_type_id, shape = FlowSOM_ids), alpha=0.25, size=10) +
    scale_shape_manual(name = "Cell type", values=c(21,24,23,22)) +
    scale_color_manual(name = "Severity", values = color4) +
    scale_fill_manual(name = "Severity", values = color4) +
    #xlim(0,2.5) +
    #ylim(0,1.25) +
    scale_y_continuous(breaks = seq(0,1.5,0.5), limits = c(0,1.5)) +
    scale_x_continuous(breaks = seq(0,2.5,0.5), limits = c(0,2.5)) +
    coord_fixed(ratio=2) +
    theme_bw() +
    theme1a +
    labs(x="tangle-thread Positive Ratio", y="amyloid plaque Positive Ratio", title="DG to CA1: Ratio of cells in proximity to disease objects", color="Severity", shape="Cell type") +
    #theme1 +
    facet_grid(cols = vars(anat_id), scales = "fixed", space = "free_x") + 
    theme(panel.spacing.x = unit(2, "lines"))
  
  E1
  ggsave(filename = 'scatter_ratio_DiseasePOS_DGtoCA1_facetANATOMY_VERT_v2.pdf', plot = E1, width = 30, height = 15, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )#15,15 HORZ
  
  # scatterplot of cell types and object positivity ratios by severity
  E1b <- ggplot(mga_simp, aes(x = ratio_t, y = ratio_a)) +
    geom_point(aes(color = anat_id, shape = FlowSOM_ids), size=7, alpha=0.75) +
    scale_shape_manual(values=c(15,17,19,18)) +
    scale_color_manual(values = color2) +
    #xlim(0,2.5) +
    #ylim(0,1.25) +
    scale_y_continuous(breaks = seq(0,1.5,0.5)) +
    scale_x_continuous(breaks = seq(0,3,0.5)) +
    coord_fixed(ratio=2) +
    theme_bw() +
    theme1a +
    labs(x="tangle-thread Positive Ratio", y="amyloid plaque Positive Ratio", title="DG to CA1: Ratio of cells in proximity to disease objects", color="Severity", shape="Cell type") +
    #theme1 +
    facet_grid(rows = vars(run_type_id), scales = "fixed", space = "free_x")
  
  E1b
  ggsave(filename = 'scatter_ratio_DiseasePOS_DGtoCA1_facetSEVERITY_HORZ.pdf', plot = E1b, width = 15, height = 15, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  
  # barplot with just neurons or all cells for tangle ONLY
  E2t <- ggplot(filter(mga_simp, FlowSOM_ids == "neurons")) +
    aes(x = anat_id, fill = anat_id, weight = ratio_t) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = color2) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    labs(x="Neurons", y="tangle-thread positivity", title="DG to CA1: Ratio of cells in proximity to disease objects", fill="Anatomical region", shape="Cell type") +
    facet_grid(cols = vars(run_type_id))
  
  E2t
  ggsave(filename = 'bar_ratio_DiseasePOS_DGtoCA1_tauONLY_neuronsONLY.pdf', plot = E2t, width = 40, height = 40, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  # barplot with just neurons or all cells for amyloid ONLY
  E2a <- ggplot(filter(mga_simp, FlowSOM_ids == "neurons")) +
    aes(x = anat_id, fill = anat_id, weight = ratio_a) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = color2) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    labs(x="Neurons", y="amyloid plaque positivity", title="DG to CA1: Ratio of cells in proximity to disease objects", fill="Anatomical region", shape="Cell type") +
    facet_grid(cols = vars(run_type_id))
  
  E2a
  ggsave(filename = 'bar_ratio_DiseasePOS_DGtoCA1_amyloidONLY_neuronsONLY.pdf', plot = E2a, width = 40, height = 40, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )

  # dotplot - tangle-thread positivity
  E3t <- ggplot(mga_simp) +
    aes(x = anat_id, y = ratio_t, fill = FlowSOM_ids) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.75, dotsize=3) +
    scale_fill_manual(values = colorKV1[3:6]) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    scale_y_continuous(breaks = seq(0,2.0,0.5), limits=c(0,2.25)) +
    theme_bw() +
    theme1a +
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="Anatomical region", y="tangle-thread positivity", title="DG to CA1: Ratio of cells in proximity to disease objects", fill="Cell type") +
    facet_wrap(vars(run_type_id))
  
  E3t
  ggsave(filename = 'dotplot_ratio_DiseasePOS_DGtoCA1_tauONLY.pdf', plot = E3t, width = 15, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  # dotplot - plaque positivity
  E3a <- ggplot(mga_simp) +
    aes(x = anat_id, y = ratio_a, fill = FlowSOM_ids) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.75, dotsize=3) +
    scale_fill_manual(values = colorKV1[3:6]) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    scale_y_continuous(breaks = seq(0,2.0,0.5), limits=c(0,2.25)) +
    theme_bw() +
    theme1a +
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="Anatomical region", y="amyloid plaque positivity", title="DG to CA1: Ratio of cells in proximity to disease objects", fill="Cell type") +
    facet_wrap(vars(run_type_id))
  
  E3a
  ggsave(filename = 'dotplot_ratio_DiseasePOS_DGtoCA1_amyloidONLY.pdf', plot = E3a, width = 15, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  
#### F) Heatmap: show distribution of marker differences in different cell types between CA1 and CA2 across severity ####
  master_gated_anat2 <- inner_join(disease_cell_overlap, j_gates, by='spatial_id')
  #create vector for mean normalization and subtract from each cell
  single_mean_norm_vector <- apply(master_gated_anat2[,panel], 2, mean)
  master_gated_anat3 <- master_gated_anat2
  master_gated_anat3[,panel] <- master_gated_anat3[,panel] - single_mean_norm_vector
  #create dataframe of mean values for groups
  master_gated_anat3 %>%
    filter(anat_id == c('CA1','CA2')) %>%
    dplyr::group_by(run_type_id, anat_id, FlowSOM_ids) %>%
    dplyr::summarize_if(is.numeric, funs(mean)) %>%
    ungroup() ->
    heat_gated_regions
  #for heatmap colors - annotation
  mycolors <- list(
    run_type_id = c(Ctrl = "#46B2AC", MedAD = "#D1EF54", HiAD = "#E22D43"),
    anat_id = c(CA2 = "#CAB2D6", CA1 = "palegreen2"))
  #colors for consistent expression
    breaksList = seq(-0.1, .25, by = 0.0026) #Sets the minimum (-0.1), the maximum (0.25), and the increasing steps (.1) for the color scale
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)) # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
  
  #neurons
    heat_gated_regions_2 <- filter(heat_gated_regions, FlowSOM_ids == 'neurons')
    #organize naming, ordering schema
    hgr_panel = c('CD47','TotalTau','Calretinin','CD56Lyo','Parvalbumin','MAP2','PolyubiK48','MFN2','PanGAD6567','VGAT','PolyubiK63','Synaptophysin','VGLUT1','VGLUT2','Calbindin','PSD95')
    heat_gated_regions_2$anat_id = droplevels(heat_gated_regions_2$anat_id)
    heat_gated_regions_2 <- as.data.frame(heat_gated_regions_2)
    rownames(heat_gated_regions_2) <- c(1:nrow(heat_gated_regions_2))
    anno_data = heat_gated_regions_2[, c('run_type_id','anat_id')]
    #make heatmap
    F1 <- pheatmap(heat_gated_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Neurons", 
                   annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, gaps_row = c(2,4),
                   col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 10, height = 10,
                   filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3/hm_neurons_CA1CA2_meannorm_inferno.pdf')
    F1
    
  #non-immune glia
    heat_gated_regions_2 <- filter(heat_gated_regions, FlowSOM_ids == 'non-immune glia')
    #organize naming, ordering schema
    hgr_panel = c('GFAP','MBP','MAG','ApoE4')
    heat_gated_regions_2$anat_id = droplevels(heat_gated_regions_2$anat_id)
    heat_gated_regions_2 <- as.data.frame(heat_gated_regions_2)
    rownames(heat_gated_regions_2) <- c(1:nrow(heat_gated_regions_2))
    anno_data = heat_gated_regions_2[, c('run_type_id','anat_id')]
    #make heatmap
    F2 <- pheatmap(heat_gated_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Non-immune Glia", 
                   annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, gaps_row = c(2,4),
                   col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 10, height = 10, 
                   filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3/hm_nig_CA1CA2_meannorm_inferno.pdf')
    F2
    
  #endothelial
    heat_gated_regions_2 <- filter(heat_gated_regions, FlowSOM_ids == 'endothelial')
    #organize naming, ordering schema
    hgr_panel = c('CD31','CD105','MCT1')
    heat_gated_regions_2$anat_id = droplevels(heat_gated_regions_2$anat_id)
    heat_gated_regions_2 <- as.data.frame(heat_gated_regions_2)
    rownames(heat_gated_regions_2) <- c(1:nrow(heat_gated_regions_2))
    anno_data = heat_gated_regions_2[, c('run_type_id','anat_id')]
    #make heatmap
    F3 <- pheatmap(heat_gated_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Endothelial", 
                   annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, gaps_row = c(2,4),
                   col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 10, height = 10,
                   filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3/hm_endo_CA1CA2_meannorm_inferno.pdf')
    F3
    
  #microglia
    heat_gated_regions_2 <- filter(heat_gated_regions, FlowSOM_ids == 'microglia')
    #organize naming, ordering, coloring schema  
    hgr_panel = setdiff(c(glia_struct3), c('X8OHGuano','pTDP43','EEA1','CD47'))
    heat_gated_regions_2$anat_id = droplevels(heat_gated_regions_2$anat_id)
    heat_gated_regions_2 <- as.data.frame(heat_gated_regions_2)
    rownames(heat_gated_regions_2) <- c(1:nrow(heat_gated_regions_2))
    anno_data = heat_gated_regions_2[, c('run_type_id','anat_id')]
    #make heatmap  
    F4 <- pheatmap(heat_gated_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Microglia", 
             annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, gaps_row = c(2,4),
             col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 10, height = 10,
             filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3/hm_mgla_CA1CA2_meannorm_inferno.pdf')
    F4

#### G1,2) Data Wrangling ####
    hgr_panel_n = c('CD47','TotalTau','Calretinin','CD56Lyo','Parvalbumin','MAP2','PolyubiK48','MFN2','PanGAD6567','VGAT','PolyubiK63','Synaptophysin','VGLUT1','VGLUT2','Calbindin','PSD95')
    hgr_panel_nig = c('GFAP','MBP','MAG','ApoE4')
    hgr_panel_e = c('CD31','CD105','MCT1')
    hgr_panel_m = setdiff(c(glia_struct3), c('X8OHGuano','pTDP43','EEA1','CD47'))
    hgr_panel_neutral = c('C12', 'HistonerH3Lyo')
    
    # NOTE: can adjust data to account for unique severity (Hi, Med, Ctrl) by including run_type_id in the initial filtering down from master_gated_anat2 to mga2_diseaseObj_POS
    
    #### Tau objects ####
      master_gated_anat2 %>%
        filter(!(amyloid_obj_POS == 1)) %>% #filters out any cells having JUST amyloid plaques.
        dplyr::group_by(FlowSOM_ids, anat_id, tau_obj_POS, run_type_id) %>%
        summarize_if(is.numeric, funs(mean)) -> mga2_tauPOS
      
      mga2_tauPOS_CA1 <- filter(mga2_tauPOS, tau_obj_POS %in% c(1,0), anat_id %in% c('CA1'))
      mga2_tauPOS_CA2 <- filter(mga2_tauPOS, tau_obj_POS %in% c(1,0), anat_id %in% c('CA2'))
      mga2_tauPOS_1 <- filter(mga2_tauPOS, tau_obj_POS %in% c(1), anat_id %in% c('CA1','CA2'))
      mga2_tauPOS_0 <- filter(mga2_tauPOS, tau_obj_POS %in% c(0), anat_id %in% c('CA1','CA2'))
      
      mga2_tauPOS_markers <- mga2_tauPOS_CA2[,markers] / mga2_tauPOS_CA1[,markers]
      mga2_tauPOS_markers <- mga2_tauPOS_1[,markers] / mga2_tauPOS_0[,markers]
      #use if computing CA2 / CA1 ratio
      mga2_tauPOS_markers$run_type_id <- mga2_tauPOS_CA1$run_type_id
      mga2_tauPOS_markers$tau_obj_POS <- mga2_tauPOS_CA1$tau_obj_POS
      mga2_tauPOS_markers$FlowSOM_ids <- mga2_tauPOS_CA1$FlowSOM_ids
      #use if computing disease object + / - ratio
      mga2_tauPOS_markers$run_type_id <- mga2_tauPOS_1$run_type_id
      mga2_tauPOS_markers$anat_id <- mga2_tauPOS_1$anat_id
      mga2_tauPOS_markers$FlowSOM_ids <- mga2_tauPOS_1$FlowSOM_ids
      
    #### Amyloid objects ####
      master_gated_anat2 %>%
        filter(!(tau_obj_POS == 1 )) %>% #filters out any cells having JUST tangles/threads.
        dplyr::group_by(FlowSOM_ids, anat_id, amyloid_obj_POS, run_type_id) %>%
        summarize_if(is.numeric, funs(mean)) -> mga2_amyPOS  
      
      mga2_amyPOS_CA1 <- filter(mga2_amyPOS, amyloid_obj_POS %in% c(1,0), anat_id %in% c('CA1'))
      mga2_amyPOS_CA2 <- filter(mga2_amyPOS, amyloid_obj_POS %in% c(1,0), anat_id %in% c('CA2'))
      mga2_amyPOS_1 <- filter(mga2_amyPOS, amyloid_obj_POS %in% c(1), anat_id %in% c('CA1','CA2'))
      mga2_amyPOS_0 <- filter(mga2_amyPOS, amyloid_obj_POS %in% c(0), anat_id %in% c('CA1','CA2'))
      
      mga2_amyPOS_markers <- mga2_amyPOS_CA2[,markers] / mga2_amyPOS_CA1[,markers]
      mga2_amyPOS_markers <- mga2_amyPOS_1[,markers] / mga2_amyPOS_0[,markers]
      #use if computing CA2 / CA1 ratio
      mga2_amyPOS_markers$run_type_id <- mga2_amyPOS_CA1$run_type_id
      mga2_amyPOS_markers$amyloid_obj_POS <- mga2_amyPOS_CA1$amyloid_obj_POS
      mga2_amyPOS_markers$FlowSOM_ids <- mga2_amyPOS_CA1$FlowSOM_ids
      #use if computing disease object + / - ratio
      mga2_amyPOS_markers$run_type_id <- mga2_amyPOS_1$run_type_id
      mga2_amyPOS_markers$anat_id <- mga2_amyPOS_1$anat_id
      mga2_amyPOS_markers$FlowSOM_ids <- mga2_amyPOS_1$FlowSOM_ids
      
#### G1) Lineplot: ####
  
  # Compare tau+ and tau- object mean expressions across disease severity within CA2 / CA1 ratio
  filtered_data <- filter(mga2_tauPOS_markers, FlowSOM_ids == 'microglia')
  filtered_data <- filtered_data %>% mutate(group_name = paste0(tau_obj_POS,'_'))
   
  ggplot(filtered_data, aes(x = run_type_id, y = Iba1, group = tau_obj_POS)) +
    geom_line(aes(color = tau_obj_POS)) +
    geom_point(aes(color = tau_obj_POS)) +
    geom_hline(yintercept=1, linetype="dashed", color="black") +
    ylim(0,2) +
    theme_minimal() 
  
  # Compare disease object mean expressions across anatomical regions 
  filtered_data <- filter(mga2_tauPOS, FlowSOM_ids == 'neurons', anat_id %in% c('CA1','CA2'))
  filtered_data <- filtered_data %>% mutate(group_name = paste0(tau_obj_POS,'_'))
  
  ggplot(filtered_data, aes(x = anat_id, y = VGAT, group = tau_obj_POS)) +
    geom_boxplot(aes(color = tau_obj_POS)) +
    theme_minimal() 
 
#### G2) *Ratio-plot: ####
  #### Tau objects ####
  # Calculate CA2 to CA1 ratio of cell marker expression means, compare differences of this metric in cells with disease objects and those without (each sample a data point)
  filter(mga2_tauPOS_markers, FlowSOM_ids == 'neurons') %>% select(c(hgr_panel_n,'anat_id','FlowSOM_ids')) -> filtered_data
  gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(filtered_data), c(hgr_panel_n))))
  gathered_data$expr <- log(gathered_data$expr, base = 2)
  
  G1 <- ggplot(gathered_data, aes(x = expr, y = markers)) +
    geom_point(aes(color = anat_id, shape = run_type_id, alpha = 0.85), size = 10) +
    scale_color_manual(labels = c("CA2", "CA1"), values = color2[c(4,5)]) +
    geom_vline(xintercept=0, linetype='dashed', color="black") +
    labs(x="Mean tau-diseased neuron expression / Mean healthy neuron expression (log2)", y="Markers", title="Differences of neuronal marker expression in tau-disease", color="Anatomical region", shape="Sample") +
    theme_bw() +
    theme1a +
    facet_wrap(vars(anat_id))
  
  G1
  ggsave(filename = 'scatter_ratio_CA2CA1_tauPOS_neuron_markers_expr_anat.pdf', plot = F1, width = 20, height = 20, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/paper1_final/Fig4/v2/'
  )
  
  #### Amyloid objects ####
  # Calculate CA2 to CA1 ratio of cell marker expression means, compare differences of this metric in cells with disease objects and those without (each sample a data point)
  filter(mga2_amyPOS_markers, FlowSOM_ids == 'neurons') %>% select(c(hgr_panel_n,'run_type_id','anat_id','FlowSOM_ids')) -> filtered_data
  gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(filtered_data), c(hgr_panel_n))))
  gathered_data$expr <- log(gathered_data$expr, base = 2)
  
  G2 <- ggplot(gathered_data, aes(x = expr, y = markers)) +
    geom_point(aes(color = anat_id, shape = run_type_id, alpha = 0.85), size = 10) +
    scale_color_manual(labels = c("CA2", "CA1"), values = color2[c(4,5)]) +
    geom_vline(xintercept=0, linetype='dashed', color="black") +
    labs(x="Mean amyloid-diseased neuron expression / Mean healthy neuron expression (log2)", y="Markers", title="Differences of neuronal marker expression in amyloid-disease", color="Anatomical region", shape="Sample") +
    theme_bw() +
    theme1a +
    facet_wrap(vars(anat_id))
  
  G2
  ggsave(filename = 'scatter_ratio_CA2CA1_amyloidPOS_neuron_markers_expr_anat2.pdf', plot = F2, width = 20, height = 20, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/paper1_final/Fig4/v2/'
  )
  
  
#### G3) Boxplot: #### 
  protect_markers <- c("PolyubiK48","PolyubiK63","MFN2")
  # more specific distro of markers - neurons with tangles vs those without in selected markers between CA1 and CA2 across severity (amyloid in supplement)
  #DG->CA1: 
  filtered_data <- filter(master_gated_anat2, anat_id %in% c('CA1','CA2') & FlowSOM_ids == 'neurons' )
  gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(master_gated_anat2), markers)))
  G3box <- ggplot(filter(gathered_data, markers %in% c(protect_markers, "PHF1Tau","Calbindin","Calretinin","HistoneH3Lyo"))) +
    aes(x = anat_id, y = expr, fill = amyloid_obj_POS) +
    geom_violin(adjust = 1L, scale = "area") +
    #scale_y_continuous(trans = "log") +
    scale_fill_hue() +
    labs(title = "", subtitle = "") +
    theme_minimal() +
    labs(fill = "Sample Type") +
    labs(x = "Region", y = "Expression", title = "") +
    theme1a +
    facet_grid(rows = vars(markers), cols = vars(run_type_id), scales = "free_y")
  
  G3box
  ggsave(filename = 'box_ratio_tauPOS_neuron_markers_expr_FACET_anat_severity.pdf', plot = G3box, width = 40, height = 40, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  
  # more specific distro of markers - microglia with tangles vs those without in selected markers between CA1 and CA2 across severity (amyloid in supplement)
  #DG->CA1: 
  filtered_data <- filter(master_gated_anat2, anat_id %in% c('CA1','CA2') & FlowSOM_ids == 'microglia' )
  gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(master_gated_anat2), nu_panel)))
  ggplot(filter(gathered_data, markers %in% c('PHF1Tau','Iba1','CD47','ApoE4'))) +
    aes(x = run_type_id, y = expr, fill = anat_id, color = tau_obj_POS) +
    geom_boxplot(adjust = 1L, scale = "area") +
    #scale_y_continuous(trans = "log") +
    scale_fill_hue() +
    labs(title = "", subtitle = "") +
    theme_minimal() +
    labs(fill = "Sample Type") +
    labs(x = "Region", y = "Expression", title = "") +
    theme1 +
    facet_wrap(vars(markers), scales = "free_y")

#### G4) Barplot w/ Error bars ####
  #### Tau objects ####
  # Calculate CA2 to CA1 ratio of cell marker expression means, compare differences of this metric in cells with disease objects and those without (each sample a data point)
  filter(mga2_tauPOS_markers, FlowSOM_ids == 'neurons') %>% select(c(hgr_panel_n,'HistoneH3Lyo','anat_id','FlowSOM_ids','run_type_id')) -> filtered_data
  gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(filtered_data), c(hgr_panel_n,'HistoneH3Lyo'))))
  gathered_data$expr <- log(gathered_data$expr, base = 2)
  
  G1 <- ggplot(gathered_data) +
    aes(y = markers, weight = expr, fill = anat_id) +
    geom_bar(position = "dodge") +
    scale_fill_manual(labels = c("CA2", "CA1"), values = color2[c(4,5)]) +
    labs(x="Mean tau-diseased neuron expression / Mean healthy neuron expression (log2)", y="Markers", title="Regional differences of tau-diseased neuronal marker expression", fill="Anatomical region") +
    theme_minimal() +
    theme1a +
    facet_grid(rows = vars(anat_id), cols = vars(run_type_id))
  
  G1
  ggsave(filename = 'bar_ratio_tauPOS_neuron_markers_expr_FACET_anat_severity.pdf', plot = G1, width = 20, height = 20, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
  
  #### Amyloid objects ####
  # Calculate CA2 to CA1 ratio of cell marker expression means, compare differences of this metric in cells with disease objects and those without (each sample a data point)
  filter(mga2_amyPOS_markers, FlowSOM_ids == 'neurons') %>% select(c(hgr_panel_n,'HistoneH3Lyo','anat_id','FlowSOM_ids','run_type_id')) -> filtered_data
  gathered_data <- gather(filtered_data, key = "markers", value = "expr", -c(setdiff(names(filtered_data), c(hgr_panel_n, 'HistoneH3Lyo'))))
  gathered_data$expr <- log(gathered_data$expr, base = 2)
  
  G2 <- ggplot(gathered_data) +
    aes(y = markers, weight = expr, fill = anat_id) +
    geom_bar(position = "dodge") +
    scale_fill_manual(labels = c("CA2", "CA1"), values = color2[c(4,5)]) +
    labs(x="Mean amyloid-diseased neuron expression / Mean healthy neuron expression (log2)", y="Markers", title="Regional differences of amyloid-diseased neuronal marker expression", fill="Anatomical region") +
    theme_minimal() +
    theme1a +
    facet_grid(rows = vars(anat_id), cols = vars(run_type_id))
  
  G2
  ggsave(filename = 'bar_ratio_amyloidPOS_neuron_markers_expr_FACET_anat_severity.pdf', plot = F2, width = 20, height = 20, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
  )
#### G5) Biaxials #####
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
    anat_dive_data <- master_gated_anat2 %>% filter(FlowSOM_ids %in% "neurons" & anat_id %in% c("CA1","CA2"))
    anat_dive_subset <- filter(anat_dive_data, PolyubiK48>0 & PHF1Tau>0)
    anat_dive_subset_inverse <- filter(anat_dive_data, PolyubiK48==0 | PHF1Tau==0)
    anat_dive_percents <- table(anat_dive_subset$run_type_id, anat_dive_subset$anat_id) / table(anat_dive_data$run_type_id, anat_dive_data$anat_id) * 100
    anat_dive_percents
    anat_dive_percents <- as.data.frame(anat_dive_percents)[10:15,] #grab just CA1, CA2 rows
    names(anat_dive_percents) <- c("run_type_id", "anat_id", "Dbl_Pos_Percent")
    round(anat_dive_percents["Dbl_Pos_Percent"], digits=2) %>% mutate(textPercent = paste0(Dbl_Pos_Percent, "%")) %>% select(textPercent) -> anat_dive_percents["textPercent"] #convert to text percents and string for figures
    #plot
    G5 <- anat_dive_data %>%
          ggplot() +
          aes(x = PHF1Tau, y = PolyubiK48, color = run_type_id) +
          geom_point(size = 1L) +
          scale_color_manual(values = color4) +
          stat_smooth(data = filter(anat_dive_subset, run_type_id != "Ctrl"), method='lm', formula= y~x, color = "black") +
          geom_rect(xmin=0.03,xmax=1.16,ymin=0.01,ymax=0.42,color="black",alpha=0) +
          geom_point(data=anat_dive_subset_inverse, size=1L, color="gray") +
          theme_bw() +
          theme1a +
          labs(title="Polyubiquitins & Disease within Neurons", color="Severity") +
          facet_grid(rows=vars(anat_id), cols=vars(run_type_id)) +
          stat_cor(data=anat_dive_subset, label.x=0.3, label.y=0.4) + #calculate correlation of each dbl positive subpopulation
          stat_regline_equation(data=anat_dive_subset, aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y~x, label.x=0.3, label.y=0.38) + #caluclate and add information on corr, regression line to each subplot
          geom_text(data=anat_dive_percents, aes(x=0.9, y=0.35, label=textPercent, size=20), inherit.aes = FALSE) + #add percent dbl positive to each subplot
          theme(legend.position="none")
    
    G5
    ggsave(filename = 'biaxial_Pk48_PHF1Tau_DGtoCA1.pdf', plot = G5, width = 15, height = 10, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3'
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
    
    
    
#### -- -- -- Replacement UMAPs --> Microglia UMAPs with marker overlays (instead of heatmaops) -- -- -- ####
#### F) Microglia UMAPs ####
    #### A) UMAP of subclustered neurons ####  
    # Severity  
    p1_U <- filter(microglia_umap_data) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color=run_type_id)) +
      geom_density_2d(color="#CCCCCC") +
      geom_point_rast(size = 0.2) + # original size = 0.5
      scale_color_manual(values = color4) +
      coord_fixed(ratio = 1) +
      xlim(-20,30) +
      ylim(-20,30) +
      labs(x="UMAP1", y="UMAP2", title="Microglia", color="") +
      theme_bw() + 
      theme1a +
      guides(colour = guide_legend(override.aes = list(size=4)))
    #facet_wrap(vars(markers))
    p1_U
    ggsave(filename = 'UMAP_neuron_subclusters_Severity_blank_pt3.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigSubNeurons/v2'
    )
    # markers
    for (run in c("Ctrl","MedAD","HiAD")) {
      for (marker in c("Iba1","CD45","ApoE4","CD33Lyo")) {
        p1_U <- filter(subsetted_norm_data, run_type_id == run) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
          ggplot(aes(x = UMAP1, y = UMAP2, color=eval(parse(text=marker)))) +
          geom_density_2d(data = subsetted_norm_data, color="#CCCCCC") +
          geom_point_rast(size = 5) + # original size = 0.5
          scale_color_gradientn(colours = inferno(100)) +
          coord_fixed(ratio = 1) +
          xlim(-20,30) +
          ylim(-20,30) +
          labs(x="UMAP1", y="UMAP2", title="Microglia", color="") +
          theme_bw() + 
          theme4 + #4
          guides(colour = guide_legend(override.aes = list(size=4))) +
          #facet_grid(cols = vars(run_type_id)) +
          guides(fill = guide_colourbar())
        p1_U
        ggsave(filename = paste0('UMAP_microglia_',run,'_',marker,'_severity_DGtoCA1_ONLY_pt5_v2.pdf'), plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
               path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigAnat/v3/UMAPs/separate_microglia_UMAPs_v2'
        )
      }
    }
    
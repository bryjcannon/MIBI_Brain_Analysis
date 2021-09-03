#### Figure DeNovo -> Identifying cell and disease object compositional differences in data-driven segregated human brain by AD status ####

#### resetting de_novo_regions order ####
master_regions_c$de_novo_regions <- factor(master_regions_c$de_novo_regions, levels = c(1,2,5,3,4))
disease_cell_overlap$de_novo_regions <- factor(disease_cell_overlap$de_novo_regions, levels = c(1,2,5,3,4))
cell_disease_overlap$de_novo_regions <- factor(cell_disease_overlap$de_novo_regions, levels = c(1,2,5,3,4))
syn_percent$`De novo region` <- factor(syn_percent$`De novo region`, levels = c("ND","GD","MD","PD","TD"))

#### A) Segmenting workflows: ez_segmenter and deepcell processes for identifying and extracting disease objects and cells, respectively ####  
#### B) Image overlays: DG->CA1 (use object_image_writer scripts, clustering on cell/object type, anatomical regions) ####
#### C) Barplot: region distribution in entirety of tissue samples ####

  #cell/obj distribution of each region as percentage of total cells in region
  C1 <- ggplot(master_regions_c) +
    aes(x = run_type_id, fill = de_novo_regions) +
    geom_bar(position = "fill", width = 0.8) +
    scale_fill_manual(name = "Region", labels = dnr_names, values = color3) +
    theme_minimal() +
    labs(fill = "Region") +
    scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
    labs(x = "Severity", y = "Regional proportion", title = "Hippocampus - de novo region composition") +
    theme1a +
    theme(axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.5)) +#rotates yaxis label
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #gets rid of lines
    
  C1
  ggsave(filename = 'bar_composition_fill_denovoregions.pdf', plot = C1, width = 10, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  
#### D) Barplot: cell + object distribution in regions ####
  # proportional cells and objects, by severity and de_novo_region
  D1 <- ggplot(master_regions_c) +
    aes(x = run_type_id, fill = FlowSOM_ids) +
    geom_bar(position = "fill", width = 0.9) +
    scale_fill_manual(labels = c("tangles-threads", "amyloid plaques", "microglia", "endothelial", "neurons", "non-immune glia"), values = colorKV1) +
    theme_minimal() +
    labs(fill = "Cell / Object Type") +
    scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
    labs(x = "Severity", y = "Proportion of cells or objects / region", title = "Cell and object composition by Region") +
    theme1a +
    theme(axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.5)) +#rotates yaxis label
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    facet_wrap(vars(de_novo_regions), nrow = 1, scales = "fixed", labeller = labeller(de_novo_regions = dnr_names)) +
    theme(panel.spacing = unit(2, "lines"))
  
  D1
  ggsave(filename = 'bar_composition_fill_cells_objects_denovoregions_fill.pdf', plot = D1, width = 30, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )

  # total cells and objects, by severity and de_novo_region
  D2 <- ggplot(master_regions_c) +
    aes(x = run_type_id, fill = FlowSOM_ids) +
    geom_bar(position = position_dodge2(width = 3, preserve = "single")) +
    scale_fill_manual(labels = c("tangles-threads", "amyloid plaques", "microglia", "endothelial", "neurons", "non-immune glia"), values = colorKV1) +
    theme_minimal() +
    labs(fill = "Cell / Object Type") +
    scale_x_discrete(limits=c("Ctrl", "MedAD", "HiAD")) +
    theme_minimal() +
    labs(x = "Severity", y = "Total cells or objects / region", title = "Cell and object composition by Region") +
    theme1a +
    theme(axis.title.y = element_text(angle = 90, vjust = 1, hjust=0.5)) +#rotates yaxis label
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=20, angle=45, hjust=1)) + #gets rid of lines
    facet_wrap(vars(de_novo_regions), nrow = 1, scales = "fixed", labeller = labeller(de_novo_regions = dnr_names))
  
  D2
  ggsave(filename = 'bar_composition_fill_cells_objects_denovoregions_dodge.pdf', plot = D2, width = 20, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  
#### E) Scatterplot: cell + disease object overlap distribution in de novo regions (by regions) ####
  disease_cell_overlap %>%
    dplyr::group_by(run_type_id, de_novo_regions, FlowSOM_ids) %>%
    dplyr::summarise(ratio_t = sum(tau_obj_POS==1)/sum(tau_obj_POS==0), ratio_a = sum(amyloid_obj_POS==1)/sum(amyloid_obj_POS==0), n = n()) -> dco2
  
  #scatterplot original
  E <- ggplot(dco2, aes(x = ratio_t, y = ratio_a)) +
    geom_point(aes(color = run_type_id, shape = FlowSOM_ids), size=10) +
    geom_point(aes(fill = run_type_id, shape = FlowSOM_ids), alpha=0.25, size=10) +
    scale_shape_manual(name = "Cell type", values=c(21,24,23,22)) +
    scale_color_manual(name = "Severity", values = color4) +
    scale_fill_manual(name = "Severity", values = color4) +
    scale_y_continuous(breaks = seq(0,3,0.5), limits = c(0,3)) +
    scale_x_continuous(breaks = seq(0,2,0.5), limits = c(0,2)) +
    coord_fixed(ratio=2) +
    theme_bw() +
    theme1a +
    labs(x="tangle-thread Positive Ratio", y="amyloid plaque Positive Ratio", title="De novo regions: Ratio of cells in proximity to disease objects", color="De novo Region", shape="Cell type") +
    facet_grid(cols = vars(de_novo_regions), scales = "fixed", space = "free_x", labeller = labeller(de_novo_regions = dnr_names)) +
    theme(panel.spacing.x = unit(2, "lines"))
  
  E
  ggsave(filename = 'scatter_ratio_DiseasePOS_denovoregions_fixed_facetREGIONS_HORZ.pdf', plot = E, width = 20, height = 20, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )#15,30 for VERT

  
  # dotplot - tangle-thread positivity
  E3t <- ggplot(dco2) +
    aes(x = de_novo_regions, y = ratio_t, fill = FlowSOM_ids) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.75, dotsize=3) +
    scale_fill_manual(values = colorKV1[3:6]) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    scale_y_continuous(breaks = seq(0,3.0,0.5), limits=c(0,3.0)) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="De novo Region", y="tangle-thread positivity", title="De novo Region: Ratio of cells in proximity to disease objects", fill="Cell type") +
    facet_wrap(vars(run_type_id))
  
  E3t
  ggsave(filename = 'dotplot_ratio_DiseasePOS_denovoregions_tauONLY_square5.pdf', plot = E3t, width = 15, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  ) #w=20,h=5 in v1
  # dotplot - plaque positivity
  E3a <- ggplot(dco2) +
    aes(x = de_novo_regions, y = ratio_a, fill = FlowSOM_ids) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.75, dotsize=3) +
    scale_fill_manual(values = colorKV1[3:6]) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    scale_y_continuous(breaks = seq(0,3.0,0.5), limits=c(0,3.0)) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="De novo Region", y="amyloid plaque positivity", title="De novo Region: Ratio of cells in proximity to disease objects", fill="Cell type") +
    facet_wrap(vars(run_type_id))
  
  E3a
  ggsave(filename = 'dotplot_ratio_DiseasePOS_denovoregions_amyloidONLY_square5.pdf', plot = E3a, width = 15, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )#w=20,h=5 in v1
  
  # barplot with just neurons or all cells for tangle ONLY
  E2t <- ggplot(disease_cell_overlap) +
    aes(x = de_novo_regions, fill = tau_obj_POS) +
    geom_bar(position = "dodge") +
    scale_fill_manual(labels=c("No", "Yes"), values = color5) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    labs(x="De novo Region", y="Cells", title="De novo Region: Cells in proximity to disease objects", fill="Tangle-thread associated?") +
    facet_grid(rows = vars(FlowSOM_ids), cols = vars(run_type_id))
  
  E2t
  ggsave(filename = 'bar_counts_DiseasePOS_denovorregions_tau.pdf', plot = E2t, width = 20, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  # barplot with just neurons or all cells for amyloid ONLY
  E2a <- ggplot(disease_cell_overlap) +
    aes(x = de_novo_regions, fill = amyloid_obj_POS) +
    geom_bar(position = "dodge") +
    scale_fill_manual(labels=c("No", "Yes"), values = color5) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    labs(x="De novo Region", y="Cells", title="De novo Region: Cells in proximity to disease objects", fill="Amyloid plaque associated?") +
    facet_grid(rows = vars(FlowSOM_ids), cols = vars(run_type_id))
  
  E2a
  ggsave(filename = 'bar_counts_DiseasePOS_denovoregions_amyloid.pdf', plot = E2a, width = 20, height = 10, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  
  
#### F) Scatterplot: disease object overlap (in/out of cells) distribution in de novo regions (by regions) ####
  cell_disease_overlap %>%
    dplyr::group_by(run_type_id, de_novo_regions, FlowSOM_ids) %>%
    dplyr::summarise(ratio_ct = sum(cells_tau_obj_POS==1)/sum(cells_tau_obj_POS==0), ratio_ca = sum(cells_amyloid_obj_POS==1)/sum(cells_amyloid_obj_POS==0), n_c = n()) -> cdo2
  
  # dotplot - tangle-thread positivity
  E3t <- filter(cdo2, FlowSOM_ids == "tangles") %>% ggplot() +
    aes(x = run_type_id, y = ratio_ct, fill = de_novo_regions) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.75, dotsize=3) +
    scale_fill_manual(values = color3) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    scale_y_continuous(breaks = seq(0,3.0,0.5), limits=c(0,3.0)) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="De novo Region", y="tangle-thread positivity", title="De novo Region: Ratio of cells in proximity to disease objects", fill="Cell type") +
    facet_wrap(vars(run_type_id))
  
  E3t
  ggsave(filename = 'dotplot_ratio_DiseasePOS_denovoregions_tauONLY.pdf', plot = E3t, width = 20, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  # dotplot - plaque positivity
  E3a <- ggplot(dco2) +
    aes(x = de_novo_regions, y = ratio_ca, fill = FlowSOM_ids) +
    geom_dotplot(binaxis='y', stackdir='center', alpha=0.75, dotsize=3) +
    scale_fill_manual(values = colorKV1[3:6]) +
    geom_hline(yintercept=1.0, linetype="dashed", color = "black") +
    scale_y_continuous(breaks = seq(0,3.0,0.5), limits=c(0,3.0)) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="De novo Region", y="amyloid plaque positivity", title="De novo Region: Ratio of cells in proximity to disease objects", fill="Cell type") +
    facet_wrap(vars(run_type_id))
  
  E3a
  ggsave(filename = 'dotplot_ratio_DiseasePOS_denovoregions_amyloidONLY.pdf', plot = E3a, width = 20, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  
  # barplot with just neurons or all cells for tangle ONLY
  E2t <- ggplot(filter(cell_disease_overlap, FlowSOM_ids == "tangles")) +
    aes(x = de_novo_regions, fill = cells_tau_obj_POS) +
    geom_bar(position = "dodge") +
    scale_fill_manual(labels=c("No", "Yes"), values = color5) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    labs(x="De novo Region", y="tangle-threads", title="De novo Region: Disease objects by cell soma association", fill="Cell soma associated?") +
    facet_grid(cols = vars(run_type_id))
  
  E2t
  ggsave(filename = 'bar_counts_DiseasePOS_denovorregions_tauSOMA.pdf', plot = E2t, width = 20, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  # barplot with just neurons or all cells for amyloid ONLY
  E2a <- ggplot(filter(cell_disease_overlap, FlowSOM_ids == "plaques")) +
    aes(x = de_novo_regions, fill = cells_amyloid_obj_POS) +
    geom_bar(position = "dodge") +
    scale_fill_manual(labels=c("No", "Yes"), values = color5) +
    scale_x_discrete(labels=dnr_names) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #gets rid of lines
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(face="bold", size=14, angle=45, hjust=1)) + #gets rid of lines
    labs(x="De novo Region", y="amyloid plaques", title="De novo Region: Disease objects by cell soma association", fill="Cell soma associated?") +
    facet_grid(cols = vars(run_type_id))
  
  E2a
  ggsave(filename = 'bar_counts_DiseasePOS_denovoregions_amyloidSOMA.pdf', plot = E2a, width = 20, height = 5, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
  
#### G) UMAP: Figure Overlays - Regions, Cell Types, & Callouts ####
  #data to be used for umaps
  mod_subsampled <- filter(obj_umap_plot_all_data, UMAP1 > -20 & UMAP2 < 10)
  
  #### general ####
  umap_clustering_markers <-  c("CD47","TotalTau","ApoE4","Amyloidbeta140","Amyloidbeta142","Calretinin",    
                                "Parvalbumin","PanAmyloidbeta1724","CD31","MAP2","PolyubiK48","MFN2",           
                                "PanGAD6567","VGAT", "PolyubiK63","pTDP43","Synaptophysin","CD45",            
                                "PanApoE2E3E4","X8OHGuano","Presenilin1NTF","MAG","VGLUT1","VGLUT2",            
                                "Calbindin","Iba1","PSD95","MBP","MCT1","CD33Lyo",           
                                "GFAP","CD105","PHF1Tau")
  # umap cell + object overlays: De Novo Regions
    p1_U <-
      ggplot(mod_subsampled, aes(x = UMAP1, y = UMAP2)) +
      geom_density2d(color="#CCCCCC") +
      geom_point_rast(size = 5, stroke = 0, aes(color=de_novo_regions)) + # original size = 0.5
      coord_fixed(ratio = 1) + 
      
      xlim(min(mod_subsampled$UMAP1)-1, max(mod_subsampled$UMAP1)+1) +
      ylim(min(mod_subsampled$UMAP2)-1.5, max(mod_subsampled$UMAP2)+1.5) +
      scale_color_manual(values=color3) + 
      labs(x="UMAP1", y="UMAP2", title="De Novo Spatial Regions", color="") +
      theme_bw() + 
      theme3 +
      theme(legend.key.size = unit(2, 'cm'))
       #+ 
    #guides(color = guide_legend(override.aes = list(size = 5))) #exclude this line if using empty template
    p1_U
    ggsave(filename = 'UMAP_denovoregions_blank_pt5.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
    )
    
  # umap cell + object overlays: Cell and object types
    p1_U <-
      ggplot(mod_subsampled, aes(x = UMAP1, y = UMAP2)) +
      geom_density2d(color="#CCCCCC") +
      geom_point_rast(size = 3, stroke = 0, aes(color=FlowSOM_ids)) + # original size = 0.5
      coord_fixed(ratio = 1) + 
      
      xlim(min(mod_subsampled$UMAP1)-1, max(mod_subsampled$UMAP1)+1) +
      ylim(min(mod_subsampled$UMAP2)-1.5, max(mod_subsampled$UMAP2)+1.5) +
      scale_color_manual(values=colorKV1) + 
      labs(x="UMAP1", y="UMAP2", title="De Novo Spatial Regions", color="") +
      theme_bw() + 
      theme3 +
      theme(legend.key.size = unit(2, 'cm'))
    #+ 
    #guides(color = guide_legend(override.aes = list(size = 5))) #exclude this line if using empty template
    p1_U
    ggsave(filename = 'UMAP_cells_disease_objects_blank_pt3.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
    )
    
  # umap cell + object overlays: Severity
    p1_U <-
      ggplot(mod_subsampled, aes(x = UMAP1, y = UMAP2)) +
      geom_density2d(color="#CCCCCC") +
      geom_point_rast(size = 3, stroke = 0, aes(color=run_type_id)) + # original size = 0.5
      coord_fixed(ratio = 1) + 
      
      xlim(min(mod_subsampled$UMAP1)-1, max(mod_subsampled$UMAP1)+1) +
      ylim(min(mod_subsampled$UMAP2)-1.5, max(mod_subsampled$UMAP2)+1.5) +
      scale_color_manual(values=color4) + 
      labs(x="UMAP1", y="UMAP2", title="De Novo Spatial Regions", color="") +
      theme_bw() + 
      theme3 +
      theme(legend.key.size = unit(2, 'cm'))
    #+ 
    #guides(color = guide_legend(override.aes = list(size = 5))) #exclude this line if using empty template
    p1_U
    ggsave(filename = 'UMAP_severity_blank_pt3.pdf', plot = p1_U, width = 20, height = 20, dpi = 300, device = 'pdf',
           path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
    )
  #### umap marker overlays (tangles) ####
  path <- '/Volumes/BryJC_Stanford/paper1_final/Fig4/v2/UMAP/endothelial'
  dir.create(path, showWarnings = FALSE)
  for (marker in markers) {
    d1 <- filter(obj_umap_plot_all_data, UMAP1 > -20 & UMAP2 < 10)
    d2 <- filter(obj_umap_plot_all_data, UMAP1 > -20 & UMAP2 < 10 & de_novo_regions %in% c(1,2,3,4,5) & FlowSOM_ids %in% c('endothelial')) #FlowSOM_ids %in% c('neurons','tangles'))
    p1_U <- 
      ggplot(d1, aes(x = UMAP1, y = UMAP2)) +
      geom_density2d(color="#CCCCCC") +
      geom_point(data=d1, size = 0.5, color='ivory3') + 
      geom_point(data=d2, size = 1.5, stroke = 0, aes(color=eval(parse(text=marker)))) + 
      coord_fixed(ratio = 1) +
      theme_bw() + labs(x="UMAP1", y="UMAP2", title=paste0('endothelial', ':', marker), color="") +
      xlim(min(mod_subsampled$UMAP1)-1, max(mod_subsampled$UMAP1)+1) +
      ylim(min(mod_subsampled$UMAP2)-1.5, max(mod_subsampled$UMAP2)+1.5) +
      scale_color_gradient(low = colorExp[1], high = colorExp[100]) +
      theme1a
    p1_U
    picname <- paste0(path,'/',marker)
    ggsave(p1_U, file=paste(picname,"png",sep="."), width = 7, height = 7, dpi = 320)
  }
  
  #### umap marker overlays (neurons) ####
#### H) Carpet Markers based on voronoi masks ####
  
  #create vector for mean normalization and subtract from each cell
    single_mean_norm_vector <- apply(pixel_regions_info[setdiff(names(pixel_regions_info), c('Severity','De novo regions'))], 2, mean)
    pixel_regions_info[,setdiff(names(pixel_regions_info), c('Severity','De novo regions'))] <- pixel_regions_info[setdiff(names(pixel_regions_info), c('Severity','De novo regions'))] - single_mean_norm_vector
  #create dataframe of mean values for groups
  pixel_regions_info %>%
    dplyr::group_by(`De novo regions`, Severity) %>%
    dplyr::summarize_if(is.numeric, funs(mean)) %>%
    ungroup() ->
    heat_pixels_regions
  # pheatmap(heat_pixels_regions_mat[, carpet_markers], cluster_rows = T, cluster_cols = T, border_color = FALSE, scale = "column", main = "Mean expression per de novo region", 
  #          annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row =T, gaps_row = c(7,8),
  #          col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 25, height = 10,
  #          filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v3/carpet_hm.pdf')
  # Convert data.frame to matrix for making heatmap with pheatmap
  heat_pixels_regions_mat <- data.matrix(heat_pixels_regions)
  
  # Make the heatmap
  carpet_markers = c("PHF1Tau", "Amyloidbeta140", "Amyloidbeta142", "PanAmyloidbeta1724", "CD47", "PanGAD6567", "Parvalbumin", 
                     "PSD95", "Synaptophysin", "VGAT", "VGLUT1", "VGLUT2", "Calretinin", "Calbindin", "MAP2", "TotalTau", "CD56Lyo", 
                     "MBP", "MAG", "GFAP")
  nu_panel =  setdiff(markers, c('HistoneH3Lyo','EEA1','SERT','TH','Reelin','X8OHGuano'))
  colnames(heat_pixels_regions_mat) = names(heat_pixels_regions)
  rownames(heat_pixels_regions_mat) = 1:nrow(heat_pixels_regions_mat)
  
  heat_pixels_regions_mat <- as.data.frame(heat_pixels_regions_mat)
  anno_data = heat_pixels_regions[, c('De novo regions', 'Severity')]
  anno_data = as.data.frame(anno_data)
  
  #for heatmap colors
  mycolors <- list(
    `De novo regions` = c(ND = "blue1", GD = "yellow3", PD = "green1", TD = "brown", MD = "orchid1"),
    Severity = c(Ctrl = "dodgerblue2", MedAD = 'green4', HiAD = '#E31A1C'))
  lgl <- c('Ctrl','MedAD','HiAD','ND','GD','PD','TD','MD')
  pheatmap(heat_pixels_regions_mat[, carpet_markers], cluster_rows = T, cluster_cols = T, border_color = FALSE, scale = 'column', main = "Mean expression per de novo region", 
           legend_labels = lgl, annotation_row = anno_data, annotation_colors = mycolors, show_rownames = T, annotation_names_row =T,
           col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 25, width = 25, height = 15,
           filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v3/hm_carpet_scaled.pdf')
#### I) Heatmap: show distribution of marker differences in different cell types between de novo regions across severity ####
  
  #create vector for mean normalization and subtract from each cell
  single_mean_norm_vector <- apply(master_regions_c[,panel], 2, mean)
  master_regions_c2 <- master_regions_c
  master_regions_c2[,panel] <- master_regions_c2[,panel] - single_mean_norm_vector
  levels(master_regions_c2$de_novo_regions) <- dnr_names # factorize region
  master_regions_c2 %>% rename(Severity = run_type_id, `De novo regions` = de_novo_regions) -> master_regions_c2
  #create dataframe of mean values for groups
  master_regions_c2 %>%
    #filter(anat_id == c()) %>%
    dplyr::group_by(`De novo regions`, FlowSOM_ids) %>%
    dplyr::summarize_if(is.numeric, funs(mean)) %>%
    ungroup() ->
    heat_denovo_regions
  #for heatmap colors
  mycolors <- list(
    `De novo regions` = c(ND = "blue1", GD = "yellow3", MD = "orchid1", APD = "green1", TTD = "brown"),
    Severity = c(Ctrl = "dodgerblue2", MedAD = 'green4', HiAD = '#E31A1C'))
  # original colors for non-z scored expr
  #col = colorRampPalette(c("#202020","#00CCCC"))(100)
  
  #neurons
  heat_denovo_regions_2 <- filter(heat_denovo_regions, FlowSOM_ids == 'neurons')
  #organize naming, ordering schema
  hgr_panel = c('CD47','TotalTau','Calretinin','CD56Lyo','Parvalbumin','MAP2','PolyubiK48','MFN2','PanGAD6567','VGAT','PolyubiK63','Synaptophysin','VGLUT1','VGLUT2','Calbindin','PSD95')
    hgr_panel_mini = c('PolyubiK48','MFN2','PolyubiK63')
    #heat_denovo_regions_2$anat_id = droplevels(heat_denovo_regions_2$anat_id)
  heat_denovo_regions_2 <- as.data.frame(heat_denovo_regions_2)
  rownames(heat_denovo_regions_2) <- c(1:nrow(heat_denovo_regions_2))
  anno_data = heat_denovo_regions_2[c('De novo regions')]
  #make heatmap
  I1 <- pheatmap(heat_denovo_regions_2[,hgr_panel_mini], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Neurons", 
                 annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, annotation_legend = FALSE,
                 col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 25, width = 15, height = 15,
                 filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4/hm_neurons_mini_denovo_meannorm_inferno_noSeverity.pdf')
  I1
  dev.off()
  
  #non-immune glia
  heat_denovo_regions_2 <- filter(heat_denovo_regions, FlowSOM_ids == 'non-immune glia')
  #organize naming, ordering schema
  hgr_panel = c('GFAP','MBP','MAG','ApoE4')
    #heat_denovo_regions_2$anat_id = droplevels(heat_denovo_regions_2$anat_id)
  heat_denovo_regions_2 <- as.data.frame(heat_denovo_regions_2)
  rownames(heat_denovo_regions_2) <- c(1:nrow(heat_denovo_regions_2))
  anno_data = heat_denovo_regions_2[c('De novo regions')]
  #make heatmap
  I2 <- pheatmap(heat_denovo_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Non-immune Glia", 
                 annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, annotation_legend = FALSE,
                 col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 25, width = 15, height = 15, 
                 filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4/hm_nig_denovo_meannorm_inferno_noSeverity.pdf')
  I2
  dev.off()
  
  #endothelial
  heat_denovo_regions_2 <- filter(heat_denovo_regions, FlowSOM_ids == 'endothelial')
  #organize naming, ordering schema
  hgr_panel = c('CD31','CD105','MCT1')
    #heat_denovo_regions_2$anat_id = droplevels(heat_denovo_regions_2$anat_id)
  heat_denovo_regions_2 <- as.data.frame(heat_denovo_regions_2)
  rownames(heat_denovo_regions_2) <- c(1:nrow(heat_denovo_regions_2))
  anno_data = heat_denovo_regions_2[c('De novo regions')]
  #make heatmap
  I3 <- pheatmap(heat_denovo_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Endothelial", 
                 annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F, annotation_legend = FALSE,
                 col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 25, width = 15, height = 15,
                 filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4/hm_endo_denovo_meannorm_inferno_noSeverity.pdf')
  I3
  dev.off()
  
  #microglia
  heat_denovo_regions_2 <- filter(heat_denovo_regions, FlowSOM_ids == 'microglia')
  #organize naming, ordering, coloring schema  
  hgr_panel = setdiff(c(glia_struct3), c('X8OHGuano','pTDP43','EEA1','CD47'))
    #heat_denovo_regions_2$anat_id = droplevels(heat_denovo_regions_2$anat_id)
  heat_denovo_regions_2 <- as.data.frame(heat_denovo_regions_2)
  rownames(heat_denovo_regions_2) <- c(1:nrow(heat_denovo_regions_2))
  anno_data = heat_denovo_regions_2[c('De novo regions')]
  #make heatmap  
  I4 <- pheatmap(heat_denovo_regions_2[,hgr_panel], cluster_rows = F, cluster_cols = T, border_color = FALSE, main = "Microglia", 
                 annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = F,
                 col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 25, width = 15, height = 15,
                 filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4/hm_mgla_denovo_meannorm_inferno_noSeverity.pdf')
  I4
  dev.off()
  
  
#### J) Barplot: synapse percent distributions ####
  # show increase in positive tau, excitiory later, inhibitory first?
  J1 <- ggplot(filter(syn_percent, `Disease type` == 'PHF1Tau positive' & `De novo region` %in% c("ND","GD","APD","TTD","MD"))) +
    aes(x = `Synapse type`, fill = `Synapse type`, weight = `Percent positive`) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = c('blue','purple')) +
    theme_bw() +
    theme1a +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(vjust=0.5, angle=35)) + #gets rid of lines
    theme(legend.key.size = unit(2, 'cm')) +
    labs(x="Synapse type", y="Synaptic pixels - PHF1Tau positivity (%)", title="Proteopathy colocalizes with synaptic signal", fill="Synapse type") +
    theme(panel.spacing = unit(2, "lines")) +
    facet_grid(cols = vars(`De novo region`))
  
  J1
  ggsave(filename = 'barplot_ratio_synpixels_denovoregions_tauONLY_revised.pdf', plot = J1, width = 20, height = 15, dpi = 300, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4'
  )
#### FlowSOM Clustering Neurons ONLY ####
  #filter for just neurons
  fsom_clustering_data <- filter(disease_cell_overlap, FlowSOM_ids == "neurons")
  #convert all factors to numeric                                  
  fsom_clustering_data <- mutate_if(fsom_clustering_data, is.factor, ~ as.numeric(.x))
  
  flowFrame_cluster <- flowCore::flowFrame(as.matrix(fsom_clustering_data))

    # initial cluster number  
    num_clusters <- 100
    # color palette of choice, here using c25 from brain_data_pkg
    palette <- c25
    
    #clustering_markers <- c('TH', 'Calretinin', 'Calbindin', 'Parvalbumin', 'PanGAD6567', 'VGAT', 'VGLUT2', 'MAP2', 'MBP', 'MAG', 'GFAP', 'Iba1', 'CD45', 'MCT1', 'CD31', 'CD105')
    clustering_markers <- c("CD47",               "Calretinin",         "Parvalbumin",        "PolyubiK48",         "MFN2",              
                            "PanGAD6567",         "VGAT",               "PolyubiK63",         "Synaptophysin",      "VGLUT1",             "VGLUT2",             "Calbindin",          "PSD95",             
                            "Presenilin1NTF",     "PHF1Tau",            "Amyloidbeta142",     "Amyloidbeta140",     "PanAmyloidbeta1724")
    exclude <- c("C12", "Na23", "Si28", "Ca40", "Background", "Ta181", "Au197", "empty113" )
    markers <- setdiff(panel, exclude)
  
  # run FlowSOM on the above flowFrame which has already had transformation performed 
  # to adjust initial cluster numbers: xdim=10, ydim=10, maxMeta = 50
  fSOM2 <- FlowSOM(flowFrame_cluster,
                   compensate = F, scale = F, colsToUse = clustering_markers, nClus = 30, seed = 123)
  
  table(fSOM2$metaclustering) 
  metaClustering_consensus <- metaClustering_consensus(fSOM2$FlowSOM$map$codes,k=5)
  metaClustering2 <- fSOM2$metaclustering
  
  # visualize manual gating over metaclusters
  PlotStars(fSOM2[[1]], backgroundValues = as.factor(metaClustering_consensus))
  PlotPies(fSOM2$FlowSOM, cellTypes=cluster_comparisons$og_cid, backgroundValues = as.factor(metaClustering2))
  PlotNumbers(fSOM2$FlowSOM, backgroundValues = as.factor(metaClustering2))
  PlotMarker(fSOM2[[1]], "TotalTau")
  
  # isolate cluster and metaclustering data
  fSOM2_clustering <- data.frame(fSOM2$FlowSOM$map$mapping)
  colnames(fSOM2_clustering) <- c("Cluster", "Value")
  fSOM2_clustering$Meta <- as.numeric(fSOM2$metaclustering[fSOM2_clustering[,1]])
  # attach above data to manual gated data to compare
  fsom_clustering_data$FlowSOM_ids <- fSOM2_clustering$Cluster
  
  # calculate heatmap from all clusters - by Meta
  fsom_clustering_data %>%
    dplyr::group_by(Meta) %>%
    dplyr::summarize_if(is.numeric, funs(median)) %>%
    ungroup() ->
    heat_f
  rownames(heat_f) <- heat_f$Meta
  heatmaply(heat_f[,clustering_markers], labRow = rownames(heat_f), main="Meta")
  
  # calculate heatmap from all clusters - by Meta
  fsom_clustering_data %>%
    filter(Meta %in% c(10)) %>%
    dplyr::group_by(FlowSOM_ids) %>%
    dplyr::summarize_if(is.numeric, funs(median)) %>%
    ungroup() ->
    heat_f
  rownames(heat_f) <- heat_f$FlowSOM_ids
  heatmaply(heat_f[,clustering_markers], labRow = rownames(heat_f), main="Meta = 10")
  
  
  ### reassign metaclusters based upon plot inspection ###
  fsom_clustering_data[, "Meta"] = 0
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(1,12,23,42,32,21,11))) 
  fsom_clustering_data[data_indexes, "Meta"] = 1
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(22,51,2,24)))
  fsom_clustering_data[data_indexes, "Meta"] = 2
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(49,40,59,50,60)))
  fsom_clustering_data[data_indexes, "Meta"] = 3 
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(39,30,20,29,28,18,34,8,9,10,19)))
  fsom_clustering_data[data_indexes, "Meta"] = 4 
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(91,61,71,81)))
  fsom_clustering_data[data_indexes, "Meta"] = 5 
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(14,46,15,25,35,45,36,26,55,65,66,56,72,4,43,24)))
  fsom_clustering_data[data_indexes, "Meta"] = 6 
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(27,38,37,16,17,82,69,67,58,48,47,92,93,83,84,73,74,94)))
  fsom_clustering_data[data_indexes, "Meta"] = 7 
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(75,41,31)))
  fsom_clustering_data[data_indexes, "Meta"] = 8 
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(70,80,90,88,77,87,76,85,95,96,97,78,79,99,100,98,86,89,57,68)))
  fsom_clustering_data[data_indexes, "Meta"] = 9
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(13,5,7,33,3)))
  fsom_clustering_data[data_indexes, "Meta"] = 10
  # first-pass cluster
  data_indexes = which(with(fsom_clustering_data, FlowSOM_ids %in% c(53,64,54,63,44,6,43,62,52)))
  fsom_clustering_data[data_indexes, "Meta"] = 11
  
  
  ### create new data.frame 'fully_FlowSOM' with combined newly labeled cells and ez_objects ###
  fully_FlowSOM_neurons <- fsom_clustering_data

#### UMAP neurons only ####  
  # sample (sample_n) for later use by object type
  n_sub_fraction <- 0.15
  
  subsetted_norm_data <- data.frame()
  for (id in c('tangles','plaques','microglia','endothelial','neurons','non-immune glia')) {
    typed_obj <- filter(master_regions_c, FlowSOM_ids==id)
    sampled_obj <- dplyr::sample_frac(typed_obj, n_sub_fraction)
    subsetted_norm_data <- rbind(subsetted_norm_data, sampled_obj)
  }
  
  #umap testing
  set.seed(123)
  umap_markers <- setdiff(markers, c('HistoneH3Lyo','CD56Lyo','EEA1','SERT','TH','Reelin'))
  #c('Iba1', 'CD45', 'GFAP', 'CD31', 'CD105', 'MCT1', 'MAP2')
  #setdiff(markers, c('HistoneH3Lyo','CD56Lyo',))
  obj_umap <- fully_FlowSOM_neurons[clustering_markers]
  obj_umap <- as.matrix(obj_umap)
  
  custom_config <- umap.defaults # Run umap.defaults to see original settings
  custom_config$n_neighbors <- 25 # knn
  #custom_config$min_dist <- 0.7 # distance betweeen 2 closet cells.  default 0.1, at 0.05 farway, at 0.5 close together
  custom_config$local_connectivity <- 1 #
  custom_config$negative_sample_rate <- 25 # choose 5 random points
  custom_config$knn_repeats <- 3 # 3 repeats for each iteration or epoch of 200 attempts
  custom_config$a = 5 # higher numbers brings similar clusters closer
  custom_config$b = 3 # higher numbers brings different clusters closer
  #custom_config$metric <- 'manhattan'
  
  #make the umap
  out_obj_umap <- umap(obj_umap, custom_config)
  
  obj_umap_plot <- as.data.frame(out_obj_umap$layout)
  colnames(obj_umap_plot) <- c("UMAP1", "UMAP2")
  obj_umap_plot_all_data <- obj_umap_plot
  obj_umap_plot_all_data <- cbind(obj_umap_plot, fully_FlowSOM_neurons)
  obj_umap_plot_all_data$Meta <- as.factor(obj_umap_plot_all_data$Meta)
  obj_umap_plot_all_data$run_type_id <- as.factor(obj_umap_plot_all_data$run_type_id)
  obj_umap_plot_all_data$de_novo_regions <- as.factor(obj_umap_plot_all_data$de_novo_regions)
  #
  #### subsetting and markers ####
    gathered_umap_data <- gather(obj_umap_plot_all_data, key = "markers", value = "expr", -c(setdiff(names(obj_umap_plot_all_data), markers)))
    #gathered_umap_data[which(gathered_umap_data$expr > 0.5), 'expr'] <- 0.5
    
    # sample (sample_n) for later use by ggplot
    n_sub_fraction <- 0.2
    
    subsetted_plot_data <- data.frame()
    for (id in c('tangles','plaques','microglia','endothelial','neurons','non-immune glia')) {
      typed_obj <- filter(obj_umap_plot_all_data, FlowSOM_ids==id)
      sampled_obj <- dplyr::sample_frac(typed_obj, n_sub_fraction)
      subsetted_plot_data <- rbind(subsetted_plot_data, sampled_obj)
    }
    gathered_umap_data_sample <- gather(subsetted_plot_data, key = "markers", value = "expr", -c(setdiff(names(obj_umap_plot_all_data), markers)))
    #gathered_umap_data_sample[which(gathered_umap_data$expr > 0.5), 'expr'] <- 0.5
    
  #### PLOT REPRESENTATION ####
    p1_U <- filter(obj_umap_plot_all_data) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
             ggplot(aes(x = UMAP1, y = UMAP2, color=run_type_id)) +
             geom_point(size = 0.5) + 
             scale_color_manual(values = color4) +
             coord_fixed(ratio = 1) +
             theme_bw(20) +
             guides(colour = guide_legend(override.aes = list(size=4)))
            #facet_wrap(vars(markers))
    p1_U
    
#### Test neuronal subpops ####
  # adjust dataset
  fsom_neurons <- filter(disease_cell_overlap, FlowSOM_ids == "neurons")
  fsom_neurons <- cbind(obj_umap_plot, fsom_neurons)
  fsom_neurons <- inner_join(fsom_neurons, fsom_clustering_data[,c('spatial_id','Meta')], by='spatial_id')
  fsom_neurons <- left_join(fsom_neurons, master_gated_anat2[,c("spatial_id","anat_id")])
  fsom_neurons$Meta <- as.factor(fsom_neurons$Meta)
  fsom_neurons %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Disease clusters", "Non-disease clusters")) -> fsom_neurons
  fsom_neurons$MetaDisease <- as.factor(fsom_neurons$MetaDisease)
  # viz umap
  d1 <- filter(obj_umap_plot_all_data, UMAP1 > -20 & UMAP2 > -20)
  d2 <- filter(obj_umap_plot_all_data, UMAP1 > -20 & UMAP2 > -20 & de_novo_regions %in% c(1,2,3,5,4) & FlowSOM_ids %in% c('plaques')) #FlowSOM_ids %in% c('neurons','tangles'))
  p1_U <- 
    ggplot(d1, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data=d1, size = 0.5, color='ivory3') + 
    geom_point(data=d2, size=2, aes(color=Amyloidbeta142)) + 
    scale_color_gradient(low='blue', high='red', limits=c(0,1)) +
    coord_fixed(ratio = 1) +
    facet_wrap(vars(de_novo_regions))
  p1_U
  # Additional UMAPS #
  p1_U <- filter(fsom_neurons) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color=run_type_id)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = color4) +
    coord_fixed(ratio = 1) +
    theme_bw(20) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  #facet_wrap(vars(markers))
  p1_U
  
  
#### Proportions of neuron subpops in severity/regions ####
  #color 6 for cluster colors, color 7 for patho disease/n disease
  #collect neuron counts and normalize by area
  fsom_neurons %>%
     dplyr::group_by(run_type_id, Meta) %>%
     dplyr::summarise(n = n()) %>%
     rowwise() %>%
     mutate(meta_area_norm = n/run_areas[run_type_id]) -> fsom_neurons_meta_norm
  fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10,11), "Disease clusters", "Non-disease clusters")) -> fsom_neurons_meta_norm
  #Severity
  ggplot(fsom_neurons_meta_norm) +
    aes(x = run_type_id, fill = Meta, weight = meta_area_norm) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = color6) +
    theme_bw(20)
  ggplot(fsom_neurons_meta_norm) +
         aes(x = run_type_id, fill = Meta, weight = meta_area_norm) +
         geom_bar(position = "dodge") +
         scale_fill_manual(values = color6) +
         theme_bw(20) +
         facet_wrap(vars(MetaDisease))
  
  #collect neuron counts and normalize by num of neurons in de novo region
  fsom_neurons_total_denovo <- t(fsom_neurons %>% dplyr::group_by(de_novo_regions) %>% summarise(n = n()) %>% select(n))
  fsom_neurons %>%
    dplyr::group_by(de_novo_regions, Meta) %>%
    dplyr::summarise(n = n()) %>%
    rowwise() %>%
    mutate(meta_n_norm = n/fsom_neurons_total_denovo[de_novo_regions]) -> fsom_neurons_meta_norm
  fsom_neurons_meta_norm %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10,11), "Disease clusters", "Non-disease clusters")) -> fsom_neurons_meta_norm
  #De novo regions
  filter(fsom_neurons_meta_norm, de_novo_regions %in% c(1,4,5)) %>%
  ggplot() +
    aes(x = factor(de_novo_regions, level = c(1,5,4)), fill = Meta, weight = meta_n_norm) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = color6) +
    theme_bw(20)
  ggplot(fsom_neurons_meta_norm) +
         aes(x = MetaDisease, fill = Meta, weight = meta_n_norm) +
         geom_bar(position = "dodge") +
         scale_fill_manual(values = color6) +
         theme_bw(20) +
         facet_grid(cols = vars(factor(de_novo_regions, level = c(1,2,5,3,4))))
  
#### Anat subsets ####
  fsom_neurons <- left_join(fsom_neurons, master_gated_anat2[,c("spatial_id","anat_id")])
  
  
  
#### Marker distros ####
  #heatmap
  fsom_neurons %>%
    dplyr::group_by(Meta) %>%
    dplyr::summarize_if(is.numeric, funs(median)) %>%
    ungroup() ->
    heat_f
  
  heat_f <- as.data.frame(heat_f)
  rownames(heat_f) <- heat_f$Meta
  anno_data = heat_f[,c('Meta')]
  anno_data = as.data.frame(anno_data)
  names(anno_data) = 'Meta'
  mycolors <- list(
    Meta = c(`1`='dodgerblue2',`2`='gold1',`3`='skyblue2',`4`="#FB9A99",`5`="palegreen2",`6`="#CAB2D6",`7`="#FDBF6F",`8`="darkorange4",`9`="khaki2",`10`="maroon",`11`="orchid1"))
  
  h1u <- pheatmap(heat_f[,clustering_markers], cluster_rows = T, cluster_cols = T, border_color = FALSE, main = "Neurons", 
                 annotation_row = anno_data, annotation_colors = mycolors, show_rownames = T, annotation_names_row = F, col = colorRampPalette(c("#202020","#00CCCC"))(100),
                 cellwidth = 30, cellheight = 30, fontsize = 10, width = 10, height = 10,
                 filename = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v3/meta_neurons_hm.pdf')
  h1u
  
  #boxplot --> MFN2+
  filter(fsom_neurons, de_novo_regions %in% c(1,2,3,4,5) & Meta %in% c(1,2,3,4,10)) %>%
  ggplot(fsom_neurons) +
    aes(x = MetaDisease, y = MFN2, fill = Meta) +
    geom_boxplot() +
    scale_fill_hue() +
    theme_minimal() +
    facet_grid(cols = vars(factor(de_novo_regions, level = c(1,2,5,3,4))))
  #line plot - MFN2 (MFN2(+) neurons vs MFN2(-) or respective clusters --> look at progression in severity, path disease)
      # fsom_neurons <- mutate(fsom_neurons, MFN2_Pos = ifelse(MFN2 > 0, 1, 0))
      # fsom_neurons$MFN2_Pos <- as.factor(fsom_neurons$MFN2_Pos)
      # fsom_neurons <- mutate(fsom_neurons, MFN2_Meta_Pos = ifelse(Meta %in% c(1,2,4,10,3), 1, 0))
      # fsom_neurons$MFN2_Meta_Pos <- as.factor(fsom_neurons$MFN2_Meta_Pos)
    fsom_neurons %>%
      dplyr::group_by(run_type_id, de_novo_regions, Meta) %>% 
      dplyr::summarize_if(is.numeric, funs(median)) %>%
      ungroup() -> fsom_neurons_markers
    
    filter(fsom_neurons, de_novo_regions %in% c(1,2,3,4,5) )%>%#& Meta %in% c(1,2,3,4,10)) %>%
    ggplot() +
      aes(x = factor(de_novo_regions, level = c(1,2,3,5,4)), y = MFN2, group = Meta, color = Meta) +
      geom_point(stat='summary', fun=mean) +
      stat_summary(fun=mean, geom="line") +
      scale_color_manual(values = color7) +
      theme_minimal() +
      facet_wrap(vars())
    
    filter(fsom_neurons, de_novo_regions %in% c(1,2,3,4,5) )%>%#& Meta %in% c(1,2,3,4,10)) %>%
      ggplot() +
      aes(x = run_type_id, y = MFN2, group = Meta, color = Meta) +
      geom_point(stat='summary', fun=mean) +
      stat_summary(fun=mean, geom="line") +
      scale_color_manual(values = color6) +
      theme_minimal() +
      facet_wrap(vars())
    
  #stat tests
    test <- compare_means(MFN2~run_type_id, data=fsom_neurons, group.by="Meta")
    ggboxplot(ToothGrowth, x = "dose", y = "len",
              color = "dose", palette = "jco")+
      stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
      stat_compare_means(label = "p.signif", method = "t.test",
                         ref.group = ".all.")
    
    

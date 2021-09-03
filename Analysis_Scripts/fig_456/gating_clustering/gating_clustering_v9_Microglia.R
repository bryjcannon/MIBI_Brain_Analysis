microglia_umap_data <- filter(disease_cell_overlap, FlowSOM_ids == "microglia")
microglia_umap_data <- left_join(microglia_umap_data, master_gated_anat2[,c("spatial_id","anat_id")])
microglia_umap_data <- filter(microglia_umap_data, anat_id %in% c("DG","CA4","CA3","CA2","CA1"))

n_sub_fraction <- 554 #total Ctrl microglia amount

subsetted_norm_data <- data.frame()
for (id in c("Ctrl","MedAD","HiAD")) {
  typed_obj <- filter(microglia_umap_data2, run_type_id==id)
  sampled_obj <- dplyr::sample_n(typed_obj, n_sub_fraction)
  subsetted_norm_data <- rbind(subsetted_norm_data, sampled_obj)
}
subsetted_norm_data$UMAP1 <-NULL
subsetted_norm_data$UMAP2 <-NULL

#umap testing
set.seed(123)
umap_markers <- c("Iba1","CD45","ApoE4","CD33Lyo")
#c('Iba1', 'CD45', 'GFAP', 'CD31', 'CD105', 'MCT1', 'MAP2')
#setdiff(markers, c('HistoneH3Lyo','CD56Lyo',))
obj_umap <- subsetted_norm_data[umap_markers]
obj_umap <- as.matrix(obj_umap)

custom_config <- umap.defaults # Run umap.defaults to see original settings
custom_config$n_neighbors <- 5 # knn
#custom_config$min_dist <- 0.7 # distance betweeen 2 closet cells.  default 0.1, at 0.05 farway, at 0.5 close together
custom_config$local_connectivity <- 1 #
custom_config$negative_sample_rate <- 25 # choose 5 random points
custom_config$knn_repeats <- 1 # 3 repeats for each iteration or epoch of 200 attempts
custom_config$a = 7 # higher numbers brings similar clusters closer 
custom_config$b = 5 # higher numbers brings different clusters closer 
#custom_config$metric <- 'manhattan'

#make the umap
out_obj_umap <- umap(obj_umap, custom_config)

obj_umap_plot <- as.data.frame(out_obj_umap$layout)
colnames(obj_umap_plot) <- c("UMAP1", "UMAP2")
obj_umap_plot_all_data <- obj_umap_plot
obj_umap_plot_all_data <- cbind(obj_umap_plot, subsetted_norm_data)
#obj_umap_plot_all_data$Meta <- as.factor(obj_umap_plot_all_data$Meta)
obj_umap_plot_all_data$run_type_id <- as.factor(obj_umap_plot_all_data$run_type_id)
obj_umap_plot_all_data$de_novo_regions <- as.factor(obj_umap_plot_all_data$de_novo_regions)
#

#### additional alterations ####
  n_sub_fraction <- 554 #total Ctrl microglia amount
  
  subsetted_norm_data <- data.frame()
  for (id in c("Ctrl","MedAD","HiAD")) {
    typed_obj <- filter(microglia_umap_data2, run_type_id==id)
    sampled_obj <- dplyr::sample_n(typed_obj, n_sub_fraction)
    subsetted_norm_data <- rbind(subsetted_norm_data, sampled_obj)
  }
  
  p1_U <- filter(microglia_umap_data) %>% #, UMAP1 > -20 & UMAP2 > -20) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color=run_type_id)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = color4) +
    coord_fixed(ratio = 1) +
    theme_bw(20) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  #facet_wrap(vars(markers))
  p1_U


#### adjust dataset ####
#microglia_umap_data <- filter(disease_cell_overlap, FlowSOM_ids == "microglia")
microglia_umap_data <- cbind(obj_umap_plot, microglia_umap_data)
subsetted_norm_data <- cbind(obj_umap_plot, subsetted_norm_data)
#microglia_umap_data <- inner_join(microglia_umap_data, fsom_clustering_data[,c('spatial_id','Meta')], by='spatial_id')
#microglia_umap_data <- left_join(microglia_umap_data, master_gated_anat2[,c("spatial_id","anat_id")])
#microglia_umap_data$Meta <- as.factor(microglia_umap_data$Meta)
#microglia_umap_data %>% mutate(MetaDisease = ifelse(Meta %in% c(1,2,6,8,10), "Disease clusters", "Non-disease clusters")) -> fsom_neurons
#microglia_umap_data$MetaDisease <- as.factor(microglia_umap_data$MetaDisease)
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


# mean normalize
single_mean_norm_vector <- apply(microglia_umap_data[,c("Iba1","CD45","ApoE4","CD33Lyo","PHF1Tau")], 2, mean)
microglia_umap_data2 <- microglia_umap_data
microglia_umap_data2[,c("Iba1","CD45","ApoE4","CD33Lyo","PHF1Tau")] <- microglia_umap_data2[,c("Iba1","CD45","ApoE4","CD33Lyo","PHF1Tau")] - single_mean_norm_vector
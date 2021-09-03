
#### analysis ####
# combine all three datasets - Ctrl, MedAD, HiAD
master_obj_data_Ctrl_together <- master_obj_data_Ctrl
levels(master_obj_data_Ctrl_together$run_type_id) <- c('Healthy Control', 'Healthy Control', 'Healthy Control')
master_obj_data <- rbind(master_obj_data_Ctrl_together, master_obj_data_Med, master_obj_data_Hi)

# rename object_types (NOTE: move and use in analysis, not here)
nu_obj_types <- c("cells", "amyloid plaques", "tau tangles", "microglia", "vessels")
levels(master_data$obj_type_id) <- nu_obj_types

# assign and standardize panels (if needed if you havent done this already in transformation)
panel_start = 8
panel_end = 48
# use names(data) to figure out the indices of which columns to ignore (metals, empty, background, etc), keeping only the markers in panel
panel <- names(master_data[ ,panel_start:panel_end][, -c(4,5)]) # -c(#,#,#): remove metals, composites, other labels after already paring down columns (i.e. numbers will change)

# rename panel with shorter names (may need to reload library -> data.table)
nu_panel <- c('H3', 'CD47', 'Tau', 'ApoE4', 'ABeta40', 'ABeta42', 'Calretinin', 'CD56', 'Parvalbumin', 'PanABeta', 'CD31', 'MAP2', 'PolyubiK48', 'MFN2', 'PanGAD', 'VGAT', 'SERT', 'PolyubiK63', 'Reelin', 'EEA1', 'pTDP43', 'Synaptophysin', 'CD45', 'PanApoE', 'X8OHGuano', 'TH', 'Presenilin1NTF', 'MAG', 'VGLUT1', 'VGLUT2', 'Calbindin', 'Iba1', 'PSD95', 'MBP', 'MCT1', 'CD33', 'GFAP', 'CD105', 'PHF1Tau')
setnames(master_data, old = panel, new = nu_panel)

# SAMPLING #
# sample (sample_n) for later use by object type
n_sub_fraction <- 1500

object_types <- c(object_types, 'cell')
subsetted_norm_data <- data.frame()
for (obj_type in object_types) {
  typed_obj <- filter(data_normalized, obj_type_id == obj_type)
  sampled_data <- sample_n(typed_obj, n_sub_fraction)
  subsetted_norm_data <- rbind(subsetted_norm_data, sampled_data)
}

#### heatmaps ####
# heatmap of median signals (after processing, transformation, and normalization)

hmp_data <- as.data.table(subset(master_obj_data, select = panel_obj))
hmp_data_m <- as.matrix(subset(hmp_data, select = panel), rownames = hmp_data$obj_type_id)
hmp_matrix <- hmp_data[, lapply(.SD, median(na.rm = T)), .SDcols=panel, by=obj_type_id] %>% # get medians of channels by gates
              as.matrix(rownames="obj_type_id") # turn it into a matrix
heatmaply(hmp_matrix)

# heatmap of median signals (zero removed, either normalized or percentized in hmp output)

  # plus1_data <- master_obj_data
  # plus1_data[,panel_start:panel_end] <- plus1_data[,panel_start:panel_end] + 1
hmp_data <- panel_special_median(master_obj_data, rm_zeros = T, nu_obj_types, panel)
hmp_data_matrix <- as.matrix(hmp_data)
heatmaply(normalize(hmp_data_matrix))
heatmaply(percentize(hmp_data_matrix))

# violin plot of object size distributions (log10 scaled)
ggplot(master_obj_data) +
  aes(x = run_type_id, y = obj_size, fill = obj_type_id) +
  geom_violin(adjust = 1L, scale = "area") +
  scale_fill_hue() +
  scale_y_continuous(trans = "log10") +
  labs(x = "Region", y = "Size (pixels / object)", title = "Size distribution of objects across brain regions", fill = "Object Type") +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20))
#### ####

# expr data distributions for object types, by region
panel_obj_region <- c(nu_panel, 'obj_type_id', 'run_type_id', 'Subregions')
gathered_data <- gather(subset(master_obj_data, select = panel_obj_region), key = "markers", value = "expr", -c(run_type_id, obj_type_id, Subregions))
gathered_data$expr <- gathered_data$expr+1 # for log+1 transformations
gathered_data_disease_struct <- filter(gathered_data, markers == c('ABeta40','ABeta42','ApoE4', 'PanABeta', 'PanApoE', 'PHF1Tau', 'X8OHGuano', 'pTDP43'))

ggplot(gathered_data) +
  aes(y = "expr", fill = obj_type_id) + 
  geom_violin(position = "dodge", alpha = 0.5) +
  scale_fill_hue() +
  #scale_y_continuous(trans = "log10") +
  facet_wrap("markers", scales = 'free') + theme_bw()


ggplot(gathered_data) +
  aes(x = obj_type_id, y = expr, fill = obj_type_id) +
  geom_violin(adjust = 1L, scale = "area") +
  scale_fill_hue() +
  scale_y_continuous(trans = "log10") +
  labs(x = "Objects", y = "Expression (log 10)", title = "Marker expression distribution in non-cellular objects", subtitle = "Hi AD Hippocampus") +
  facet_wrap(vars(markers), scales = 'free')


ggplot(gathered_data_filt) +
  aes(x = markers, y = expr, fill = markers) +
  geom_violin(adjust = 1L, scale = "area") +
  scale_fill_hue() +
  scale_y_continuous(trans = "log10") +
  theme_minimal() +
  facet_wrap(vars(obj_type_id))


ggplot(gathered_data) +
  aes(x = markers, fill = markers, weight = expr) +
  geom_bar() +
  scale_fill_hue() +
  theme_minimal() +
  labs(x = "markers", y = "Expression", title = "Total disease marker expression in non-cellular objects", subtitle = "Ctrl Hippocampus") +
  facet_wrap(vars(obj_type_id), scales = "free")


ggplot(gathered_data) +
  aes(x = obj_type_id, y = expr, fill = run_type_id) +
  geom_violin(adjust = 1L, scale = "area") +
  scale_fill_hue() +
  scale_y_continuous(trans = "log") +
  theme_minimal() +
  labs(x = "markers", y = "Expression", title = "Total disease marker expression in non-cellular objects", subtitle = "Across Hippocampus samples") +
  facet_wrap(vars(markers), scales = "free")

ggplot(gathered_data) +
  aes(x = run_type_id, fill = obj_type_id, weight = expr) +
  geom_bar() +
  scale_fill_hue() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
             axis.title.x = element_text(size = 5),
             axis.text.x = element_text(size = 5),
             axis.title.y = element_text(size = 5),
             axis.text.y = element_text(size = 5)) +
  facet_wrap(vars(markers))

# histone_cells_gating
master_data_all_cuberoot %>%
  filter(HistoneH3Lyo >= 0 & obj_size <= 250) %>%
  
  filter(obj_type_id %in% "cells") %>%
  ggplot() +
  aes(x = obj_size) +
  geom_histogram(bins = 30L, fill = "#0c4c8a") +
  theme_minimal() +
  facet_wrap(vars(run_type_id))

########################################################################################################################################################
########################################################################################################################################################

##### 5) VISUALISE (DIM REDUC) #####

# tSNE analysis
# prepare object data - pick a transformed data.frame and revert back to matrix for RtSNE
obj_rtsne <- obj_normalized[, panel]
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
obj_tsne_plot_all_data <- cbind(obj_tsne_plot, obj_normalized)

# plot tSNE
p1 <- ggplot(obj_tsne_plot_all_data, aes(x = tSNE1, y = tSNE2, color = obj_type_id)) +
  geom_point(size = 1) + 
  coord_fixed(ratio = 1)
p1

# UMAP analysis
obj_umap <- as.numeric(fsom_clustering_data[panel][-c(1:4,44:45)])
obj_umap <- as.matrix(obj_umap)
out_obj_umap <- umap(obj_umap)

obj_umap_plot <- as.data.frame(out_obj_umap$layout)
colnames(obj_umap_plot) <- c("UMAP1", "UMAP2")
obj_umap_plot_all_data <- obj_umap_plot
obj_umap_plot_all_data <- cbind(obj_umap_plot, fsom_clustering_data)

p1_U <- ggplot(obj_umap_plot_all_data, aes(x = UMAP1, y = UMAP2, color = CD47)) +
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

            
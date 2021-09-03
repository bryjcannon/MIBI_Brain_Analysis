# Read in RGB mask TIF(s) and plot a heatmap of mean expression values per region (and per sample)

#### Import All ####
  # Read in csv with pixel per mask, run info (generated in MATLAB)
  pixel_regions_info <- read_csv('/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/data_outputs/mask_expansion_pixel_capture/pixel_data_denovoregions_voronoi/voronoi_full/DeNovo_pixel_info.csv',  col_names = TRUE)
  pixel_regions_info$Run <- factor(pixel_regions_info$Run, levels = c(1,2,3)) # factorize run
  levels(pixel_regions_info$Run) <- c('HiAD', 'MedAD', 'Ctrl')
  pixel_regions_info$Mask <- factor(pixel_regions_info$Mask, levels = c(1,2,3,4,5)) # factorize region
  levels(pixel_regions_info$Mask) <- dnr_names # factorize region
  pixel_regions_info %>% rename(Severity = Run, `De novo regions` = Mask) -> pixel_regions_info

#### Import Tangles & Plaques ####
  # Read in csv with pixel per mask, run info (generated in MATLAB)
  pixel_regions_info_v2 <- c();
  for (i in c('plaques')) {
    temp <- read_csv(paste0('/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi_',i,'/separated_masks_voronoi/pixel_info.csv'),  col_names = TRUE)
    temp$Run <- factor(temp$Run, levels = c(1,2,3)) # factorize run
    levels(temp$Run) <- c('HiAD', 'MedAD', 'Ctrl')
    temp$Mask <- factor(temp$Mask) # factorize region
    type_factor <- rep(i, dim(temp)[1]) # create column of length object number with type info
    temp <- cbind(temp, type_factor)
    pixel_regions_info_v2 <- rbind(pixel_regions_info_v2, temp)
    rm(temp)
  }

# Create an alternative data.frame where all pixels with no + value are excluded
  pixel_regions_info_pos <- filter_all(pixel_regions_info, any_vars(. > 0))

#### Marker Summary Analysis Code ####
  # Gather data by Region (mask) and Severity (run), grab mean value, correlation, or other metrics to plot
  #notes: To fix the heatmap, all you need to do is run
    # plot_data = as.data.frame(plot_data)
    # do that before you make the annotation data frame. It seems pheatmap does not understand tibbles
  pixel_regions_info %>%
    dplyr::group_by(`De novo regions`) %>%
    dplyr::summarize_if(is.numeric, funs(mean)) %>%
    ungroup() ->
    heat_pixels_regions
  
  heat_pixels_regions_OG <- heat_pixels_regions <- heat_pixels_regions_OG # use instead of re-running heat_pixel_regions grouping
  
  # Make the heatmap
  carpet_markers = c("PHF1Tau", "Amyloidbeta140", "Amyloidbeta142", "PanAmyloidbeta1724", "CD47", "PanGAD6567", "Parvalbumin", 
                     "PSD95", "Synaptophysin", "VGAT", "VGLUT1", "VGLUT2", "Calretinin", "Calbindin", "MAP2", "TotalTau", "CD56Lyo", 
                     "MBP", "MAG", "GFAP")
  nu_panel =  setdiff(markers, c('HistoneH3Lyo','EEA1','SERT','TH','Reelin','X8OHGuano'))
  
  colnames(heat_pixels_regions) = names(heat_pixels_regions)
  rownames(heat_pixels_regions) = 1:nrow(heat_pixels_regions)
  heat_pixels_regions<- as.data.frame(heat_pixels_regions)
  anno_data = heat_pixels_regions[, c("De novo regions")]
  anno_data = as.data.frame(anno_data)
  names(anno_data) = c("De novo regions")
  
  # names for heatmap
  hm_names <- names(heat_pixels_regions)[-c(1,2)]
  
  #for heatmap colors
  mycolors <- list(
    `De novo regions` = c(ND = "#EE00FF", GD = "#00CD25", APD = "#00D1FF", TTD = "#6D2BA4", MD = "#FFCC84"),
    Severity = c(Ctrl = "#46B2AC", MedAD = '#D1EF54', HiAD = '#E22D43'))
  
  mycolors <- list(
    `De novo regions` = c(ND = "#EE00FF", GD = "#00CD25", APD = "#00D1FF", TTD = "#6D2BA4", MD = "#FFCC84"))
  
  h2u <- pheatmap(heat_pixels_regions[, nu_panel], cluster_rows = T, cluster_cols = T, border_color = FALSE, scale = "column", main = "Mean expression per de novo region", 
           annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row = T, cutree_rows = 2, cutree_cols = 2, treeheight_row = 100,
           col = inferno(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 25, height = 10,
           filename = '/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Final_Figures/FigDeNovo/v4/carpet_hm_noSeverity_colScale.pdf')
  h2u
  
  #test
  filtered_data <- filter(heat_pixels_regions)
  heatmaply(filtered_data[,c(carpet_marker)], scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_pixels_regions)[c(1:5)], 
            col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
  
  heatmaply(filtered_data[,c('Run','Mask',carpet_markers)], scale = 'column', Rowv = TRUE, Colv = TRUE, labRow = rownames(heat_pixels_regions)[c(1:15)], 
                         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
  
   #scale_fill_gradient_fun = colorRampPalette(c("#202020","#00CCCC"))(100),
  #limit = c(-2,2),

  ##########
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
  ##########
  
  # violin plot
  gathered_pixel_data <- gather(pixel_regions_info_pos, key = "markers", value = "expr", -c("Run", "Mask"))
  ggplot(pixel_regions_info_pos) +
    aes(x = Mask, y = CD47, color = Run) +
    geom_violin(adjust = 1L, scale = "area") +
    #scale_y_continuous(trans = "log") +
    scale_fill_gradient() +
    theme_minimal() +
    facet_wrap(vars(Mask))


  # corr plot
  for (i in c('Ctrl','MedAD','HiAD')) {
         filtered_data_corr <- filter(pixel_regions_info, Run == i, Mask == 4)
         m <- cor(filtered_data_corr[3:ncol(pixel_regions_info)])
         corrplot(m, type="upper", tl.col="black", tl.srt=45, tl.cex = 0.5, insig = "blank", title = paste0(i,5))
     }

#### Synaptic Density Analysis Code ####
  ex_synaptic_markers <- c("PSD95", "Synaptophysin", "VGLUT1", "VGLUT2")
  in_synaptic_markers <- c("PanGAD6567", "VGAT")
  disease_synaptic_markers <- c("PHF1Tau", "Amyloidbeta140", "Amyloidbeta142", "PanAmyloidbeta1724", "CD47") 
  #### basic percents / quants ####
    #synapse subdivisions (excitatory, inhibitory, with or w/o disease markers)
    ex_synapses <- filter(pixel_regions_info, (Synaptophysin > 0 & (VGLUT1 > 0 | VGLUT2 > 0) & (PSD95 > 0)))
    ex_synapses_t <- filter(ex_synapses, PHF1Tau > 0)
    ex_synapses_a <- filter(ex_synapses, Amyloidbeta140 > 0 | Amyloidbeta142 > 0 | PanAmyloidbeta1724 > 0)
    in_synapses <- filter(pixel_regions_info, (Synaptophysin > 0 & (PanGAD6567 > 0 | VGAT > 0)))
    in_synapses_t <- filter(in_synapses, PHF1Tau > 0)
    in_synapses_a <- filter(in_synapses, Amyloidbeta140 > 0 | Amyloidbeta142 > 0 | PanAmyloidbeta1724 > 0)
    
    ex_syn_percent <- table(ex_synapses$Severity) / table(pixel_regions_info$Severity) * 100
    in_syn_percent <- table(in_synapses$Severity) / table(pixel_regions_info$Severity) * 100
    
    ex_syn_t_percent <- table(ex_synapses_t$Severity, ex_synapses_t$`De novo regions`) / table(ex_synapses$Severity, ex_synapses$`De novo regions`) * 100
    in_syn_t_percent <- table(in_synapses_t$Severity, in_synapses_t$`De novo regions`) / table(in_synapses$Severity, in_synapses$`De novo regions`) * 100
    
    ex_syn_a_percent <- table(ex_synapses_a$Severity, ex_synapses_a$`De novo regions`) / table(ex_synapses$Severity, ex_synapses$`De novo regions`) * 100
    in_syn_a_percent <- table(in_synapses_a$Severity, in_synapses_a$`De novo regions`) / table(in_synapses$Severity, in_synapses$`De novo regions`) * 100
    
    syn_percent <- rbind(ex_syn_t_percent, ex_syn_a_percent, in_syn_t_percent, in_syn_a_percent)
    syn_percent <- as.data.frame(syn_percent)
    syn_percent <- cbind(syn_percent, c(rep(c('HiAD','MedAD','Ctrl'), 4)))
    syn_percent <- cbind(syn_percent, c(rep('Excitatory', 6), rep('Inhibitory', 6)))
    syn_percent <- cbind(syn_percent, rep(c(rep('PHF1Tau positive', 3), rep('Amyloid positive', 3)),2))
    names(syn_percent)[6:8] <- c("Severity","Synapse type","Disease type")
    syn_percent$Severity <- ordered(syn_percent$Severity, levels = c("Ctrl","MedAD","HiAD"))
    rownames(syn_percent) <- c(1:nrow(syn_percent))
    syn_percent <- gather(syn_percent, key=`De novo region`, value = 'Percent positive', -c("Severity", "Synapse type", "Disease type"))
    syn_percent$`De novo region` <- factor(syn_percent$`De novo region`, levels = c("ND","GD","MD","PD","TD"))
    levels(syn_percent$`De novo region`) <- list(ND="ND",GD="GD",MD="MD",APD="PD",TTD="TD")
    
    #if just handling deNovo
      ex_syn_percent <- table(ex_synapses$`De novo regions`) / table(pixel_regions_info$`De novo regions`) * 100
      in_syn_percent <- table(in_synapses$`De novo regions`) / table(pixel_regions_info$`De novo regions`) * 100
      
      ex_syn_t_percent <- table(ex_synapses_t$`De novo regions`) / table(ex_synapses$`De novo regions`) * 100
      in_syn_t_percent <- table(in_synapses_t$`De novo regions`) / table(in_synapses$`De novo regions`) * 100
      
      ex_syn_a_percent <- table(ex_synapses_a$`De novo regions`) / table(ex_synapses$`De novo regions`) * 100
      in_syn_a_percent <- table(in_synapses_a$`De novo regions`) / table(in_synapses$`De novo regions`) * 100
      
      syn_percent <- rbind(ex_syn_t_percent, ex_syn_a_percent, in_syn_t_percent, in_syn_a_percent)
      syn_percent <- as.data.frame(syn_percent)
      syn_percent <- cbind(syn_percent, c(rep('Excitatory', 2), rep('Inhibitory', 2)))
      syn_percent <- cbind(syn_percent, rep(c('PHF1Tau positive', 'Amyloid positive'),2))
      names(syn_percent)[6:7] <- c("Synapse type","Disease type")
      rownames(syn_percent) <- c(1:nrow(syn_percent))
      syn_percent <- gather(syn_percent, key=`De novo region`, value = 'Percent positive', -c("Synapse type", "Disease type"))
      syn_percent$`De novo region` <- factor(syn_percent$`De novo region`, levels = c("ND","GD","MD","PD","TD"))
      levels(syn_percent$`De novo region`) <- list(ND="ND",GD="GD",MD="MD",APD="PD",TTD="TD")
    
    #disease positive synaptic pixel percents #USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE#v#USE#
    ggplot(syn_percent) +
      aes(x = `Synapse type`, fill = Severity, weight = `Percent Positive`) +
      geom_bar(position = "dodge") +
      scale_fill_manual(values = color4) +
      theme_minimal(20) +
      facet_wrap(vars(`Disease type`))
    
    ggplot(filter(syn_percent, `Disease type` == 'PHF1Tau positive' & `De novo region` %in% c("ND","MD","TD"))) +
      aes(x = `Synapse type`, fill = `Synapse type`, weight = `Percent positive`) +
      geom_bar(position = "dodge") +
      scale_fill_manual(values = c('blue','purple')) +
      theme_minimal() +
      facet_grid(vars(Severity), vars(`De novo region`))
    
    #try subsampling synaptic pixels
    ex1 <- sample_n(ex_synapses, size = 50)
    ex1_mat <- data.matrix(ex1[c('Mask', 'Run', ex_synaptic_markers, in_synaptic_markers, disease_synaptic_markers)])
    colnames(ex1_mat) = c('Mask', 'Run', ex_synaptic_markers, in_synaptic_markers, disease_synaptic_markers)
    rownames(ex1_mat) = 1:nrow(ex1_mat)
    
    ex1 <- as.data.frame(ex1_mat)
    anno_data = ex1[, c('Mask', 'Run')]
    anno_data = as.data.frame(anno_data)
    #attempt to heatmap all the synpases (probably need to subsample and UMAP them)
    pheatmap(ex1[,c(ex_synaptic_markers,in_synaptic_markers,disease_synaptic_markers)], cluster_rows = T, cluster_cols = T, border_color = FALSE, scale = "column", main = "Synapse expression", 
             annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row =T,
             col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 20, fontsize = 10, width = 25, height = 25,
             filename = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v2/exsynapse_hm.pdf')
    #attempt to plot histograms of synapses --> amount of PHF1Tau
    ex1 <- ex_synapses
    ex1[,3:38] <- ex1[,3:38]+1
    ggplot(ex1) +
      aes(x = PHF1Tau) +
      geom_histogram(bins = 1000L, aes(fill = `De novo regions`)) +
      #scale_fill_manual(values = color3) +
      scale_x_continuous(trans = "log10") +
      #scale_y_continuous(trans = "log10") +
      theme_minimal() +
      facet_grid(row = vars(Severity), col = vars(`De novo regions`))
    
    # use percentages of combo 'synapses' --> heatmap of small regions as put out by neurons in matlab
    
  #### correlation ####
    for (i in c("Ctrl","MedAD","HiAD")) {
      cor_sample <- filter(pixel_regions_info, Severity == i)
      syn_cor <- cor(cor_sample[c(ex_synaptic_markers, in_synaptic_markers, disease_synaptic_markers)])
      # syn_cor_c <- syn_cor
      # syn_cor_c[which(syn_cor == 1)] <- 0
      corrplot(syn_cor, type = "upper", order = "hclust", 
               tl.col = "black", tl.srt = 45, title = i)
    }
    for (i in dnr_names) {
      cor_sample <- filter(pixel_regions_info, `De novo regions` == i)
      syn_cor <- cor(cor_sample[c(ex_synaptic_markers, in_synaptic_markers, disease_synaptic_markers)])
      # syn_cor_c <- syn_cor
      # syn_cor_c[which(syn_cor == 1)] <- 0
      corrplot(syn_cor, type = "upper", order = "hclust", 
               tl.col = "black", tl.srt = 45, title = i)
    }
    
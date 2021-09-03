# Read in RGB mask TIF(s) and plot a heatmap of mean expression values per region (and per sample)

#### Import All ####
  # Read in csv with pixel per mask, run info (generated in MATLAB)
  pixel_regions_info <- read_csv('/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi/separated_masks_voronoi/pixel_info.csv',  col_names = TRUE)
  pixel_regions_info$Run <- factor(pixel_regions_info$Run, levels = c(1,2,3)) # factorize run
  levels(pixel_regions_info$Run) <- c('HiAD', 'MedAD', 'Ctrl')
  pixel_regions_info$Mask <- factor(pixel_regions_info$Mask) # factorize region

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
    dplyr::group_by(Mask, Run) %>%
    dplyr::summarize_if(is.numeric, funs(mean)) %>%
    ungroup() ->
    heat_pixels_regions
  
  heat_pixels_regions_OG <- heat_pixels_regions
  norm_means <- function(x){(x-min(x))/(max(x)-min(x))}
  heat_pixels_regions[3:ncol(heat_pixels_regions)] <- norm_means(heat_pixels_regions[3:ncol(heat_pixels_regions)])
  heat_pixels_regions[3:ncol(heat_pixels_regions)] <- log(heat_pixels_regions[3:ncol(heat_pixels_regions)])

  # Convert data.frame to matrix for making heatmap with pheatmap
  heat_pixels_regions_mat <- data.matrix(heat_pixels_regions[c('Mask', 'Run', nu_panel)])

  # Make the heatmap
  carpet_markers = c("PHF1Tau", "Amyloidbeta140", "Amyloidbeta142", "PanAmyloidbeta1724", "CD47", "PanGAD6567", "Parvalbumin", 
                     "PSD95", "Synaptophysin", "VGAT", "VGLUT1", "VGLUT2", "Calretinin", "Calbindin", "MAP2", "TotalTau", "CD56Lyo", 
                     "MBP", "MAG", "GFAP")
  nu_panel =  setdiff(markers, c('HistoneH3Lyo','EEA1','SERT','TH','Reelin','X8OHGuano'))
  colnames(heat_pixels_regions_mat) = c('Mask', 'Run', nu_panel)
  rownames(heat_pixels_regions_mat) = 1:nrow(heat_pixels_regions_mat)
  
  heat_pixels_regions_mat <- as.data.frame(heat_pixels_regions_mat)
  anno_data = heat_pixels_regions[, c('Mask', 'Run')]
  anno_data = as.data.frame(anno_data)
  
  #for heatmap colors
  mycolors <- list(
    Mask = c('1' = "blue1", '2' = "yellow3", '3' = "green1", '4' = "brown", '5' = "orchid1"),
    Run = c(Ctrl = "dodgerblue2", MedAD = 'green4', HiAD = '#E31A1C'))
  
  pheatmap(heat_pixels_regions_mat[, carpet_markers], cluster_rows = T, cluster_cols = T, border_color = FALSE, scale = "column", main = "Mean expression per de novo region", 
           annotation_row = anno_data, annotation_colors = mycolors, show_rownames = F, annotation_names_row =T, gaps_row = c(7,8),
           col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 30, cellheight = 30, fontsize = 10, width = 25, height = 10,
           filename = '/Volumes/BryJC_Stanford/paper1_final/Fig5/v2/carpet_hm.pdf')
  
  #test
  filtered_data <- filter(heat_pixels_regions)
  heatmaply(filtered_data[,c(carpet_marker)], scale = 'column', Rowv = FALSE, Colv = FALSE, labRow = rownames(heat_pixels_regions)[c(1:5)], 
            col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
  
  heatmaply(filtered_data[,c('Run','Mask',carpet_markers)], scale = 'column', Rowv = TRUE, Colv = TRUE, labRow = rownames(heat_pixels_regions)[c(1:15)], 
                         col= colorRampPalette(c("#202020","#00CCCC"))(100), cellwidth = 10, cellheight = 10)
  
   #scale_fill_gradient_fun = colorRampPalette(c("#202020","#00CCCC"))(100),
  #limit = c(-2,2),


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
  
  #synapse subdivisions (excitatory, inhibitory, with or w/o disease markers)
  ex_synapses <- filter(pixel_regions_info, PSD95 > 0 & (Synaptophysin > 0 | VGLUT1 > 0 | VGLUT2 > 0))
  ex_synapses_t <- filter(pixel_regions_info, PHF1Tau > 0 & (PSD95 > 0 & (Synaptophysin > 0 | VGLUT1 > 0 | VGLUT2 > 0)))
  ex_synapses_a <- filter(pixel_regions_info, (Amyloidbeta140 > 0 | Amyloidbeta142 > 0 | PanAmyloidbeta1724 > 0) & (PSD95 > 0 & (Synaptophysin > 0 | VGLUT1 > 0 | VGLUT2 > 0)))
  in_synapses <- filter(pixel_regions_info, PanGAD6567 > 0 & VGAT > 0)
  in_synapses_t <- filter(pixel_regions_info, PHF1Tau > 0 & (PanGAD6567 > 0 & VGAT > 0))
  in_synapses_a <- filter(pixel_regions_info, (Amyloidbeta140 > 0 | Amyloidbeta142 > 0 | PanAmyloidbeta1724 > 0) & (PanGAD6567 > 0 & VGAT > 0))
  
  ex_syn_percent <- table(ex_synapses$Run) / table(pixel_regions_info$Run) * 100
  in_syn_percent <- table(in_synapses$Run) / table(pixel_regions_info$Run) * 100
  
  ex_syn_t_percent <- table(ex_synapses_t$Run) / table(ex_synapses$Run) * 100
  in_syn_t_percent <- table(in_synapses_t$Run) / table(in_synapses$Run) * 100
  
  ex_syn_a_percent <- table(ex_synapses_a$Run) / table(ex_synapses$Run) * 100
  in_syn_a_percent <- table(in_synapses_a$Run) / table(in_synapses$Run) * 100
  
  syn_percent <- rbind(ex_syn_t_percent, ex_syn_a_percent, in_syn_t_percent, in_syn_a_percent)
  rownames(syn_percent) = c('Excitatory synaptic pixels - Tau Positive', 'Excitatory synaptic pixels - Amyloid Positive', 'Inhibitory synaptic pixels - Tau Positive', 'Inhibitory synaptic pixels - Amyloid Positive')
  syn_percent <- as.data.frame(syn_percent)
  syn_percent <- gather(syn_percent)
  syn_percent <- cbind(syn_percent, c(rep(c('Excitatory synaptic pixels', 'Excitatory synaptic pixels', 'Inhibitory synaptic pixels', 'Inhibitory synaptic pixels'), 3)))
  syn_percent <- cbind(syn_percent, c(rep(c('Tau Positive', 'Amyloid Positive'), 6)))
  colnames(syn_percent) <- c('Severity','Percent Positive','Synapse type','Disease type')
  syn_percent$Severity <- ordered(syn_percent$Severity, levels = c("Ctrl","MedAD","HiAD"))
  #disease positive synaptic pixel percents #USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE##USE#v#USE#
  ggplot(syn_percent) +
    aes(x = `Synapse type`, fill = Severity, weight = `Percent Positive`) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = color4) +
    theme_minimal(20) +
    facet_wrap(vars(`Disease type`))
  
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
  ex1 <- ex_synapses+1
  ggplot(ex_synapses) +
    aes(x = PHF1Tau) +
    geom_histogram(bins = 100L, aes(fill = Mask)) +
    scale_fill_hue() +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    theme_minimal() +
    facet_wrap(vars(Run))
  
  # use percentages of combo 'synapses' --> heatmap of small regions as put out by neurons in matlab
  
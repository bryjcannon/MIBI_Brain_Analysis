install.packages("renv")
library("renv")
renv::init("/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/single_object_pipeline")
setwd("/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/single_object_pipeline")

install_ez_packages(T)
load_ez_packages(T)

#Seed
seed <- 123 
set.seed(seed)

#reset plots
dev.off()

#full panel
panel <- c("C12",                "Na23",              
           "Si28",               "Ca40",               "HistoneH3Lyo",       "CD47",               "TotalTau",           "Background",         "empty139",           "ApoE4",             
           "Amyloidbeta140",     "Amyloidbeta142",     "Calretinin",         "CD56Lyo",            "Parvalbumin",        "PanAmyloidbeta1724", "CD31",               "MAP2",              
           "PolyubiK48",         "MFN2",               "PanGAD6567",         "VGAT",               "SERT",               "PolyubiK63",         "Reelin",             "EEA1",              
           "pTDP43",             "Synaptophysin",      "CD45",               "PanApoE2E3E4",       "X8OHGuano",          "TH",                 "Presenilin1NTF",     "MAG",               
           "VGLUT1",             "VGLUT2",             "Calbindin",          "Iba1",               "PSD95",              "MBP",                "MCT1",               "CD33Lyo",           
           "GFAP",               "CD105",              "PHF1Tau",            "Ta181",              "Au197")

#panel_notes <- [, 4:44][, -c(4,5,7,9,11,15,16,21,29,33,38)]) 

# name, x, y, direction --> for tiling script
  #tile_dimensions <- c(('HiResHi',5,17,1,'CA2 then DG'), ('Ctrl',15,20,1), ('Hi',14,14,1) ('Med',20,13,1))
  #number_pts_per_run <- c(('Hi',1:196), ('Med',0:259), ('Ctrl1',1:100), ('Ctrl2',1:105), ('Ctrl3',1:70)) 

# De Novo Region Names --> c('Neuron Dominant','Glia Dominant','Plaque Dominant','Tangle/Thread Dominant','Mixed Disease')
  dnr_names <- c('ND','GD','MD','APD','TTD')
  names(dnr_names) <- c(1,2,5,3,4)

#AD3 [, -c(6,8,9,11,13,15,19,20,25,33,37,42)]) 

# Pixel information for medRes runs
pix_l_medRes <- 500/512 #length in um
pix_a_medRes <- pix_l_medRes**2 #area in um**2
pixel_nums <- table(pixel_regions_info$Severity) #total pixels for each run (using only masked pixels from cleaned voronoi expansion)
pixel_anat_area <- as.data.frame(table(pixel_anat_info$Run, pixel_anat_info$Mask) * pix_a_medRes)  #total area for each subanatomical region within each run, in mm**2
  names(pixel_anat_area) <- c('Run', 'Mask', 'area')
  pixel_anat_area$area <- pixel_anat_area$area/1e6
run_areas <- (pixel_nums*pix_a_medRes) / 1e6 # calculates mm**2 area
run_areas <- run_areas[c('Ctrl', 'MedAD', 'HiAD')]

#### Figure display Schemes ####
  
  #Color palette information

  install.packages("viridis")
  library(viridis)

  # color palette taken fro stackOverflow <- https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
  c25 <- c(
    "#E31A1C", # red
    "#FF7F00", # orange
    "green4",
    "#6A3D9A", # purple
    "dodgerblue2", 
    "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown", "black")
  pie(rep(1, 25), col = palette)

  #Cell/Object Types -> c('tangles', 'plaques', 'microglia', 'endothelial', 'neurons', 'non-immune glia') 
  color1 <- c25[c(1:6)]
  #Cell/Object Types KV -> c('tangles-threads', 'plaques', 'microglia', 'endothelial', 'neurons', 'non-immune glia')
    # old tangle color --> #6409DB
  colorKV1 <-  c("#00FFFF", "#FF007F", "#FFBF00", "#0028D4", "#980A0A", "#106900")
  #46B2AC, CIND - #D1EF54, CID - #E22D43 #Anatomical Regions -> c('DG', 'CA4', 'CA3', 'CA2', 'CA1') 
  color2 <- c25[c(14,13,11,10,9)]
  names(color2) <- c('DG', 'CA4', 'CA3', 'CA2', 'CA1') 
  #De novo regions -> c(1,2,3,4,5) == c(ND,GD,APD,TTD,MD) == c(N,G,AP,TT,M)
    # old colors -> c25[c(17,22,20,24,15)] 
  color3 <- c("#EE00FF","#00CD25","#00D1FF","#6D2BA4","#FFCC84")
  names(color3) <- c(1,2,3,4,5)
  names(color3) <- c("ND","GD","APD","TTD","MD")
  #Samples -> c('Ctrl', 'Med', 'HiAD')
    # old sample colors c25[c(5,3,1)]
  color4 <- c("#46B2AC", "#D1EF54", "#E22D43")
  names(color4) <- c('Ctrl','MedAD','HiAD')
  #Expression color spectrum
  colorExp <- colorRampPalette(c("#202020","#00CCCC"))(100)
  #Color disease positive cells
  color5 <- c('blue', 'red')
  names(color5) <- c(0,1)
  #color neuronal subpops
  color6 <- c(c25[c(5:11)],"#002B03",c25[c(14:16)])
  names(color6) <- c(1,2,3,4,5,6,7,8,9,10,11)
  #color neuronal subprop groups ('patho disease' vs 'no-patho disease')
  color7 <- c('red','red','blue','blue','blue','red','blue','red','blue','red','red')
  names(color7) <- c(1,2,3,4,5,6,7,8,9,10,11)
  # marker colors
    # PHF1-TAU (#00FFFF), MFN2 (#0080FF), PolyubiK48 (#FFFF00), Iba1 (#FFBF00), ApoE (#7F3F98) and Ab42 (#FA02A5), CD45 (#3db088) K48 be difficult to see with the K48 yellow, so you can use (#E960FF) for the K48..
  
  #ggplot themes for bar, linegraphs
  #broad names
  theme1 <- theme(plot.title = element_text(size=22, hjust = 0.5),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"),
          strip.text = element_text(size=14),
          axis.text.x = element_text(face="bold", size=14), 
          axis.text.y = element_text(face="bold", size=14),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12, face="bold"))
  
  theme1a <- theme(plot.title = element_text(size=34, hjust = 0.5),
                  axis.text = element_text(size=26),
                  axis.title = element_text(size=26,face="bold"),
                  strip.text = element_text(size=26),
                  axis.text.x = element_text(face="bold", size=26), 
                  axis.text.y = element_text(face="bold", size=26),
                  axis.ticks = element_line(),
                  legend.title = element_text(size=26, face="bold"),
                  legend.text = element_text(size=26, face="bold"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  #blank slate with some titles
  theme2 <- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), 
                  plot.title=element_text(colour="black", size=10, face="bold"),
                  #axis.ticks=element_blank(), axis.text=element_blank(), 
                  axis.title=element_text(colour="black", size=10, face="bold"), legend.text = element_text(colour="black", size=10, face="bold"), 
                  legend.position = "top")
  #blank slate with almost nothing
  theme3 <- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), 
                  plot.title=element_blank(),
                  axis.ticks=element_blank(), axis.text=element_blank(), 
                  axis.title=element_blank(), legend.text = element_blank())
  #blank slate with almost nothing
  theme4 <- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
                  plot.title=element_blank(),
                  axis.ticks=element_blank(), axis.text=element_blank(), 
                  axis.title=element_blank(), legend.text = element_blank())
  
  #ggsave template
  ggsave(filename = 'line_composition_abs_DGtoCA1', plot = B2v2, device = 'pdf',
         path = '/Volumes/BryJC_Stanford/paper1_analysis/Fig6/plots/fig_plots_final_set1/linegraphs',
         scale = 1, width = 3000, height = 2055, limitsize = FALSE,
         dpi = "retina")
  
  #for setting fonts completely up
  theme_bw(25)
  
    
#### plotting ####
  plotMarkers <- function(origdat, markerset, path = './'){
    dat <- filter(origdat, UMAP1 > -20 & UMAP2 > -20)
    
    ggplot(dat, aes(x = UMAP1, y = UMAP2, color=de_novo_regions)) +
      geom_point(size = 0.01) + 
      coord_fixed(ratio = 1)
    picname <- paste0(path,'de_novo_regions')
    ggsave(file=paste(picname,"png",sep="."),width = 7, height = 7, dpi = 320)
    
    ggplot(dat, aes(x = UMAP1, y = UMAP2, color=FlowSOM_ids)) +
      geom_point(size = 0.01) + 
      coord_fixed(ratio = 1)
    picname <- paste0(path,'FlowSOM_ids')
    ggsave(file=paste(picname,"png",sep="."),width = 7, height = 7, dpi = 320)
    
    for (marker in markerset){
      curr <- ggplot(dat, aes_string(x = 'UMAP1', y = 'UMAP2', color=marker)) +
          geom_point(size = 0.01) + 
          scale_colour_gradient(low='blue', high='red') +
          coord_fixed(ratio = 1)
      picname <- paste0(path,marker)
      ggsave(curr, file=paste(picname,"png",sep="."),width = 7, height = 7, dpi = 320)
    }
  } 
    
#### converting rgb to single integer ####
  col2rgb() #for converting hex/name to rgb
  
  rgb2int <- function(color){
    x = (256^2)*col2rgb(color)[1] + 256*col2rgb(color)[2] + col2rgb(color)[3]
    return(x)
  }
  
  #grab a collection of color conversions
  for (i in MFG_masks){
    print(rgb2int(i))
    }
  
  
#### gg rastering (need ggrastr) ####
  theme <- theme(strip.background = element_blank(),
                 panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  ggplot(data_myeloid.melt.rand, aes(x = umap1, y = umap2, color = value)) +
    geom_point_rast(size=3) +
    scale_color_viridis(option = 'magma', limits=c(0, 1)) +
    facet_wrap(~ variable, ncol = 4) +
    theme +
    theme(legend.position = 'none')
  
#### Instructions for adding overlapping / function for additional gating based on overlapping####
  # need to produce RGB images with unique spatial ids, import into matlab and convert back to unique spatial ids after tiling
  # overlay the objects of interest and use associated script to produce map and gating for cells with object overlaps
  # import into R using tiled_obj_import and then use below to add additional gates
  
  #adding additional gates upon overlaps
  
  #for each id or phenotype column
  #if id or phenotype asked to be covered, if tau+ or amy+, replace current id or phenotype with 
#### Plotting notes
  #dpi originally height = 6 for set 2 variants (DG to CA1)
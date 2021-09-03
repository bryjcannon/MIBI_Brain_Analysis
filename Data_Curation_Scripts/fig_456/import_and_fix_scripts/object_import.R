##### object_import ##### - Repeat for Hi, Med, and Ctrl Conditions
  # Import post-ez_segmenter data or deepcell data (MATLAB generated) for single object + cell analysis, respectively
  # Author: Bryan Cannon 2020 (multiple code snippets taken or built from DM, FH, EFM, DT)

##### INSTALL / LOAD PKGS #####
  # install and / or load ez_pkgs if you haven't already done so (will also install color scheme)
  # to do this make sure you go to ez_lib script and source it. below functions will then work
  source('/Volumes/BryJC_Stanford/ez_segmenter_analysis/packages_and_functions/ez_lib.R')
  install_ez_packages(F)
  load_ez_packages(T)
  
    # set seed for downstream analysis
  seed <- 123 
  set.seed(seed)
  
  # move to location of your data, load in csv and .mat files (spatial information stored in .mat)
  data_folder <- '/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_MedADCase/NoAuBGFFtDenoised_New'
  setwd(data_folder)

##### 1) LOAD ezSeg DATA ##### 
  # specify run folders
  ez_runs <- c('ezSegResults_MedRes_MedAD')
  run_names = c("MedAD") # needs to be same length as ez_runs
  object_types <- c('amyloid-pathies', 'tau-pathies', 'microglia-processes', 'vessels')
  
  # create data.frame containing single object data from csv files - for each run and each object type
  master_obj_data <- data.frame()
  
  for (run_index in 1:length(ez_runs)){
    obj_data_raw_all <- data.frame()
    
    for (obj_type in object_types) {
      csv_names <- list.files(path = paste0(ez_runs[run_index],'/','objects_points'), recursive = T, full.names = T, pattern = paste0(obj_type, "_dataScaleSize.csv")) # read in csv files for object type
      csv_names <- mixedsort(csv_names)
      
      obj_data_raw <- lapply(csv_names, read.csv) %>% bind_rows() # grab data from csv's then convert data to data.frame
      
      obj_type_id <- rep(obj_type, dim(obj_data_raw)[1]) # create column of length object number with obj_type info
      obj_data_raw <- cbind(obj_data_raw, obj_type_id) # add obj_type_id to data
      
      obj_data_raw_all <- rbind(obj_data_raw_all, obj_data_raw) # collate to run data.frame
    }
    
    run_type_id <- rep(run_names[run_index], dim(obj_data_raw_all)[1]) # create column of length object number with run_type info, name run names to reflect actual regions
    obj_data_raw_all <- cbind(obj_data_raw_all, run_type_id) # add run_type_id to data
    
    master_obj_data <- rbind(master_obj_data, obj_data_raw_all) # collate to master data.frame
    
    rm(obj_data_raw)
    rm(obj_data_raw_all)
  }
  
##### 2) LOAD Deep Cell DATA #####
  # move to location of your data, load in fcs and .mat files (spatial information retained here)
  deepcell_runs <- c('segmentationDC')
  
  # create FlowSet from fcs files (eventually will want to adapt for csv read-in as well)
  
  master_cell_data <- data.frame()
  
  for (run_index in 1:length(deepcell_runs)) {
    fcs_names <- c(list.files(path = paste0(deepcell_runs[run_index],'/','fcs_cells'), full.names = T, pattern = "(dataScaleSizeFCS.*fcs)"))
    cell_flowSet <- read.flowSet(files = fcs_names, alter.names = T, transformation = FALSE, emptyValue = FALSE, truncate_max_range = FALSE)
    cell_data_matrix <- fsApply(cell_flowSet, Biobase::exprs) #https://support.bioconductor.org/p/109128/ --> explains why use Biobase::exprs
    cell_data_raw <- as.data.frame(cell_data_matrix)
    
    # add columns to each file denoting point number, type of tissue, name of objects
    point_names <- sampleNames(cell_flowSet)
    point_numbers <- as.numeric(gsub("[^0-9]*", '', gsub(".*p[^0-9]*", '', point_names))) # grabs point numbers from filenames and converts to numeric
    point_id <- as.vector(x = NULL)  # create a list 'point_source' that assigns a point_id to each object
    
    dim = fsApply(cell_flowSet, dim)
    dim = as.numeric(dim[,1])
    for(i in 1:length(dim)) { #loop creates the actual vector with point id's for each object
      temp_point_id <- rep(point_numbers[i], dim[i])
      point_id <- c(point_id, temp_point_id)
    }
    cell_data_raw <- cbind(cell_data_raw, point_id) # add point_id to data
    
    
    
    master_cell_data <- rbind(master_cell_data, cell_data_raw)
    
    rm(cell_data_matrix)
    rm(cell_flowSet)
    rm(cell_data_raw)
  }
  
##### 3) COMBINE Deep Cell and ezSeg Data #####
  # rename deepCell labels to match ez_segmenter labels
  master_cell_data <- plyr::rename(master_cell_data, c('cellLabelInImage' = 'obj_id', 'cellSize' = 'obj_size', 'cell_type_id' = 'obj_type_id')) # rename deep cell columns to match ez names
  master_cell_data <- master_cell_data %>% select(point_id, everything()) # move point_id to start of deep cell data to match ez column order
  
  master_cell_data <- plyr::rename(master_cell_data, c('MFN2.' = 'MFN2', 'VGAT.' = 'VGAT')) # rename markers to handle error from original csv in cell extraction for these runs
  
  # row bind the deepCell dataFrame data with ezSeg dataFrame data
  master_data <- rbind(master_cell_data, master_obj_data[,-c(51:53)]) # excluding composite channels from imported ez data
  
##### 4) TRANSFORMATION options #####
  
  # ALWAYS assign and standardize panels
  panel_start = 4
  panel_end = 50
  # use names(data) to figure out the indices of which columns to ignore (metals, empty, background, etc), keeping only the markers in panel
  panel <- names(master_data[ ,panel_start:panel_end]) # -c(#,#,#): remove (!metals), composites, other labels after already paring down columns (i.e. numbers will change)
  
  # Normalize all objects by total counts within each cell or object
  master_data_all_ionNorm <- master_data_all %>% mutate(totalIon = rowSums(.[panel_start:panel_end]))
  master_data_all_ionNorm[,panel_start:panel_end] <- master_data_all_ionNorm[,panel_start:panel_end] / master_data_all_ionNorm[,"totalIon"]
  # using temporarily for ionNorm look
  master_data_ion_cuberoot <- master_data_all_ionNorm
  master_data_ion_cuberoot[,panel_start:panel_end] <- (master_data_ion_cuberoot[,panel_start:panel_end]*100000)^(1/3)
    
  # 10k multiply + cuberoot transform (need to standardize) --> FINAL CHOICE part1
  master_data_all_cuberoot <- master_data
  master_data_all_cuberoot[panel] <- (master_data_all_cuberoot[panel]*10000)^(1/3)
  
  # log+1 transformation
  data_log1_transform <- master_data
  data_log1_transform[,panel_start:panel_end] = log(data_log1_transform[,panel_start:panel_end]+1)
  
  # linear transformation
  data_linear_transform <- data_log1_transform
  data_linear_transform[,panel_start:panel_end] <- data_linear_transform[,panel_start:panel_end]*100 # multiply all counts by 100 (linear transform)
  
  # do arcsinh transformation only for the clustering.channels
  data_lin_asinh_transf <- data_linear_transform
  asinh_scale <- 5
  data_lin_asinh_transf[,panel_start:panel_end] <- asinh(data_lin_asinh_transf[,panel_start:panel_end] / asinh_scale)
  
  
  
  
##### 5) FINALLY use this collated dataset for downstream work
  master_data_Hi <- master_data
  master_data_Hi_cuberoot <- data_log1_transform
  rm(master_data)
  rm(master_cell_data)
  rm(master_obj_data)
  rm(data_log1_transform)
  
  # run only after all three datasets have been collected
  master_data_all <- rbind(master_data_Hi, master_data_Med, master_data_Ctrl)
  master_data_all_cuberoot <- rbind(master_data_Hi_cuberoot, master_data_Med_cuberoot, master_data_Ctrl_cuberoot)
  
  # PERCENTILE normalize expression values from 0 to 1 --> FINAL CHOICE part2
  data_normalized <- master_data_all_cuberoot
  normalization_vector <- apply(master_data_all_cuberoot[,panel_start:panel_end], 2, function(x) quantile(x, 0.9999, names = F))
  data_normalized[,panel_start:panel_end] <- t(t(data_normalized[,panel_start:panel_end]) / as.numeric(normalization_vector))
  # check whether you adjusted the range approximately from 0 to 1
  apply(data_normalized[,panel_start:panel_end], 2, max)
  

    
  
  
#### create data.frame containing neighborhood data from csv files - for each run ####
cytomap_files_path <- "/Volumes/BryJC_Stanford/paper1_analysis/Fig6/single_object_pipeline/cytomap_files/Analysis_FSOM3_v3"
setwd(cytomap_files_path)
conditions <- c('Ctrl', 'MedAD', 'HiAD')

master_neighborhood <- data.frame()

for (run_index in 1:length(conditions)){
  neighborhood_raw <- data.frame()
  
  csv_names <- list.files(full.names = T, pattern = paste0(conditions[run_index], "_neighborhoods.csv")) # read in csv files for condition
  
  neighborhood_raw <- lapply(csv_names, read.csv) %>% bind_rows() # grab data from csv's then convert data to data.frame
  
  run_type_name <- rep(conditions[run_index], dim(neighborhood_raw)[1]) # create column of length object number with run_type info, name run names to reflect actual regions
  neighborhood_raw <- cbind(neighborhood_raw, run_type_name) # add run_type_id to data
  
  master_neighborhood <- rbind(master_neighborhood, neighborhood_raw) # collate to master data.frame
  
  rm(csv_names)
  rm(run_type_name)
  rm(neighborhood_raw)
}
  #get rid of Ch_ abbreviation
  names(master_neighborhood)[7:72] <- substring(names(master_neighborhood)[7:72], 4)
  
  #alter gate names - rename gated to cell and object ids
  names(master_neighborhood)[74:79] <- c('endothelial','microglia','neurons','non-immune glia','plaques','tangles')

#### import region annotations for each cell into new data.frame ####
master_regions <- data.frame()

for (run_index in 1:length(conditions)){
  obj_raw <- data.frame()
  
  csv_names <- list.files(full.names = T, pattern = paste0(conditions[run_index], "_obj_w_regions.csv")) # read in csv files for condition
  obj_raw <- lapply(csv_names, read.csv) %>% bind_rows() # grab data from csv's then convert data to data.frame
  
  master_regions <- rbind(master_regions, obj_raw) # collate to master data.frame
  
  rm(csv_names)
  rm(obj_raw)
}
  #get rid of Ch_ abbreviation
  names(master_regions)[4:69] <- substring(names(master_regions)[4:69], 4)
  
  #adjust gate names
  # c('Ctrl_DG','Ctrl_CA4','Ctrl_CA3','Ctrl_CA2','Ctrl_CA1')
  # c('MedAD_DG','MedAD_CA4','MedAD_CA3','MedAD_CA2','MedAD_CA1')
  # c('HiAD_CA1','HiAD_DG','HiAD_CA4','HiAD_CA2','HiAD_CA3')
  names(master_regions)[78:82] <- c('DG','CA4','CA3','CA2','CA1')
  
  #check names
  names(master_regions)
  
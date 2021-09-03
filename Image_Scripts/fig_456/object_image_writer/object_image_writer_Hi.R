##### object_image_writer ##### - Repeat for Hi, Med, and Ctrl Conditions
  # Working with post-ez_segmenter data ann deepcell (MATLAB GUI generated) for single object image creation
  # Author: Bryan Cannon 2020 (multiple code snippets taken or built from DM, FH, EFM, DT)

  # install and / or load ez_pkgs if you haven't already done so (will also install color scheme)
  # to do this make sure you go to ez_lib script and source it. below functions will then work
  install_ez_packages(F)
  load_ez_packages(T)

##### FILE & DATA SETTINGS - user must adjust these! #####
  # run on data that has been already read and transformed into R using ezAnalysis script (script will output all_data)
  data_for_images <- master_gated_regions
  # enter head folder where your ez data is stored
  data_folder = '/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_HiADCase/'
  setwd(data_folder)
  # name to give this specific image creation run when saving
  save_name <- 'DGtoCA1'
  # shorthand for regions or runs, used in run_type_id
  run_names <- c('HiAD')
  # object types you plan to cluster)
  object_types <- c('tangles', 'plaques', 'cells')
  # boolean denoting if clustering was performed on the objects (T) or not (F)
  objects_clustered <- F
  ignore_clusters <- c()
  # boolean denoting if anantomical region labeling was performed on the objects (T) or not (F)
  anatomy_clustered <- T
  ignore_anatomy <- c()
  # boolean denoting if this image is to be used for creating an image that will be used for tieing together
  # unique objects to spatial locations in tiled image (or labelled annotations in the future)
  objects_for_spacing <- F
  # color palete of choice
  palette <- c25
  # the resolution of your images, e.g. 512 x 512, 1024 x 1024, etc.
  resolution_dim <- c(512, 512)
  # list of actual point numbers (i.e. which points you actually want to image, e.g. c(4,6:9) for Points 4 and 6 through 9)
  point_list <- list(1:196)
  
  ## object location information - path information below. ##
  # runs, i.e. regions scanned
  ez_runs <- c('NoAuBGFFtDenoised_New/ezSegResults_MedRes_HiAD')
  # where the objects / data / csv / mat files are stored
  ez_data_container <- 'objects_points'
  # deepcell info, similar structure to counterparts above.
  deepcell_runs <- c('NoAuBGFFtDenoised_New/segmentationDC')
  # where the cells / data / csv / mat files are stored
  deepcell_data_container <- 'single_cell_dynamic_expansion'
  

##### IMAGE CONSTRCTION PROCESS - run entire code block #####
  # for each run, for each point, for each object, create and save colored images
  for (i in 1:length(run_names)) {
    region = run_names[i]
    
    for (p in point_list[[i]]) {
      print(paste0("Processing region: ", region, ", Point: ", p))
      
      for (obj_index in 1:length(object_types)) {  
        # read matlab matrix of specific object type from specific point into R, if no entry create empty mask (note: cell locations obtained from deep cell newLmod.csv)
        tryCatch({
          if (object_types[obj_index] == 'cells') {
            obj_properties <- as.matrix(read.csv(paste0(data_folder, '/', deepcell_runs[i], '/', deepcell_data_container, '/Point', p, '/', 'newLmod.csv'), header = FALSE))
            # pull out mapped object ids as a matrix (defaults to image dimensions)
            obj_mask <- obj_properties
            # set 1's (used in earlier watershedding) to 0 (the actual background value)
            obj_mask[obj_mask == 1] = 0
          }
          else {
            obj_properties <- readMat(paste0(data_folder, '/', ez_runs[i], '/', ez_data_container, '/Point', p, '/', object_types[obj_index], '_objData.mat'))
            # pull out mapped object ids as a matrix (defaults to image dimensions)
            obj_mask <- obj_properties$mapped.obj.ids
          }
        }, error = function(err) {
          obj_mask <- zeros(resolution_dim[1], resolution_dim[2])
        }) # end of tryCatch loop
        
        # create three channels (red, green, blue) for the mask to use in later Image creation - inits with empty matrix size of the image
        obj_mask_r <- zeros(resolution_dim[1], resolution_dim[2])
        obj_mask_g <- zeros(resolution_dim[1], resolution_dim[2])
        obj_mask_b <- zeros(resolution_dim[1], resolution_dim[2])
        
        ##### for each object in point p we replace object ID with colour based on meta number. #####
        pointdata <- filter(data_for_images, run_type_id == region & point_id == p & obj_type_id == object_types[obj_index])
        for (uniq_object in unique(pointdata$obj_id)){
          if (objects_clustered == T) {pointdata <- filter(pointdata, !cluster_id %in% ignore_clusters)}
          if (anatomy_clustered == T) {pointdata <- filter(pointdata, !anat_id %in% ignore_anat)}
          
          # color objects based on clustering ids
          if (objects_clustered == T) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, subset(pointdata, obj_id == uniq_object)$cluster_id]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, subset(pointdata, obj_id == uniq_object)$cluster_id]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, subset(pointdata, obj_id == uniq_object)$cluster_id]
          }
          # produce grayscale masks of cells/objects by anatomical region
          else if (anatomy_clustered == T) {
            obj_mask_r[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$anat_id
          }
          # color objects based on unique RGB_ids calculated earlier
          else if (objects_for_spacing == T) {
            obj_mask_r[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$R
            obj_mask_g[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$G
            obj_mask_b[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$B
          }
          # if no clustering, then use default palette numbers to color objects
          else if ((objects_clustered == F) & (objects_for_spacing == F)) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, obj_index]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, obj_index]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, obj_index]
          }
        }
        # convert palette numbers to 8bit values 
        obj_mask_r = obj_mask_r / 255
        obj_mask_g = obj_mask_g / 255
        obj_mask_b = obj_mask_b / 255
        
        ##### create actual image and save #####
        type_of_image <- ''
        if (objects_clustered == T) {type_of_image <- 'clustered'}
        else if (objects_for_spacing == T) {type_of_image <- 'spacing_unique_ids'}
        else if (regions_clustered == T) {type_of_image <- 'anatomy_labels'}
        else {type_of_image <- 'basic'}
        
        # convert arrays to RGB image if needed
        if (regions_clustered == F) {
          img <- EBImage::transpose(rgbImage(Image(obj_mask_r), Image(obj_mask_g), Image(obj_mask_b)))
        }
        else {
          img <- Image(obj_mask_r)
        }
        
        dir.create(paste0(data_folder, "obj_images"), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i]), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image, "/Point", p), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image, "/Point", p, "/TIFs/"), showWarnings = FALSE)
        writeImage(img, paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image, "/Point", p, "/TIFs/", object_types[obj_index],".tif"))
      }
    }
  }

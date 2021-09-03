##### object_image_writer ##### - Repeat for Hi, Med, and Ctrl Conditions
  # Working with post-ez_segmenter data ann deepcell (MATLAB GUI generated) for single object image creation
  # Author: Bryan Cannon 2020 (multiple code snippets taken or built from DM, FH, EFM, DT)

  # install and / or load ez_pkgs if you haven't already done so (will also install color scheme)
  # to do this make sure you go to ez_lib script and source it. below functions will then work
  install_ez_packages(F)
  load_ez_packages(T)

##### FILE & DATA SETTINGS - user must adjust these! #####
  # run on data that has been already read and transformed into R using ezAnalysis script (script will output all_data)
  data_for_images <- master_regions_img
  # enter head folder where your ez data is stored
  data_folder = '/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/'
  setwd(data_folder)
  # name to give this specific image creation run when saving
  save_name <- 'gray_regions'
  # shorthand for regions or runs, used in run_type_id
  run_names <- c('Ctrl1', 'Ctrl2', 'Ctrl3')
  # rgb or grayscale
  grayscale <- T
  # object types you plan to cluster)
  object_types <- c('cells', 'amyloid-pathies', 'tau-pathies')
  # boolean denoting if clustering was performed on the objects (T) or not (F)
  objects_clustered <- T
  clustering_variable <- 'FlowSOM_ids'
  include_clusters <- c(1, 2, 3, 4, 5) #represents c('tangles','plaques','microglia','endothelial','neurons','non-immune glia') OR de_novo_regions 1-5
  # boolean denoting if anantomical region labeling was performed on the objects (T) or not (F)
  anatomy_clustered <- F
  anantomy_variable <- 'anat_ids'
  include_anatomy <- c(1, 2, 3, 4, 5) #represents c('DG', 'CA4', 'CA3', 'CA2', 'CA1')
  # boolean denoting if this image is to be used for creating an image that will be used for tieing together
  # unique objects to spatial locations in tiled image (or labelled annotations in the future)
  spatial_id_assign <- F
  # color palete of choice, use seond line to see color choices
  palette <- colorKV1
  #To check colors use -> 
  if (grayscale == F) {pie(rep(1, 6), col = palette)}
  # the resolution of your images, e.g. 512 x 512, 1024 x 1024, etc.
  resolution_dim <- c(512, 512)
  # list of actual point numbers (i.e. which points you actually want to image, e.g. c(4,6:9) for Points 4 and 6 through 9)
  point_list <- list(c(1:100), c(1:105), c(1:70))
  
  ## object location information - path information below. ##
  # runs, i.e. regions scanned
  ez_runs <- c('NoAuBGFFtDenoised_New_(Set1)/ezSegResults_Ctrl1', 'NoAuBGFFtDenoised_New_(Set2)/ezSegResults_Ctrl2', 'NoAuBGFFtDenoised_New_(Set3)/ezSegResults_Ctrl3')
  # where the objects / data / csv / mat files are stored
  ez_data_container <- 'objects_points'
  # deepcell info, similar structure to counterparts above.
  deepcell_runs <- c('NoAuBGFFtDenoised_New_(Set1)/segmentationDC', 'NoAuBGFFtDenoised_New_(Set2)/segmentationDC', 'NoAuBGFFtDenoised_New_(Set3)/segmentationDC')
  # where the cells / data / csv / mat files are stored
  deepcell_data_container <- 'single_cell_dynamic_expansion'
  

##### IMAGE CONSTRUCTION PROCESS - run entire code block #####
  # for each run, for each point, for each object, create and save colored images
  for (i in 1:length(run_names)) {
    sample = run_names[i]
    
    for (p in point_list[[i]]) {
      print(paste0("Processing sample: ", sample, ", Point: ", p))

      for (obj_index in 1:length(object_types)) { 
        ##### read in spatial coordinate data ####
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
        
        ##### create an empty matrix to fill in pixel values ####
        # create three channels (red, green, blue) for the mask to use in later Image creation - inits with empty matrix size of the image
        if (grayscale == F) {
          obj_mask_r <- zeros(resolution_dim[1], resolution_dim[2])
          obj_mask_g <- zeros(resolution_dim[1], resolution_dim[2])
          obj_mask_b <- zeros(resolution_dim[1], resolution_dim[2])
        }
        # create single channel for mask - use for grayscale image output
        else if (grayscale == T) {
          obj_mask_gray <-zeros(resolution_dim[1], resolution_dim[2])
        }
          
        ##### for each object in point p we gather relevant objects based on clustering, anatomical, and type data #####
        pointdata <- filter(data_for_images, run_type_id == sample & point_id == p & obj_type_id == object_types[obj_index])
        if (objects_clustered == T) {pointdata <- filter(pointdata, FlowSOM_ids %in% include_clusters)}
        if (anatomy_clustered == T) {pointdata <- filter(pointdata, anat_id %in% include_anatomy)}
        
        ##### for each object in point p we replace object ID with colour based on meta number. #####
        for (uniq_object in unique(pointdata$obj_id)){
          
          # color objects based on clustering ids (need to replace name of clusters column at end depedning on clustering choices)
          if (objects_clustered == T & grayscale == F) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, subset(pointdata, obj_id == uniq_object)[clustering_variable]]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, subset(pointdata, obj_id == uniq_object)[clustering_variable]]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, subset(pointdata, obj_id == uniq_object)[clustering_variable]]
          }
          # produce color masks of cells/objects by anatomical region
          else if (anatomy_clustered == T & objects_clustered == F) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, subset(pointdata, obj_id == uniq_object)[anatomy_variable]]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, subset(pointdata, obj_id == uniq_object)[anatomy_variable]]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, subset(pointdata, obj_id == uniq_object)[anatomy_variable]]
          }
          # color objects based on unique RGB_ids calculated earlier
          else if (spatial_id_assign == T & grayscale == F) {
            obj_mask_r[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$R
            obj_mask_g[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$G
            obj_mask_b[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$B
          }
          # color objects based on unique spatial_ids calculated earlier - grayscale
          else if (spatial_id_assign == T & grayscale == T) {
            obj_mask_gray[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$spatial_id
          }
          # color objects based on unique cluster ids calculated earlier
          else if (objects_clustered == T & grayscale == T) {
            obj_mask_gray[obj_mask == uniq_object] = subset(pointdata, obj_id == uniq_object)$de_novo_regions
          }
          # if no clustering, then use default palette numbers to color objects
          else if ((objects_clustered == F) & (spatial_id_assign == F) & (anatomy_clustered == F)) {
            obj_mask_r[obj_mask == uniq_object] = col2rgb(palette)[1, obj_index]
            obj_mask_g[obj_mask == uniq_object] = col2rgb(palette)[2, obj_index]
            obj_mask_b[obj_mask == uniq_object] = col2rgb(palette)[3, obj_index]
          }
        }
        
        ##### save images ####
        if (grayscale == F) {
          # convert palette numbers to 8bit values 
          obj_mask_r = obj_mask_r / 255
          obj_mask_g = obj_mask_g / 255
          obj_mask_b = obj_mask_b / 255
        }
        
        # set up naming scheme
        type_of_image <- ''
        
        if (spatial_id_assign == T) {type_of_image <- 'spacing_unique_ids'}
        else if (objects_clustered == T & anatomy_clustered == F) {type_of_image <- 'clustered'}
        else if (objects_clustered == T & anatomy_clustered == T) {type_of_image <- 'anatomy_and_clustered'}
        else if (objects_clustered == F & anatomy_clustered == T) {type_of_image <- 'anatomy_labels'}
        else {type_of_image <- 'basic'}
        if (grayscale == T) {type_of_image <- paste0(type_of_image, '_gray')}
        
        # convert arrays to RGB image if needed
        if (grayscale == F) {
          img <- EBImage::transpose(rgbImage(Image(obj_mask_r), Image(obj_mask_g), Image(obj_mask_b)))
        }
        else if (grayscale == T) {
          img <- Image(obj_mask_gray)
        }
        
        # create actual image and save
        dir.create(paste0(data_folder, "obj_images"), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i]), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image, "/Point", p), showWarnings = FALSE)
        dir.create(paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image, "/Point", p, "/TIFs/"), showWarnings = FALSE)
        writeImage(img, paste0(data_folder, "obj_images", "/", save_name, "/", run_names[i], "/", type_of_image, "/Point", p, "/TIFs/", object_types[obj_index],".tif"))
      }
    }
    rm(list = c("img", "obj_mask", "obj_mask_r", "obj_mask_g", "obj_mask_b"))
  }

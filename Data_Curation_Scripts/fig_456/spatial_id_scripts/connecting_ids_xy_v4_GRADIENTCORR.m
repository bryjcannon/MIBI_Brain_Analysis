%input values
object_types = ["Tiled_cells"];
RGB_objects_img_path = "/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_MedADCase/NoAuBGFFtDenoised_New/obj_images_GRADIENTCORR/dc1/MedAD/unique_ids/TILING";
output_folder = "/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_MedADCase/NoAuBGFFtDenoised_New/obj_images_GRADIENTCORR/dc1/MedAD/unique_ids/TILING";

for obj_i = 1:length(object_types)

    %read in images
    disp("Reading in images.")
    RGB_objects_img = imread(char(join([RGB_objects_img_path, '/', object_types(obj_i),".tif"], "")));

    %re-create single RGB_id for objects
    disp(['Processing RGB object image -> ', char(object_types(obj_i))])
    obj_RGB_map = (256^2)*double(RGB_objects_img(:,:,1))+ 256*double(RGB_objects_img(:,:,2))  + double(RGB_objects_img(:,:,3));
    obj_RGB_ids = unique(obj_RGB_map); % use to check obj collected == total objects in R dataframe
    obj_RGB_ids(1) = []; % remove 0 from id's

    % initialize empty arrays for each desired property with size of number of objects / cells
    
    spatial_id = zeros(length(obj_RGB_ids), 1);
    x_centroid = zeros(length(obj_RGB_ids), 1);
    y_centroid = zeros(length(obj_RGB_ids), 1);
    Area = zeros(length(obj_RGB_ids), 1);
    MajorAxisLength = zeros(length(obj_RGB_ids), 1);
    MinorAxisLength = zeros(length(obj_RGB_ids), 1);
    Eccentricity = zeros(length(obj_RGB_ids), 1);
    Orientation = zeros(length(obj_RGB_ids), 1);
    Circularity = zeros(length(obj_RGB_ids), 1);
    Perimeter = zeros(length(obj_RGB_ids), 1);
    
    % for each id, make a mask of only that object, regionprops that, extract
    % properties, add to separate matrix - done to avoid issues with
    % bordering objects with no separation (don't want to add in
    % watershedding atm(.
    for id = 1:length(obj_RGB_ids)
        disp(obj_RGB_ids(id))
        temp_matrix = (obj_RGB_map == obj_RGB_ids(id));
        id_properties = regionprops(temp_matrix, 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation', 'Circularity', 'Perimeter');
        % add object properties to greater arrays
        spatial_id(id) = obj_RGB_ids(id);
        x_centroid(id) = id_properties.Centroid(1);
        y_centroid(id) = id_properties.Centroid(2);
        Area(id) = id_properties.Area;
        MajorAxisLength(id) = id_properties.MajorAxisLength;
        MinorAxisLength(id) = id_properties.MinorAxisLength;
        Eccentricity(id) = id_properties.Eccentricity;
        Orientation(id) = id_properties.Orientation;
        Circularity(id) = id_properties.Circularity;
        Perimeter(id) = id_properties.Perimeter;
    end
    
    %save struct as mat and load into RStudio, add to master_obj_data (use RGB_id as
    %connector)
    % save attributes and locations (mask) of the objects extracted from ez_segmenter in point folders :)
    disp("Saving .mat")
    imwrite(uint16(obj_RGB_map), char(join([output_folder, '/Point', num2str(point), '/TIFs/', object_types(obj_i),'_unique_mask.tiff'], "")));
    save(char(join([output_folder,'/',object_types(obj_i),'_tiled_properties.mat'], "")),'spatial_id', 'x_centroid', 'y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation', 'Circularity', 'Perimeter');
    
end
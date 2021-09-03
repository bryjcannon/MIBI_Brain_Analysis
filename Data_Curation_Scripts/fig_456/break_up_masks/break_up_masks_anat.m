%% break_up_masks.m
%  Script for breaking up masks by their color, saving them separately.
clear()
%input values
run_types = ["HiAD","MedAD","Ctrl"];
%if working with images
common_name_img = "_anat_id_exp_voronoi_CORRECTED.tif";
img_path = "/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/data_outputs/mask_expansion_pixel_capture/masks_anat/";
color_map = containers.Map([11546720,15656581,16629615,13284054,9498256], ["DG","CA4","CA3","CA2","CA1"]);
%output folder
output_folder = "/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/data_outputs/mask_expansion_pixel_capture/masks_anat/separated_masks";
if ~exist(output_folder, 'dir')
   mkdir(output_folder)
end

for run = 1:length(run_types)
    
    %read in images and calculate size
    disp("Reading in image.")
    start_mask = imread(char(join([img_path, '/', run_types(run), common_name_img], "")));
    sm_size = size(start_mask);

    %re-create single RGB_id for objects
    disp(['Processing RGB object image -> ', char(run_types(run))])
    obj_RGB_map = (256^2)*double(start_mask(:,:,1)) + 256*double(start_mask(:,:,2)) + double(start_mask(:,:,3));

    %parse out unique colors (representative of regions or other identity)
    unique_colors = unique(obj_RGB_map);
    unique_colors = unique_colors(2:end);
    
    mkdir(char(join([output_folder, '/', run_types(run)], "")));
    for color = 1:length(unique_colors)
        if isKey(color_map, unique_colors(color))
            new_mask = zeros(sm_size(1), sm_size(2));
            new_mask(obj_RGB_map == unique_colors(color)) = 1;
            imwrite((uint8(new_mask)), char(join([output_folder, '/', run_types(run), '/', color_map(unique_colors(color)), '_submask.tif'], "")));
        end
    end
end

        
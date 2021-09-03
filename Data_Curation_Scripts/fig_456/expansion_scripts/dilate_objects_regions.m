%% expand_object_regions
% Script for creating masks to expand coverage of de novo or other region types
% attached to either objects maps or object centroids to grab more of
% surrounding non-object space. Currently done either by object dilation or
% by vornoi diagram (work-in-progress).

%input values
run_types = ["HiAD", "MedAD", "Ctrl"];
method = "dilate";
%if working with images
common_name_img = "_all_denovoregions_FORCARPET.tif";
img_path = "/Volumes/BryJC_Stanford/paper1_analysis/Fig6/plots/fig_plots_final/overlays/set2/expanded";
%if working with coordinates
common_name_csv = "_forVoronoi.csv";
csv_path = "/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi";
%output folder
output_folder = "/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi";

for run = 1:length(run_types)
    
    if method == "dilate"
        %read in images
        disp("Reading in image.")
        start_mask = imread(char(join([img_path, '/', run_types(run), common_name_img], "")));

        %re-create single RGB_id for objects
        disp(['Processing RGB object image -> ', char(run_types(run))])
        obj_RGB_map = (256^2)*double(start_mask(:,:,1)) + 256*double(start_mask(:,:,2)) + double(start_mask(:,:,3));

        % Expansion Step
        %dilate out objects by n pixels
        dil_obj_map = imdilate(obj_RGB_map, strel('disk',25,4));
        imagesc(dil_obj_map)
        %
    elseif method == "voronoi"
        start_coord = readmatrix(char(join([csv_path, '/', run_types(run), common_name_csv], "")));
        [vc_x, vc_y] = voronoi(start_coord(:,1), start_coord(:,2));
    end
    
    %convert map back to RGB
    tiled_size = size(obj_RGB_map);
    R = zeros(tiled_size(1), tiled_size(2));
    G = zeros(tiled_size(1), tiled_size(2));
    B = zeros(tiled_size(1), tiled_size(2));
  
    B = mod(dil_obj_map, 256);
    G = mod(((dil_obj_map-B) / 256), 256);
    R = ((dil_obj_map-B) / 256^2) - (G/256);
    
    %create RGB array and save as tif
    end_mask = cat(3,R,G,B);
    imwrite((uint8(end_mask)), char(join([output_folder, "/", run_types(run), '_denovoregions_exp.tif'], "")));
end
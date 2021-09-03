%% expand_object_regions.m
%   Script for creating masks to expand coverage of de novo or other region types
%   attached to either objects maps or object centroids to grab more of
%   surrounding non-object space. Currently done either by object dilation or
%   by vornoi diagram (work-in-progress).

% input values
run_types = ["HiAD","MedAD", "Ctrl"];
method = "voronoi";
subset = 1;
subset_term = 'tangles';
% if working with images
common_name_img = "_all_denovoregions_FORCARPET.tif";
img_path = "/Volumes/BryJC_Stanford/paper1_analysis/Fig6/plots/fig_plots_final/overlays/set2/expanded";
% if working with coordinates
common_name_csv = "_forVoronoi_expanded.csv";
csv_path = "/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi";
color_map = containers.Map([1,2,3,4,5], [255,13487360,65280,10824234,16745466]);
% output folder
output_folder = "/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi_tangles";
if ~exist(output_folder, 'dir')
       mkdir(output_folder)
end

for run = 1:length(run_types)
    
    % read in images and calculate size
        disp("Reading in image.")
        start_mask = imread(char(join([img_path, '/', run_types(run), common_name_img], "")));
        sm_size = size(start_mask);
    
    if method == "dilate"

        % re-create single RGB_id for objects
        disp(['Processing RGB object image -> ', char(run_types(run))])
        obj_RGB_map = (256^2)*double(start_mask(:,:,1)) + 256*double(start_mask(:,:,2)) + double(start_mask(:,:,3));

        % Expansion Step
        % dilate out objects by n pixels
        new_map = imdilate(obj_RGB_map, strel('disk',10,4));
        imagesc(new_map)
        
    elseif method == "voronoi"
        %
        new_map = zeros(sm_size(1), sm_size(2));
        start_coord = readmatrix(char(join([csv_path, '/', run_types(run), common_name_csv], "")));
        sc_size = size(start_coord);
        
        dt = delaunayTriangulation(start_coord(:,1:2));
        [V,R] = voronoiDiagram(dt);
                
        for i = 1:sc_size(1)
            % if subset is active
            if (subset == 1) & (start_coord(i,4) ~= subset_term)
                continue
            end
            A = V(R{i},:);
            B = A(any(~isinf(A),2),:); % omit points at infinity
            bw = poly2mask(B(:,1), B(:,2), sm_size(1), sm_size(2));
            new_map(bw == 1) = color_map(start_coord(i,3));
        end
        
    end
    
    % convert map back to RGB
    R = zeros(sm_size(1), sm_size(2));
    G = zeros(sm_size(1), sm_size(2));
    B = zeros(sm_size(1), sm_size(2));
  
    B = mod(new_map, 256);
    G = mod(((new_map-B) / 256), 256);
    R = ((new_map-B) / 256^2) - (G/256);
    
    % create RGB array and save as tif
    end_mask = cat(3,R,G,B);
    
    subset_str = "";
    if subset == TRUE
        subset_str = char(join(['_', subset_term,]));
    end
    
    if method == "dilate"
        filename = char(join([output_folder, "/", run_types(run), subset_str, '_denovoregions_exp_dilate.tif'], ""));
    elseif method == "voronoi"
        filename = char(join([output_folder, "/", run_types(run), subset_str, '_denovoregions_exp_voronoi.tif'], ""));
    end
    % write image file
    imwrite((uint8(end_mask)), filename);
end

%% extra code previously used
% f = figure()
% hold on
% B=A(any(~isinf(A),2),:); % omit points at infinity
% fill(B(:,1), B(:,2), 'g');
% fill(polyshape(B), 'FaceColor', 'green');
% 
% %[vc_x, vc_y] = voronoi(start_coord(:,1), start_coord(:,2));
% 
% bw = poly2mask(x,y,256,256);
% Display the mask, drawing a line around the polygon.
% 
% imshow(bw)
% hold on
% plot(x,y,'b','LineWidth',2)
% hold off
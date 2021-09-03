%% expand_object_regions.m
%   Script for creating masks to expand coverage of de novo or other region types
%   attached to either objects maps or object centroids to grab more of
%   surrounding non-object space. Currently done either by object dilation or
%   by vornoi diagram (work-in-progress).
clear()

% input values
run_types = ["MedAD", "Ctrl", "HiAD"],
color_name = 'anat_id';
method = "voronoi";
subset = 0;
subset_name = 'FlowSOM_ids';
subset_term = 'plaques';
% if working with images
common_name_img = "_DGtoCA1_anat_FORCARPET.tif";
img_path = "/Volumes/BryJC_Stanford/paper1_analysis/Fig6/plots/fig_plots_final/overlays/set2/expanded";
% if working with coordinates
common_name_csv = "_anat_forVoronoi_expanded.csv";
csv_path = "/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi";
color_map = containers.Map(["DG","CA4","CA3","CA2","CA1"], [11546720,15656581,16629615,13284054,9498256]);
% output folder
output_folder = "/Volumes/BryJC_Stanford/other/ForDmitry/fig6_carpet_voronoi";
if ~exist(output_folder, 'dir')
       mkdir(output_folder)
end

for run = 1:length(run_types)
    
    % read in images and calculate size
        disp(join(["Reading in image:", run_types(run)]))
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
        % pull in objects - their centroid coordinates and additional information
        new_map = zeros(sm_size(1), sm_size(2));
        start_coord = readtable(char(join([csv_path, '/', run_types(run), common_name_csv], "")));
        sc_size = size(start_coord);
        
        % Expansion Step
        % compute delaunay triangulation and form voronoi diagram, getting vertices for each polygon
        dt = delaunayTriangulation([start_coord.X, start_coord.Y]);
        [V,R] = voronoiDiagram(dt);
        
        disp(join(['Processing total polygons -> ', string(sc_size(1))]))
        for i = 1:sc_size(1)
            % if subset is active
            if subset == 1 && ~strcmp(start_coord.(genvarname(subset_name))(i), subset_term) %genvarname for getting table column name
                continue
            end
            if mod(i, 10) == 0
                disp(join([run_types(run), ':Processing polygon # ', string(i), ' out of ', sc_size(1)]))
            end
            A = V(R{i},:);
            B = A(any(~isinf(A),2),:); % omit points at infinity
            bw = poly2mask(B(:,1), B(:,2), sm_size(1), sm_size(2)); % creates mask based upon the polygon vertices
            % pull out the object's id for the variable you want to color by
            color_name_id = start_coord.(genvarname(color_name))(i);
            if ~isnumeric(color_name_id)
                color_name_id = string(color_name_id);
            end
            new_map(bw == 1) = color_map(color_name_id); % colors the polygon based upon the color_name, or variable you're using to color your images
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
    if subset == 1
        subset_str = char(join(['_', subset_term,]));
    end
    
    if method == "dilate"
        filename = char(join([output_folder, "/", run_types(run), subset_str, '_', color_name, '_exp_dilate.tif'], ""));
    elseif method == "voronoi"
        filename = char(join([output_folder, "/", run_types(run), subset_str, '_', color_name, '_exp_voronoi.tif'], ""));
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
%% pixel_mask_data :: Create and export pixel level marker data by previously assigned mask identity %%
%  Input :: data tifs, mask tifs
%  Output :: .csv that stores pixel level information by run, mask identity

clear()
% input values
run_types = ["HiAD","MedAD", "Ctrl"];
run_sizes = [8192,8192; 11264,7680; 11264,8704];
mask_types = ["DG","CA4","CA3","CA2","CA1"];
num_markers = 36;
% if working with images
common_mask_name = "_submask.tif";
mask_path = "/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/data_outputs/mask_expansion_pixel_capture/masks_anat/separated_masks";
data_path = "/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/data_outputs/mask_expansion_pixel_capture/pixel_data";
% output folder
output_folder = "/Volumes/BryJC_Stanford/Analysis/KV_BJC_MIBI_Brain_Paper_1/Fig_MedRes/data_outputs/mask_expansion_pixel_capture/pixel_data_anat/anat_only";
output_name = "Anat";
if ~exist(output_folder, 'dir')
   mkdir(output_folder)
end

% Set up array for storing pixel level info
max_pixel_num = 0; % will represent greatest possible number of rows (or pixels) within masks
runs_dim = size(run_sizes);
for row = 1:runs_dim(1) %calculate the entire number of pixels within the images based on input data
    pixel_num = run_sizes(row,1) * run_sizes(row,2);
    max_pixel_num = max_pixel_num + pixel_num;
end
% Create pixel array for storing all run, mask, and marker data per pixel (i.e. row)
pixel_info_array = zeros(max_pixel_num, num_markers + 2);
curr_pixel_num = 1; % starting point for array
run_data_names = ''; % empty variable for data tif name storage

% Iterate over samples
for run = 1:length(run_types)
    % Load in all data tifs
    run_data_names = dir(fullfile(char(join([data_path,'/',run_types(run)], "")), '/*.tif')); % names of marker tifs to load in
    data_array = zeros(run_sizes(run,1), run_sizes(run,2), length(run_data_names)); % initialize array to store all marker data
    for marker=1:length(run_data_names) %import and store all marker data
        curr_data = imread(char(join([data_path, '/', run_types(run), '/', run_data_names(marker).name], "")));
        data_array(:,:,marker) = curr_data;
    end
    % Reshape data from l x w x n dimensions to (l * w) x n :: forces each row to represent a pixel from each image
    flip_data_array = reshape(data_array, [], length(run_data_names));
    
    % Iterate over sub masks
    for mask = 1:length(mask_types)
        % Load in sub mask
        curr_mask = imread(char(join([mask_path, '/', run_types(run), '/', mask_types(mask), common_mask_name], "")));
        % For pixels in mask, enter the corresponding intensity value for each data tif pixel into an array for this mask. Each new column is a tif, each row is a mask.
            % ToDo this: reshape the mask, take all rows == 1 and grab the data from the same rows in the reshaped data tifs.
            % Essentially the original mask into an array with just mask positive pixels and their data. Transverse this array and add to pixel_info_array.
        
        % Reshape, find mask pixels, and add to array
        flip_mask = reshape(curr_mask, [], 1); % forces each row to represent a pixel from each image
        flip_all = [flip_mask, flip_data_array]; % combine mask and marker arrays
        flip_all_mask_pixels = flip_all(flip_all(:,1,:) == 1, :); % keep only pixels contained in the mask
        flip_all_mask_pixels(:,1,:) = mask; % relabel mask positive value to mask idenity values
        flip_all_mask_pixels = [repmat(run, length(flip_all_mask_pixels), 1), flip_all_mask_pixels]; % add a column with run value
        
        % Add the array to the greater table - use pre-allocated array to speed up processing, iterating through each new addition of pixels.
        stop_num = curr_pixel_num + length(flip_all_mask_pixels) - 1; % where you will add to
        pixel_info_array(curr_pixel_num:stop_num, :) = flip_all_mask_pixels; % add pixel info to greater array
        curr_pixel_num = curr_pixel_num + length(flip_all_mask_pixels); % reset current row in greater array
    end
end
% Clear unused rows from pre-allocated array
pixel_info_array = pixel_info_array(1:stop_num, :);
% Rename rows, columns in table if needed
marker_names = extractfield(run_data_names, 'name');
for name = 1:length(marker_names)
    temp_name = marker_names(name);
    temp_name = temp_name{1};
    marker_names{name} = temp_name(1:end-4); % get rid of .tif suffix
end
% Convert array to table, rename masks according to their original mask names
pixel_info_table = array2table(pixel_info_array, 'VariableNames', ['Run', 'Mask', marker_names]);
pixel_info_table.Mask = string(pixel_info_table.Mask);
for mask_num = 1:length(mask_types)
    to_name = (pixel_info_table.Mask(:) == string(mask_num));
    pixel_info_table.Mask(to_name) = mask_types(mask_num);
end
% Export table as csv for R use
writetable(pixel_info_table, char(join([output_folder, '/', output_name, '_pixel_info.csv'], "")), 'Delimiter',',')
    
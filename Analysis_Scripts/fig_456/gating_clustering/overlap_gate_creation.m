%% overlap_gate_creation :: Create gatings for spatially overlapping cells and objects in segmented, labeled images %%
%  Input :: tiled tifs of cells, objects (labeled by unique RGB spatial ids)
%  Output :: 1) container.map files as .mat, contains unique cell ids as keys and unique overlapping objects as values
%           2) .mat, contains unique cell ids and binary object gating indicators
%  NOTE: want to eventually make this so full mapping transfers over as a gate

%Denote type of run (gate for cells overlapped with acellular features
%OR acellualr features overlapped with cells)
   cell_overlap = false;
   acellular_overlap = true;

%% 1) Load images (need to be grayscale spatial ids, not RGB)

%input values
object_names = ["cells", "tau-pathies", "amyloid-pathies"];
img_paths = ["/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/obj_images/All/spacing_unique_ids/TILING",...
             "/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_MedADCase/obj_images/All/MedAD/spacing_unique_ids/TILING",...
             "/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_HiADCase/obj_images/All/HiAD/spacing_unique_ids/TILING"];
output_folder = "/Volumes/BryJC_Stanford/paper1_analysis/Fig6";
if ~exist(output_folder, 'dir')
   mkdir(output_folder)
end

if cell_overlap == true
    %initialize tau_overlapped_cells and amyloid_overlapped_cells mapping tables
    tau_overlapped_cells_map = table('Size',[1, 3],'VariableTypes',{'double','double','double'},'VariableNames',{'cell_id','tau_id','GroupCount'});
    amyloid_overlapped_cells_map = table('Size',[1, 3],'VariableTypes',{'double','double','double'},'VariableNames',{'cell_id','amyloid_id','GroupCount'});
    %initialize tau_overlapped_cells and amyloid_overlapped_cells gating matrices
    tau_overlapped_cells_gate = [1, 1];
    amyloid_overlapped_cells_gate = [1, 1];
end
if acellular_overlap == true
    %initialize cells_overlapped_tau and cells_overlapped_amyloid mapping tables
    cells_overlapped_tau_map = table('Size',[1, 3],'VariableTypes',{'double','double','double'},'VariableNames',{'tau_id','cell_id','GroupCount'});
    cells_overlapped_amyloid_map = table('Size',[1, 3],'VariableTypes',{'double','double','double'},'VariableNames',{'amyloid_id','cell_id','GroupCount'});
    %initialize cells_overlapped_tau and cells_overlapped_amyloid gating matrices
    cells_overlapped_tau_gate = [1, 1];
    cells_overlapped_amyloid_gate = [1, 1];
end


%% 2) read in image and then convert to unique graysale spatial id, concatenating images into single multi-layer matrix
for sample = 1:length(img_paths)
    %load in all objects
    for i = 1:length(object_names)
        RGB_objects_img = imread(char(join([img_paths(sample), '/', object_names(i),".tif"], "")));
        if i == 1 %initialize new matrix for new sample
            multilayer_mat = zeros(size(RGB_objects_img));
        end
        multilayer_mat(:,:,i) = (256^2)*double(RGB_objects_img(:,:,1)) + 256*double(RGB_objects_img(:,:,2)) + double(RGB_objects_img(:,:,3));
    end

    % 3) Reshape matrix from l x w x n dimensions to (l * w) x n :: forces each row to represent a pixel from each image

    pixel_mat = reshape(multilayer_mat, [], length(object_names));

    % 4) Create mat1 -> table + groupCounts

    %convert pixel rows into table with assigned names
    pixel_table_cell = array2table(pixel_mat,'VariableNames',{'cell_id','tau_id','amyloid_id'});

    if cell_overlap == true
        %group pixels with positive cell ids and tau ids
        tau_pixels = groupcounts(pixel_table_cell,{'cell_id','tau_id'});
        tau_pixels_POS = tau_pixels(tau_pixels.tau_id > 0 & tau_pixels.cell_id > 0, :);
        tau_overlapped_cells_map = [tau_overlapped_cells_map; tau_pixels_POS];

        %group pixels with positive cell ids and amyloid ids
        amyloid_pixels = groupcounts(pixel_table_cell,{'cell_id','amyloid_id'});
        amyloid_pixels_POS = amyloid_pixels(amyloid_pixels.amyloid_id > 0 & amyloid_pixels.cell_id > 0, :);
        amyloid_overlapped_cells_map = [amyloid_overlapped_cells_map; amyloid_pixels_POS];

        % 5) Create mat2 -> gating mat and save

        %isolate out unique cells associated with tau overlap and add to main tau gating matrix 
        unique_tau_cells = unique(tau_pixels_POS.cell_id);
        unique_tau_cells(:,2) = ones(length(unique_tau_cells), 1);
        tau_overlapped_cells_gate = [tau_overlapped_cells_gate; unique_tau_cells];

        %isolate out unique cells associated with amyloid overlap and add to main amyloid gating matrix 
        unique_amyloid_cells = unique(amyloid_pixels_POS.cell_id);
        unique_amyloid_cells(:,2) = ones(length(unique_amyloid_cells), 1);
        amyloid_overlapped_cells_gate = [amyloid_overlapped_cells_gate; unique_amyloid_cells];

    end

    if acellular_overlap == true
        %group pixels with positive cell ids and tau ids
        cell_tau_pixels = groupcounts(pixel_table_cell,{'tau_id','cell_id'});
        cell_tau_pixels_POS = cell_tau_pixels(cell_tau_pixels.tau_id > 0 & cell_tau_pixels.cell_id > 0, :);
        cells_overlapped_tau_map = [cells_overlapped_tau_map; cell_tau_pixels_POS];

        %group pixels with positive cell ids and amyloid ids
        cell_amyloid_pixels = groupcounts(pixel_table_cell,{'amyloid_id','cell_id'});
        cell_amyloid_pixels_POS = cell_amyloid_pixels(cell_amyloid_pixels.amyloid_id > 0 & cell_amyloid_pixels.cell_id > 0, :);
        cells_overlapped_amyloid_map = [cells_overlapped_amyloid_map; cell_amyloid_pixels_POS];

        % 5) Create mat2 -> gating mat and save

        %isolate out unique tau objects associated with cell overlap and add to main cell gating matrix 
        unique_cells_tau = unique(cell_tau_pixels_POS.tau_id);
        unique_cells_tau(:,2) = ones(length(unique_cells_tau), 1);
        cells_overlapped_tau_gate = [cells_overlapped_tau_gate; unique_cells_tau];

        %isolate out unique amyloid objects associated with cell overlap and add to main cell gating matrix 
        unique_cells_amyloid = unique(cell_amyloid_pixels_POS.amyloid_id);
        unique_cells_amyloid(:,2) = ones(length(unique_cells_amyloid), 1);
        cells_overlapped_amyloid_gate = [cells_overlapped_amyloid_gate; unique_cells_amyloid];

    end

end

% 6) Save map and gate for overlapped cells

if cell_overlap == true
    %save tau overlaps
    tau_overlapped_cells_map = tau_overlapped_cells_map(2:end,:);
    writetable(tau_overlapped_cells_map, char(join([output_folder, '/', 'tau_overlapped_cells_map.csv'], "")), 'Delimiter',',')
    tau_overlapped_cells_gate = tau_overlapped_cells_gate(2:end,:);
    writematrix(tau_overlapped_cells_gate, char(join([output_folder, '/', 'tau_overlapped_cells_gate.csv'], "")), 'Delimiter',',');

    %save amyloid overlaps

    amyloid_overlapped_cells_map = amyloid_overlapped_cells_map(2:end,:);
    writetable(amyloid_overlapped_cells_map, char(join([output_folder, '/', 'amyloid_overlapped_cells_map.csv'], "")),'Delimiter',',')
    amyloid_overlapped_cells_gate = amyloid_overlapped_cells_gate(2:end,:);
    writematrix(amyloid_overlapped_cells_gate, char(join([output_folder, '/', 'amyloid_overlapped_cells_gate.csv'], "")),'Delimiter',',')
end

if acellular_overlap == true
    %save tau overlaps
    cells_overlapped_tau_map = cells_overlapped_tau_map(2:end,:);
    writetable(cells_overlapped_tau_map, char(join([output_folder, '/', 'cells_overlapped_tau_map.csv'], "")), 'Delimiter',',')
    cells_overlapped_tau_gate = cells_overlapped_tau_gate(2:end,:);
    writematrix(cells_overlapped_tau_gate, char(join([output_folder, '/', 'cells_overlapped_tau_gate.csv'], "")), 'Delimiter',',');

    %save amyloid overlaps

    cells_overlapped_amyloid_map = cells_overlapped_amyloid_map(2:end,:);
    writetable(cells_overlapped_amyloid_map, char(join([output_folder, '/', 'cells_overlapped_amyloid_map.csv'], "")),'Delimiter',',')
    cells_overlapped_amyloid_gate = cells_overlapped_amyloid_gate(2:end,:);
    writematrix(cells_overlapped_amyloid_gate, char(join([output_folder, '/', 'cells_overlapped_amyloid_gate.csv'], "")),'Delimiter',',')
end

%% MIBIcreateObjectNeighborMatrix.m for EZSeg inputs
% Author: Erin McCaffrey
% Modified by: Dmitry & Kausalia
% Saves number of neighbours for every object in current point 
%% read in cell data
absolutePath = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellDataAnalysis/Pooled/';

% define path to where the plaque mask for each point is
% Path to Mask of Object/Cell type that is in the center of distance calculation
pathMask = [absolutePath,'masksAmyloid']; % 'masksTangles' 'masksAmyloid' 'masksEndothelia'

% Define mask name
%mask_name = 'vessel_CD31_CD105.tif'; % BBB1
%mask_name = 'vessel_MCT1.tif'; % BBB2 
%mask_name = 'tauTangles.tif'; % Tangles
mask_name = 'amyloidPlaques.tif'; % Plaques

% paths to each ezSeg or DeepCell masks folder. Path to Mask of Object/Cell type want
% to count amount present in the neighbourhood of the center object/cell
% type
paths_ez = {[absolutePath,'maskNeurons_Deepcell'],[absolutePath, 'masksAstrocyteProcess'],[absolutePath,'masksMicrogliaProcess'],[absolutePath,'masksEndothelia']};
     
names_ez = {'neurons.tif','astrocyte_process.tif','microglia_process.tif','vessel_CD31_CD105.tif',}; %, 'vessel_MCT1.tif'

resultsPath = [absolutePath,'/RealRuns/AmyloidPlaques_MatrixByIndividual_Objects'];

if ~exist(resultsPath, 'dir')
   mkdir(resultsPath)
end

points = 1:85; %cohort data to analyze

%% initiate empty matrices for cell neighborhood data
object_neighborhood_counts = [];
%object_neighborhood_freqs = [];

%% Initiate distance threshold in px and extract mask name
t = 100; % neighborhood distance in pixels

[~,name,~] = fileparts(mask_name);

%% Build data matrix for each point
ez_out_names = cell(length(names_ez), 1);

% iterate through all points
for i = 1:length(points)
    point = points(i);

    % load relevant data
    disp(['Working on point: ',num2str(point)]);

    % load data for mask
    if ~isfile([pathMask,'/Point',num2str(point),'/', mask_name])
        continue
    end
    target_mask = imread([pathMask,'/Point',num2str(point),'/', mask_name]);
    target_mask_bw = imbinarize(target_mask);
    target_boundaries = bwboundaries(target_mask_bw);
    
    % initialize results vector for the current point
    results = zeros(size(target_boundaries, 1), length(names_ez) + 2);
    
    for ez_id = 1:length(names_ez)
        % Save ez mask name for column output names
        [~, ez_out_name, ~] = fileparts(char(names_ez(ez_id)));
        ez_out_names(ez_id, 1) = {ez_out_name};
        
        % get the location of origin and target objects
        origin_mask = imread([char(paths_ez(ez_id)),'/Point',num2str(point),'/', char(names_ez(ez_id))]);
        origin_mask_bw = imbinarize(origin_mask);
        origin_boundaries = bwboundaries(origin_mask_bw);
       
        % append count results
        results(1:size(target_boundaries, 1),1) = point;

        % iterate through all objects 
        for j=1:size(target_boundaries, 1)
            distanceMat = NaN(length(origin_boundaries),1);
            
            for c=1:length(origin_boundaries)
                % get origin x and y coordinates
                origin_coords = cell2mat(origin_boundaries(c));

                % get target x and y coordinates
                target_coords = cell2mat(target_boundaries(j));
                
                p1.x = origin_coords(:,1);
                p1.y = origin_coords(:,2);
                
                p2.x = target_coords(:,1);
                p2.y = target_coords(:,2);

                distanceMat(c) = min_dist_between_two_polygons(p1, p2);
            end

            % get the indices of distances less than thresh
            neighbors = find(distanceMat <= t);
                        
            results(j,2) = j;  % Saves number of neighbours for every object in current point 
            results(j,(ez_id + 2)) = length(neighbors);
        end
        % calculates ave number of objects of a specific type for each plaque
        %results(1,(ez_id + 2)) = results(1,(ez_id + 2)) / size(target_boundaries, 1);
    end 
    
    % Save current point results
    object_neighborhood_counts = [object_neighborhood_counts; results];
    
    % produce frequency results
    %results_freqs = results;
    %results_freqs(1,3:size(results,2)) = results_freqs(1,3:size(results,2))/length(names_ez);
    %object_neighborhood_freqs = [object_neighborhood_freqs; results_freqs];
end
        
%% export as csv
channelLabels = ['Point_num';'Object_num'; ez_out_names ];
TEXT.PnS = channelLabels;
%csvwrite_with_headers([resultsPath,'/', name, '_TangleThread_counts_', num2str(t), 'px.csv'],object_neighborhood_counts,TEXT.PnS)
csvwrite_with_headers([resultsPath,'/', name, '_AmyloidPlaque_counts_', num2str(t), 'px.csv'],object_neighborhood_counts,TEXT.PnS)

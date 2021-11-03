%% Run First 
% !!!Pipeline for nuclear segmentation using pixel probabilities from deepCell

% Changed the pipeline to use the boundaries from the deepcell as boudaries
% instead of the raw nuclear intensities.

% base_dir = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSegData_uci2712J/190505HiResDG/'
% DGBand = [5,6,7,9,25,26,27,28,43,44];%DG band of DGrun

base_dir = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSegData_uci2712J/190604HiResCA2/'
DGBand =[5,9,12,23,26,30];%DG band of CA2run

%base_dir = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSegData_uci2712J/DGCA2Pooled/'



datasize=1024;

% change these
massDS = MibiReadMassData([base_dir, 'info/190503_TMAADPanel.csv']);

resultsPath = [base_dir, 'single_cell_dynamic_expansion_Smaller_no_fftnoiseDist'];
deepPath = [base_dir, 'segmentation_masks/'];
figuresPath = [base_dir, 'expansion_figures_Smaller_no_fftnoiseDist/'];
%figuresPath = [base_dir, 'expansion_figures_CA2_Smaller_no_fftnoise/'];

if ~exist(figuresPath, 'dir')
    mkdir(figuresPath);
end

% used for probability map
% this value can be changed. if want to make nuclei segmentation mask smaller 
prob_threshold = 0;
prominance_threshold = 0.2;% Used in watershedding increase values to get more cell separation. increases the interval of highest values. at 0.2 taking top 20% of intensity values
labelMaskThresh = 0.7;% if increase will set more object to BG
t=40;% cells smaller than t will get merged into closets neighbouring cell of larger cell border
objectSizeFilter = 100; % removes objects smaller than 100 pixels in area
dimObjectThresh = 0.05;% removes all objects that has a prob <= to 0.05 to be an on object. removes dim objects

S = dir([base_dir, 'no_fftnoise/']);
NumPoints = sum([S(~ismember({S.name},{'.','..'})).isdir]);



for pointNumber = 1:NumPoints %1:NumPoints [5,6,7,9,25,26,27,28,43,44] DG band of DGrun, [5,9,12,23,26,30]DG band of CA2run
    
    %% Get perimiters of nuclei for plotting
    % load data and get nuclear markers
    
    load([base_dir,'no_fftnoise/Point',num2str(pointNumber),'/dataNoFFTNoise.mat']);
    imSize = (size(countsNoNoise,1));
    % sum nuclear markers to increase contrast
    nucleiChannels = {'HistoneH3Lyo'};
    [tf loc] = ismember(nucleiChannels,massDS.Label);
    nucIm = sum(countsNoNoise(:,:,loc),3);
    maxv=25;
    rgb_image = MibiGetRGBimageFromMat(nucIm,maxv);

    % %% Get maxima from deep learning probabilities.less smooth mask
    probNuc = double(imread([deepPath, 'Point', num2str(pointNumber), 'nuc_interior_less_smoothed.tiff']));
    %probNuc = double(imread([deepPath, tif_names(tif_num).name, '/feature_1_frame_0.tif']));

    % find local maxima in probability map
    probNucCap = probNuc;
    figure;
    imagesc(probNucCap);
    probNucCap(probNucCap < prob_threshold) = 0;
    maxs = imextendedmax(probNucCap, prominance_threshold);
    rgb_image_perim_extMax = imoverlay(rgb_image , maxs, [1 0 0]);   
    figure;
    imagesc(rgb_image_perim_extMax);
    
    if ~exist([figuresPath, 'Point', num2str(pointNumber)], 'dir')
        mkdir([figuresPath, 'Point', num2str(pointNumber)],'/TIFs/');
    end
    imwrite(rgb_image_perim_extMax, [figuresPath, 'Point', num2str(pointNumber),'/TIFs/','nuc_probs.tif']);
    
    %% watershed over the deep results
    bw1 = zeros(size(probNuc));
    bw1(probNuc>dimObjectThresh) = 1;
    bw2 = bwareaopen(bw1,objectSizeFilter);
   
    %this value can be changed
    %for loop through each cell object in mask, get radius with region
    %props (or diameter)
    %and then custom define the SE for each object to get more custom
    %expansion (use regionprops on bw2 to get unique object values)    
    nuc_stats = regionprops('table',bw2,'Area','Circularity','Eccentricity','MajorAxisLength','MinorAxisLength');
    [B,L] = bwboundaries(bw2,'noholes');
    
    bw = zeros(size(bw2));
    for cell = 1:length(B)
        % If cell is neuron. Was divided by 4 for DGCA2 region
        if nuc_stats.Area(cell) > 600 
            if ismember (pointNumber,DGBand)
                cell_radius = ceil((nuc_stats.MajorAxisLength(cell) + nuc_stats.MinorAxisLength(cell)) / 40);%for DG band   
            else
                cell_radius = ceil((nuc_stats.MajorAxisLength(cell) + nuc_stats.MinorAxisLength(cell)) / 8);
            end
       
            % If cell is  ellispoid (microglia).Was divided by nothing for DGCA2 region  
        elseif nuc_stats.Eccentricity(cell) > 0.7 
            %cell_radius = ceil(nuc_stats.MinorAxisLength(cell) / 2);
            cell_radius = ceil(nuc_stats.MinorAxisLength(cell) / 2);
        
            % If cell is astrocyte.Was divided by 2 for DGCA2 region
        else 
          %cell_radius = ceil((nuc_stats.MajorAxisLength(cell) + nuc_stats.MinorAxisLength(cell)) / 4);  
          cell_radius = ceil((nuc_stats.MajorAxisLength(cell) + nuc_stats.MinorAxisLength(cell)) / 2);
        end
        SE = strel('disk',cell_radius); % calculate a disk filter
        
        bw_temp = L;
        bw_temp(bw_temp ~= cell) = 0;
        bw_temp(bw_temp == cell) = 1;
        bw = bw + imdilate(bw_temp, SE);
    end
    
    maxsFix = bw & maxs;

    % The WATERSHEDDING PART.modify the image so that the background pixels and the extended maxima pixels are forced to be the only local minima in the image.
    Jc = imcomplement(probNuc);
    I_mod = imimposemin(Jc, ~bw | maxsFix);
    L = watershed(I_mod);
    labeledImage = label2rgb(L);
    
    %labels BG as 0 and all cells as 1 = mask
    cellPerimNewMod= L;
    cellPerimNewMod(L>0) = 100;% used to be 100, set to 10000 avoid erase cells ID higher than 100
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    
    rgb_image_cellPerim = imoverlay(rgb_image, cellPerimNewMod, [1 0 0]);
    figure;
    imagesc(rgb_image_cellPerim);
    
    imwrite(rgb_image_cellPerim, [figuresPath, 'Point', num2str(pointNumber),'/TIFs/', 'cell_perim.tif']);

    %% 1. For each label, decide whether it is of a nucleus/ background
    
    %making empty matrices
    labelNum = length(unique(L(:)));
    labelIdentity = zeros (labelNum,1);
    labelPixelsPercentInNucleiMask = zeros(labelNum,1);
    labelSize = zeros(labelNum,1);

    for i=1:labelNum
        [r c] = find(L==i);
        labelSize(i) = length(r); % cal size in pixel
        labelMask = (L==i);
        labelPixelsNumInNucleiMask = sum(bw(labelMask));
        labelPixelsPercentInNucleiMask(i) = labelPixelsNumInNucleiMask / labelSize(i);
        if (labelPixelsPercentInNucleiMask(i) > labelMaskThresh)
            labelIdentity(i) = 1;
        end
    end

    % 2. Merge small regions within the nuclei mask with their neighbours
    
    keepVec = ones(labelNum,1);
    newL = L;
    for i=1:labelNum
        if (labelIdentity(i) == 1) && (labelSize(i) < t)%merge small cell to larger cell. 
            disp(['Removing label ',num2str(i),'. Size: ',num2str(labelSize(i))]);
            % get neighbour with largest border that is also in nuclear region
            [neighbourLabels , neighbouringRegionSize] = MibiGetNeighbourLabels (newL, i);
            found = 0;
            [neighbouringRegionSizeS , neighbouringRegionSizeSInd] = sort(neighbouringRegionSize,'descend');
            neighbourLabelsS = neighbourLabels(neighbouringRegionSizeSInd);
            maxInd = 1;
            while ~found
                mergeLabelId = neighbourLabelsS(maxInd);
                if (~(mergeLabelId == 0) && (labelIdentity(mergeLabelId) == 1))
                    found = 1;
                else
                    maxInd = maxInd+1;
                end
                if (maxInd >length(neighbourLabelsS)) % reached end of neighbours with no good merging candidate
                    disp (['Warning: no good merging target found for label', num2str(i), '. Keeping it.']);
                    break;
                end
            end
            % update
            if (maxInd <= length(neighbourLabelsS))
                [newL] = MibiMergeLabels (newL, i, mergeLabelId);
                keepVec(i) = 0;
            end
        end
    end

    % Update label numbers to account for deleted labels
    allLabels = [1:labelNum];
    currLabels =  allLabels(keepVec == 1);
    labelIdentityNew = zeros(length(currLabels),1); 
    newLmod = newL;
    for i = 1:length(currLabels)
        newLmod(newLmod == currLabels(i)) = i;
        labelIdentityNew(i) = labelIdentity(currLabels(i));
    end

    cellPerimNewMod= bwperim(newLmod);
    cellPerimNewMod= newLmod;
    cellPerimNewMod(newL>0) = 100;
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    rgb_image_cellPerimNewMod = imoverlay(rgb_image, cellPerimNewMod, [1 0 0]);
    figure;
    imagesc(rgb_image_cellPerimNewMod); 
    
    output_path = [resultsPath, '/Point',num2str(pointNumber)];
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    save([output_path, '/segmentationParams.mat'],'newLmod','cellPerimNewMod','labelIdentityNew');
    writematrix(newLmod,[output_path, '/newLmod.csv']);
    imwrite((uint16(newLmod)),[output_path, '/newLmod.tiff']);    
    close all;
end

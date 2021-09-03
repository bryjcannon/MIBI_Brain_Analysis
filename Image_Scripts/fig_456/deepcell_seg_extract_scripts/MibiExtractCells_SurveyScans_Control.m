
%% RUN SECOND!!!! Pipeline for nuclear segmentation using pixel probabilities from deepCell

% MibiExtractSingleCellDataFromSegmentation
% assign expression data to a cell based on segmentation
% I use regionprops to speed up. Gives same results as version2, but much
% faster!

datasize=512;

% set path to (Set1), (Set2), and (Set3) to get data from all three parts 
 path = '/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/NoAuBGFFtDenoised_New_(Set2)/';
% path = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/HiResScan_TMA/DeepCellSeg_HiResTMA/DeepcellAnalysis_new/AD_Panel/Feb2019/deepcellAD1/';
% path = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/HiResScan_TMA/DeepCellSeg_HiResTMA/DeepcellAnalysis_new/AD_Panel/Feb2019/deepcellAD2/';
% path = '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/TMA/HiResScan_TMA/DeepCellSeg_HiResTMA/DeepcellAnalysis_new/AD_Panel/Nov2018/deepcellAD3/';
 
 massDS = MibiReadMassData([path, 'info/190503_TMAADPanel.csv']);
% massDS = MibiReadMassData([path, 'info/181101_AD_Orient.csv']);
% massDS = MibiReadMassData([path, 'info/190107_Common_Panel.csv']);
% massDS = MibiReadMassData([path, 'info/181101_AD_Orient.csv']);
 
pathSegment = [path, 'segmentationDC/single_cell_dynamic_expansion'];
resultsDir = [path,'segmentationDC/fcs_cells'];
mkdir(resultsDir);
clusterChannels = {'HistoneH3Lyo'};

S = dir(pathSegment);
NumPoints = sum([S(~ismember({S.name},{'.','..'})).isdir]);
pointList = [1:105];

% For deepcellPD p = 1:NumPoints. For deepcellAD1 p = 1:NumPoints. For deepcell AD2 p = 12:22. For deepcellAD3 p = [3:9,11].
for p = 1:NumPoints
    disp(['point',num2str(p)]);
    pointNumber=p
    
    % load data and get nuclear markers
    if isfile([path,'no_fftnoise_3X/Point',num2str(pointNumber),'/dataNoFFTNoise.mat'])
       load([path,'no_fftnoise_3X/Point',num2str(pointNumber),'/dataNoFFTNoise.mat']);
       channelnames = massDS.Label;
       [~, clusterChannelsInds] = ismember(clusterChannels,massDS.Label);
       channelNum = length(massDS);
    else       
        tifsFolder= dir(fullfile([path,'no_fftnoise_3X/Point',num2str(pointNumber),'/TIFs'],'*.tif'));
        countsNoNoise=zeros(datasize,datasize,numel(tifsFolder));
        for k = 1:numel(tifsFolder)
            if tifsFolder(k).name == "totalIon.tif"
                continue;
            end
            tifpath = [tifsFolder(k).folder, filesep, tifsFolder(k).name];
            img = imread(tifpath);
            countsNoNoise(:,:,k) = img;            
        end
        channelnames=strsplit([tifsFolder.name],".tif");
        channelnames=channelnames(1:k);
        [~, clusterChannelsInds] = ismember(clusterChannels,channelnames);
        channelnames=channelnames';  
        channelNum=length(channelnames);
    end
    
    %load([path, '/no_noise', '/Point',num2str(pointNumber),'/dataDeNoise.mat']);
    load([pathSegment,'/Point',num2str(pointNumber),'/segmentationParams.mat']);
    labelNum = max(max(newLmod));
    writematrix(newLmod,[resultsDir,'/newLmod_p',num2str(pointNumber),'.csv']);
    
    stats = regionprops(newLmod,'Area','PixelIdxList');
    countsReshape= reshape(countsNoNoise,size(countsNoNoise,1)*size(countsNoNoise,2),channelNum);
    % make a data matrix the size of the number of labels x the number of markers
    data = zeros(labelNum,channelNum);
    dataScaleSize = zeros(labelNum,channelNum);
    cellSizes = zeros(labelNum,1);

    % for each label extract information
    for i=1:labelNum
        currData = countsReshape(stats(i).PixelIdxList,:);
        data(i,:) = sum(currData,1);
        dataScaleSize(i,:) = sum(currData,1) / stats(i).Area;
        cellSizes(i) = stats(i).Area;
    end

    % get the final information only for the labels with 
    % 1.positive nuclear identity (cells)
    % 2. That have enough information in the clustering channels to be
    % clustered
    labelIdentityNew2 = labelIdentityNew([1:end-1]); % fix bug resulting from previous script
    sumDataScaleSizeInClusterChannels = sum(dataScaleSize(:,clusterChannelsInds),2);
    %labelIdentityNew2(sumDataScaleSizeInClusterChannels<0.001) = 2;% 2 is not a cell, 1 is a cell
    dataCells = data(labelIdentityNew2==1,:);
    dataScaleSizeCells = dataScaleSize(labelIdentityNew2==1,:);
    labelVec=find(labelIdentityNew2==1);
    
    % get cell sizes only for cells
    cellSizesVec = cellSizes(labelIdentityNew2==1);

    % asinh transform
    dataScaleSizeCellsTrans = asinh(dataScaleSizeCells);
    dataCellsTrans = asinh(dataCells);
    
%     % standardize
%     dataScaleSizeCellsTransStd = zscore(dataScaleSizeCellsTrans);
%     dataCellsTransStd = zscore(dataCellsTrans);
%     dataCellsStd = zscore(dataCells);
%     dataScaleSizeCellsStd = zscore(dataScaleSizeCells);

    dataTransL = [labelVec,cellSizesVec,dataCellsTrans];
    dataScaleSizeTransL = [labelVec,cellSizesVec, dataScaleSizeCellsTrans];
    dataL = [labelVec,cellSizesVec, dataCells];
    dataScaleSizeL = [labelVec,cellSizesVec, dataScaleSizeCells];
%     dataTransStdL = [labelVec,dataCellsTransStd];
%     dataScaleSizeTransStdL = [labelVec,dataScaleSizeCellsTransStd];
%     dataStdL = [labelVec,dataCellsStd];
%     dataScaleSizeStdL = [labelVec,dataScaleSizeCellsStd];

    channelLabelsForFCS = ['cellLabelInImage';'cellSize';channelnames];

    %% save fcs   
    TEXT.PnS = channelLabelsForFCS;
    TEXT.PnN = channelLabelsForFCS;
    
    % handles points with 0 cells (just don't create an fcs file)
    if labelIdentityNew == 0
        continue
    end
    save([pathSegment,'/Point',num2str(pointNumber),'/cellData.mat'],'labelIdentityNew2','labelVec','cellSizesVec','dataCells','dataScaleSizeCells','dataScaleSizeCellsTrans','dataCellsTrans','channelLabelsForFCS');
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataFCS.fcs'],dataL,TEXT);
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataScaleSizeFCS.fcs'],dataScaleSizeL,TEXT);
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataTransFCS.fcs'],dataTransL,TEXT);
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataScaleSizeTransFCS.fcs'],dataScaleSizeTransL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataStdFCS.fcs'],dataStdL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataScaleSizeStdFCS.fcs'],dataScaleSizeStdL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataTransStdFCS.fcs'],dataTransStdL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataScaleSizeTransStdFCS.fcs'],dataScaleSizeTransStdL,TEXT);
    
    writeFCS([resultsDir,'/dataFCS_p',num2str(pointNumber),'.fcs'],dataL,TEXT);
    writeFCS([resultsDir,'/dataScaleSizeFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeL,TEXT);
    writeFCS([resultsDir,'/dataTransFCS_p',num2str(pointNumber),'.fcs'],dataTransL,TEXT);
    writeFCS([resultsDir,'/dataScaleSizeTransFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeTransL,TEXT);
%     writeFCS([resultsDir,'/dataStdFCS_p',num2str(pointNumber),'.fcs'],dataStdL,TEXT);
%     writeFCS([resultsDir,'/dataScaleSizeStdFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeStdL,TEXT);
%     writeFCS([resultsDir,'/dataTransStdFCS_p',num2str(pointNumber),'.fcs'],dataTransStdL,TEXT);
%     writeFCS([resultsDir,'/dataScaleSizeTransStdFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeTransStdL,TEXT);
    
end

% dlmwrite('Col_names.txt', channelLabelsForFCS)
% B = [A; dataScaleSizeL];
% dlmwrite('Scaled_Cell_Data.txt',dataScaleSizeL, '\t') 
% dlmwrite('Raw_Cell_Data.txt',dataL, '\t') 
% 
% fileID = fopen('col_names.txt','w');
% formatSpec = '%s\t%d\t%2.1f\t';
% [nrows,ncols] = size(channelLabelsForFCS);
% for row = 1:nrows
%     fprintf(fileID,formatSpec,channelLabelsForFCS{row,:});
% end
% fclose(fileID);



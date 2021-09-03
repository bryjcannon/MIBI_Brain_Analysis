% MIBI smooth stitching script
% Author: Dmitry Tebaykin
% Contact: dmitry.tebaykin@stanford.edu

%% Necessary parameters
% Provide the XML file path from the run, stitching parameters will be
% extracted automatically. Example: xmlFileName = '180501_Hip_Panel3ug_Final-2b.xml';
xmlFileName = ''; % Leave this blank for now

% Output folder for stitched images. Default: current working directory
OutputFolder = ['/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/','StitchedDeepCellSegGated'];

% Point to folder where the Point folders with TIFs are located. pwd stands for current
% working directory (path at the top of Matlab)
% Example for a different folder: TIFs_PATH = [pwd, '/extracted']
% This relies on your TIFs being inside TIFs_PATH/PointX/TIFs/ folder. To
% change this - modify MibiStitchLoopSupport around line 62
TIFs_PATH = {'/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/190604HiResCA2/expansion_figures_CA2_Gated', ...
    '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/DeepCellSeg_uci2712J/190505HiResDG/expansion_figures_DG_Gated'}; 
    
%TIFs_PATH = {'/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/190604_HiResADCA2/uci2717JCA2Final/corrected', ...
%    '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/May_June2019/HiResScans_uci2717J/190505HiResADuci2717J/190505HiResFinal1/corrected'}; 

dataSize = 1020; % Resolution of one frame, minus 4 pixels. Example: 512x512 - 4 = 508, 1024x1024-4 = 1020.

% Stitch start and end points.
startPoint = 1; % Start stitching from this point number
endPoint = 85; % End stitching with this point number. Zero means all points to the end
skipPoints = []; % These points will be skipped during stitching (The stitch should advance, leaving a blank space)

%% Set these if no XML is available
% Set stitching parameters manually if no XML file is available
xNumPoint = 5; % Set this to the number of rows
yNumPoint = 17; % Set this to the number of columns
direction = 1; % Choose starting stitch direction: 0 = left, 1 = right

% Gather all channel names
%ImageNamesP = dir([TIFs_PATH{1},'/Point', num2str(startPoint),'/TIFs/*tif']); 
ImageNamesP = dir([TIFs_PATH{1},'/Point', num2str(startPoint),'/*tif']); 
ImageNamesP = {ImageNamesP.name};
ImageNames{1, length(ImageNamesP)} = [];
for i = 1:length(ImageNamesP)
    ImageNames{1,i} = ImageNamesP{i};
end
allChannels = ImageNames;

% Overwrite the above to stitch specific channels only
%allChannels = ["cell_perim.tif"];
allChannels = ["overlay5.tif"];

% Create a weights matrix
weights = zeros(dataSize, dataSize);
for k = 1:floor(dataSize/2)
    weights([k, dataSize - k + 1], k : dataSize - k + 1) = k;
    weights(k : dataSize - k + 1, [k, dataSize - k + 1] ) = k;
end

% Start the loop Change MibiStitchLoopSupportJoinRun or
% MibiStitchLoopSupportJoinRunTextV2
for channelNonChar = allChannels
    channel = char(channelNonChar);
    disp(['Making stitched TIF for: ', channel]);
    MibiStitchLoopSupportJoinRunPaperHiRes2717JSeg;
end
disp('Finished stitching all channels');
    


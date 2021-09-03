% MIBI smooth stitching script
% Author: Dmitry Tebaykin
% Contact: dmitry.tebaykin@stanford.edu

%% Necessary parameters
% Provide the XML file path from the run, stitching parameters will be
% extracted automatically. Example: xmlFileName = '180501_Hip_Panel3ug_Final-2b.xml';
xmlFileName = ''; % Leave this blank for now

% Output folder for stitched images. Default: current working directory
OutputFolder = ['/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/obj_images/DGtoCA1_Neurons_MetaDisease/anatomy_and_clustered/TILING/'];
% 1 for grayscale, 3 for RGB image
img_color = 3;
% Point to folder where the Point folders with TIFs are located. pwd stands for current
% working directory (path at the top of Matlab)
% Example for a different folder: TIFs_PATH = [pwd, '/extracted']
% This relies on your TIFs being inside TIFs_PATH/PointX/TIFs/ folder. To
% change this - modify MibiStitchLoopSupport around line 62
% TIFs_PATH = {'/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/August2019/MedResScan_uci4516HIPR/Ctrl_EZSEG/CombinedTIFs'};
TIFs_PATH = {'/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/obj_images/DGtoCA1_Neurons_MetaDisease/Ctrl1/anatomy_and_clustered/',...
'/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/obj_images/DGtoCA1_Neurons_MetaDisease/Ctrl2/anatomy_and_clustered/',...
'/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_ControlCase/obj_images/DGtoCA1_Neurons_MetaDisease/Ctrl3/anatomy_and_clustered/'};

% TIFs_PATH = {'//Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/August2019/MedResScan_uci4516HIPR/190828_uci4516HIPR_FinalSet1/NoAuBGFFtDenoised_New/no_fftnoise_3X', ...
%     '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/August2019/MedResScan_uci4516HIPR/190828_uci4516HIPR_FinalSet2/NoAuBGFFtDenoised_New/no_fftnoise_3X', ... 
%     '/Volumes/KausaliaHD/AllProjects/BRAIN DATA/MIBI/MIBIData Matlab/August2019/MedResScan_uci4516HIPR/190828_uci4516HIPR_FinalSet3n/NoAuBGFFtDenoised_New/no_fftnoise_3X'}; 

dataSize = 512; % Resolution of one frame.

% Stitch start and end points.
startPoint = 1; % Start stitching from this point number
endPoint = 275; % End stitching with this point number. Zero means all points to the end
skipPoints = []; % These points will be skipped during stitching (The stitch should advance, leaving a blank space)

%% Set these if no XML is available
% Set stitching parameters manually if no XML file is available
xNumPoint = 20; % Set this to the number of rows
yNumPoint = 15; % Set this to the number of columns
direction = 1; % Choose starting stitch direction: 0 = left, 1 = right

% Gather all channel names
ImageNamesP = dir([TIFs_PATH{1},'/Point', num2str(startPoint), '/TIFs/*tif']);
ImageNamesP = {ImageNamesP.name};
% used to remove any hidden files (doubles of tifs)
% ImageNamesP = ImageNamesP((length(ImageNamesP)/2)+1:length(ImageNamesP));
ImageNames{1, length(ImageNamesP)} = [];
for i = 1:length(ImageNamesP)
    ImageNames{1,i} = ImageNamesP{i};
end
allChannels = ImageNames;

% Overwrite the above to stitch specific channels only
%allChannels = ["CD56Lyo.tif"];

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
    %MibiStitchLoopSupportJoinRunPaperControlEZSeg;
    MibiStitchLoopSupportText_Ctrl_noBorders;
end
disp('Finished stitching all channels');
    


% MIBI smooth stitching script
% Author: Dmitry Tebaykin
% Contact: dmitry.tebaykin@stanford.edu

%% Necessary parameters
% Provide the XML file path from the run, stitching parameters will be
% extracted automatically. Example: xmlFileName = '180501_Hip_Panel3ug_Final-2b.xml';
xmlFileName = ''; % Leave this blank for now

% Output folder for stitched images. Default: current working directory
OutputFolder = ['/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_HiADCase/NoAuBGFFtDenoised_New/no_fftnoise_3X', '/TILING/'];
% 1 for grayscale, 3 for RGB image
img_color = 1;
% Point to folder where the Point folders with TIFs are located. pwd stands for current
% working directory (path at the top of Matlab)
% Example for a different folder: TIFs_PATH = [pwd, '/extracted']
% This relies on your TIFs being inside TIFs_PATH/PointX/TIFs/ folder. To
% change this - modify MibiStitchLoopSupport around line 62
TIFs_PATH = '/Volumes/BryJC_Stanford/Data/Cleaned_Data_Kausi/MedRes_HiADCase/NoAuBGFFtDenoised_New/no_fftnoise_3X'; 
dataSize = 512; % Resolution of one frame

% Stitch start and end points.
startPoint = 1; % Start stitching from this point number
endPoint = 196; % End stitching with this point number. Zero means all points to the end
skipPoints = []; % These points will be skipped during stitching (The stitch should advance, leaving a blank space)

%% Set these if no XML is available
% Set stitching parameters manually if no XML file is available
xNumPoint = 14; % Set this to the number of rows
yNumPoint = 14; % Set this to the number of columns
direction = 1; % Choose starting stitch direction: 0 = left, 1 = right

% new stitch
% X and Y refer to pixel matrix row and column
ydRight = 0; %5 Shift this many pixels when moving right
xdRight = 0; %30 Vertical tilt of the image, shift this many pixels up each time when moving right
ydTop = 0; %15 Shift right by ydTop pixels when moving up one row. Horizontal tilt
xdTop = 0; %-25 Should be negative or 0. Controls vertical coregistration when moving up one row. Positive value would yield blank space between rows

% Gather all channel names
ImageNamesP = dir([TIFs_PATH,'/Point3/TIFs/*tif']);
ImageNamesP = {ImageNamesP.name};
%ImageNamesP = regexp(ImageNamesP, '\w*\.tif');
ImageNames{1, length(ImageNamesP)} = [];
for i = 1:length(ImageNamesP)
    ImageNames{1,i} = ImageNamesP{i};
end
allChannels = ImageNames;

% Overwrite the above to stitch specific channels only
%allChannels = ["dsDNA", "CD56"];

%% Calculating the rest of the offsets and starting the loop
ydRight = dataSize - ydRight; % Adjusting for the relative shift when going right
ydLeft = dataSize - ydRight; % do not change, similar to ydRight
xdLeft = -xdRight; % do not change, similar to xdRight

% Create a weights matrix
weights = zeros(dataSize, dataSize);
for k = 1:floor(dataSize/2)
    weights([k, dataSize - k + 1], k : dataSize - k + 1) = k;
    weights(k : dataSize - k + 1, [k, dataSize - k + 1] ) = k;
end

% Start the loop
for channelNonChar = allChannels
    channel = char(channelNonChar);
    disp(['Making stitched TIF for: ', channel]);
    MibiStitchLoopSupportText_noBorders;
end
disp('Finished stitching all channels');
    


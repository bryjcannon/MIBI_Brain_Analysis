% MIBI smooth stitching script used by the all channel stitch loop
% Author: Dmitry Tebaykin
% Contact: dmitry.tebaykin@stanford.edu
yNumPoint = 17;
if exist(xmlFileName, 'file')
    textXML = fileread(xmlFileName);
    paramNames= {'XAttrib', 'YAttrib'};
    pointsLoc = zeros(0,2);

    for i=1:length(paramNames)
        pattern=[paramNames{i},'="([\+-\w.]+)"\>'];
        [matchExp,tok,ext]= regexp(textXML, pattern, 'match','tokens','tokenExtents');

        for j=1:length(tok)
            pointsLoc(j,i) = str2double(tok{j}{1});
        end
    end

    % Calculate number of rows and cols
    if endPoint == 0
        endPoint = length(pointsLoc);
    end

else
    if endPoint == 0
        endPoint = xNumPoint * yNumPoint;
    end

    if (direction == 0)
        startPosGlobal = ([(xNumPoint - 1/2) * dataSize, (yNumPoint - 1/2) * dataSize]); % Starting point: bottom right
    else
        startPosGlobal = ([(xNumPoint - 1/2) * dataSize, (1/2) * dataSize]); % Starting point: bottom left
    end 
end

% Create list of points for this stitch
pointList = startPoint : endPoint + 1;

allDataStitch = zeros((xNumPoint + 2) * dataSize, (yNumPoint + 2) * dataSize);

% Main stitching loop, this code should not be modified on run-to-run basis
runNumber = 1;
pointsProcessed = 0;
currPoint = 0;

% new stitch
% X and Y refer to pixel matrix row and column
% will make correctin on first run
ydRight = 80; % Shift this many pixels when moving right
xdRight = 50; % Vertical tilt of the image, shift this many pixels up each time when moving right
ydTop = 25; % Shift right by ydTop pixels when moving up one row. Horizontal tilt
xdTop = -85; % Should be negative or 0. Controls vertical coregistration when moving up one row. Positive value would yield blank space between rows

%% Calculating the rest of the offsets and starting the loop
ydRight = dataSize - ydRight; % Adjusting for the relative shift when going right
ydLeft = dataSize - ydRight; % do not change, similar to ydRight
xdLeft = -xdRight; % do not change, similar to xdRight

for i=1:xNumPoint
    if currPoint > 34 %number of points in first run 
            break;
    end
    xloc = xNumPoint - i + 1;
    for j=1:yNumPoint
        yloc = yNumPoint - j + 1;
        currPoint = (i-1) * yNumPoint + j;

        if currPoint > 34 %number of points in first run 
            break;
        end
       
        % Set serpentine direction for even and odd rows
        if (mod(i,2) == 0)
            currentDirection = ~direction;
            if direction == 0
                yloc = j;
            end
        else
            currentDirection = direction; 
            if direction == 1
                yloc = j;
            end
        end
%         if currPoint == 271
%             if (mod(i,2) == 0)
%                 currentDirection = ~direction;
%                 if direction == 1
%                     yloc = j;
%                 end
%             else
%                 currentDirection = direction; 
%                 if direction == 1
%                     yloc = yNumPoint - j + 1;
%                 end
%             end
%         end
        % Get current data frame. First is the run in the bottom.
        currRun = dir([TIFs_PATH{runNumber}, '/Point*']);
        if pointList(currPoint) <= numel(currRun) + pointsProcessed % number of point inthe botom run.
            currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
        else % the consecutive run
            runNumber = runNumber + 1;
            pointsProcessed = pointsProcessed + numel(currRun);
            currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
        end
        currData = currData(3 : dataSize + 2, 3 : dataSize + 2);
        %currData=img(:,:,[1 1 1]);
        
        % Get first position
        if (i == 1) && (j == 1) % first point. No coregistering
            currPos = startPosGlobal;
            allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1) = currData;
            continue;
        end

       % Stitching starts here.
       lastPos = currPos;

       if (currentDirection == 1) && ~(yloc == 1) % registering along right movement
           currPos = ([lastPos(1) + xdRight, lastPos(2) + ydRight]);          
       elseif ((currentDirection == 0) && (yloc == yNumPoint)) || ((currentDirection == 1) && (yloc == 1)) % registration along top movement
           currPos = ([lastPos(1) - dataSize - xdTop, lastPos(2) + ydTop]); % shift coordinates. if more pos moving up or increasing gap between frames. 
       elseif (currentDirection == 0) && ~(yloc == yNumPoint) %registering along left movement
           currPos = ([lastPos(1) + xdLeft, lastPos(2) - ydRight]);
       end

       % Skip the point if needed, leaving blank space in the stitch
       if ismember(currPoint, skipPoints) 
           continue;
       end

       % Add the current point to the stitch matrix, adjusting for existing
       % pixels in the overlap
       prevData = allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1);
       currData = MibiCalcOverlap(prevData, currData, weights);
       allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1) = currData;
   end
end
% values for all point in second run
yNumPoint = 17;
runNumber = runNumber + 1;
currPoint = 0;
ydRight = 80; % Shift this many pixels when moving right. horizontal
xdRight = 50; % Vertical tilt of the image, shift this many pixels up each time when moving right
ydTop = 35; % Shift right by ydTop pixels when moving up one row. Horizontal tilt
xdTop = -110; % Should be negative or 0. Controls vertical coregistration when moving up one row. Positive value would yield blank space between rows
%% Calculating the rest of the offsets and starting the loop
ydRight = dataSize - ydRight; % Adjusting for the relative shift when going right
ydLeft = dataSize - ydRight; % do not change, similar to ydRight
xdLeft = -xdRight; % do not change, similar to xdRight


% Adjust between runs (values for first point in sencond run)
ydTop_run = 210; % Shift right by ydTop pixels when moving up one row. Horizontal tilt
xdTop_run = -135; % Should be negative or 0. Controls vertical coregistration when moving up one row. Positive value would yield blank space between rows

for i=1:xNumPoint
    if currPoint > 51 %number of points in second run 
            break;
    end
    xloc = xNumPoint - i - 1;
    for j=1:yNumPoint
        yloc = yNumPoint - j + 1;
        currPoint = (i-1) * yNumPoint + j;

        if currPoint > 51 %number of points in second run 
            break;
        end
       
        % Set serpentine direction for even and odd rows
        if (mod(i,2) == 0)
            currentDirection = ~direction;
            if direction == 0
                yloc = j;
            end
        else
            currentDirection = direction; 
            if direction == 1
                yloc = j;
            end
        end
%         if currPoint == 271
%             if (mod(i,2) == 0)
%                 currentDirection = ~direction;
%                 if direction == 1
%                     yloc = j;
%                 end
%             else
%                 currentDirection = direction; 
%                 if direction == 1
%                     yloc = yNumPoint - j + 1;
%                 end
%             end
%         end
        % Get current data frame. First is the run in the bottom.
        currRun = dir([TIFs_PATH{runNumber}, '/Point*']);
%         if pointList(currPoint) <= numel(currRun) + pointsProcessed % number of point inthe botom run.
            currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
%         else % the consecutive run
%             runNumber = runNumber + 1;
%             pointsProcessed = pointsProcessed + numel(currRun);
%             currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
%         end
        currData = currData(3 : dataSize + 2, 3 : dataSize + 2);
        %currData=img(:,:,[1 1 1]);
        
%         % Get first position
        if (i == 1) && (j == 1) % first point. No coregistering
            currPos = ([(xloc - 1/2) * dataSize, (1/2) * dataSize]);
            currPos = ([currPos(1) - xdTop_run, currPos(2) + ydTop_run]); % shift coordinates. if more pos moving up or increasing gap between frames. 
            allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1) = currData;
            continue;
        end

       % Stitching starts here.
       lastPos = currPos;

       if (currentDirection == 1) && ~(yloc == 1) % registering along right movement
           currPos = ([lastPos(1) + xdRight, lastPos(2) + ydRight]);          
       elseif ((currentDirection == 0) && (yloc == yNumPoint)) || ((currentDirection == 1) && (yloc == 1)) % registration along top movement
           currPos = ([lastPos(1) - dataSize - xdTop, lastPos(2) + ydTop]); % shift coordinates. if more pos moving up or increasing gap between frames. 
       elseif (currentDirection == 0) && ~(yloc == yNumPoint) %registering along left movement
           currPos = ([lastPos(1) + xdLeft, lastPos(2) - ydRight]);
       end

       % Skip the point if needed, leaving blank space in the stitch
       if ismember(currPoint, skipPoints) 
           continue;
       end

       % Add the current point to the stitch matrix, adjusting for existing
       % pixels in the overlap
       prevData = allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1);
       currData = MibiCalcOverlap(prevData, currData, weights);
       allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1) = currData;
   end
end
% yNumPoint = 2;
% runNumber = runNumber + 1;
% currPoint = 2;
% for i=1:xNumPoint
%     if currPoint < 1
%             break;
%     end
%     xloc = xNumPoint - i - 15;
%     for j=1:yNumPoint
%         yloc = yNumPoint - j + 1;
% %         currPoint = (i-1) * yNumPoint + j;
% 
%         if currPoint < 1
%             break;
%         end
%        
%         % Set serpentine direction for even and odd rows
%         if (mod(i,2) == 0)
%             currentDirection = ~direction;
%             if direction == 0
%                 yloc = j;
%             end
%         else
%             currentDirection = direction; 
%             if direction == 1
%                 yloc = j;
%             end
%         end
% %         if currPoint == 271
% %             if (mod(i,2) == 0)
% %                 currentDirection = ~direction;
% %                 if direction == 1
% %                     yloc = j;
% %                 end
% %             else
% %                 currentDirection = direction; 
% %                 if direction == 1
% %                     yloc = yNumPoint - j + 1;
% %                 end
% %             end
% %         end
%         % Get current data frame. First is the run in the bottom.
%         currRun = dir([TIFs_PATH{runNumber}, '/Point*']);
% %         if pointList(currPoint) <= numel(currRun) + pointsProcessed % number of point inthe botom run.
%             currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
% %         else % the consecutive run
% %             runNumber = runNumber + 1;
% %             pointsProcessed = pointsProcessed + numel(currRun);
% %             currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
% %         end
%         currData = currData(3 : dataSize + 2, 3 : dataSize + 2);
%         %currData=img(:,:,[1 1 1]);
%         
% %         % Get first position
%         if (i == 1) && (j == 1) % first point. No coregistering
%             currPos = ([(xloc + 1/2) * dataSize, (1/2) * dataSize]);
%             allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1) = currData;
%             currPoint = currPoint - 1;
%             continue;
%         end
% 
%        % Stitching starts here.
%        lastPos = currPos;
% 
%        if (currentDirection == 1) && ~(yloc == 1) % registering along right movement
%            currPos = ([lastPos(1) + xdRight, lastPos(2) + ydRight]);          
%        elseif ((currentDirection == 0) && (yloc == yNumPoint)) || ((currentDirection == 1) && (yloc == 1)) % registration along top movement
%            currPos = ([lastPos(1) - dataSize - xdTop, lastPos(2) + ydTop]); % shift coordinates. if more pos moving up or increasing gap between frames. 
%        elseif (currentDirection == 0) && ~(yloc == yNumPoint) %registering along left movement
%            currPos = ([lastPos(1) + xdLeft, lastPos(2) - ydRight]);
%        end
% 
%        % Skip the point if needed, leaving blank space in the stitch
%        if ismember(currPoint, skipPoints) 
%            continue;
%        end
% 
%        % Add the current point to the stitch matrix, adjusting for existing
%        % pixels in the overlap
%        prevData = allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1);
%        currData = MibiCalcOverlap(prevData, currData, weights);
%        allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1) = currData;
%        currPoint = currPoint - 1;
%    end
% end

% plot final image
%figure; 
data=allDataStitch;
% data(data>capImage) = capImage;
% imagesc(data);
% colormap(gray);
% colorbar;
% title(p{1}.massDS.Label{plotChannelInd});T
%set(gca,'xtick',[],'ytick',[]);

if ~exist(OutputFolder,'dir')
    mkdir(OutputFolder);
end

%Make and save images and close figures
imwrite((uint16(data)),[OutputFolder,'/StitchedSeg_',channel]);

%Scaling
capImage=max(max(max(data)));
percentile=capImage*.9;
%data(data>capImage) = capImage;
data = (data/percentile)*4;
%data2(data2(data/capImage)) = capImage;
%imwrite(data,[OutputFolder, '/StitchedScaled_',channel]);


close all;

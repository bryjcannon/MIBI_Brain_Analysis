% MIBI smooth stitching script used by the all channel stitch loop
% Author: Dmitry Tebaykin
% Contact: dmitry.tebaykin@stanford.edu
xNumPoint = 20; 
yNumPoint = 15;
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
        startPosGlobal = ([(xNumPoint - 1/2) * dataSize, (yNumPoint - 5 - 1/2) * dataSize]); % Starting point: bottom right
    else
        startPosGlobal = ([(xNumPoint - 1/2) * dataSize, (5 + 1/2) * dataSize]); % Starting point: bottom left
    end 
end

% Create list of points for this stitch
pointList = startPoint : endPoint + 1;

allDataStitch = zeros((xNumPoint + 2) * dataSize, (yNumPoint + 2) * dataSize, 3);

% Main stitching loop, this code should not be modified on run-to-run basis
runNumber = 1;
pointsProcessed = 0;
currPoint = 1;
% new stitch
% X and Y refer to pixel matrix row and column
% will make correctin on first run
ydRight = 0; % Shift this many pixels when moving right
xdRight = 0; % Vertical tilt of the image, shift this many pixels up each time when moving right
ydTop = 0; % Shift right by ydTop pixels when moving up one row. Horizontal tilt
xdTop = 0; % Should be negative or 0. Controls vertical coregistration when moving up one row. Positive value would yield blank space between rows

%% Calculating the rest of the offsets and starting the loop
ydRight = dataSize - ydRight; % Adjusting for the relative shift when going right
ydLeft = dataSize - ydRight; % do not change, similar to ydRight
xdLeft = -xdRight; % do not change, similar to xdRight
for i=1:xNumPoint
    if currPoint > 100 %number of points in first run
        runNumber = runNumber + 1;
        pointsProcessed = pointsProcessed + numel(currRun);
        break;
    end
    xloc = xNumPoint - i + 1;
    for j=1:yNumPoint
        yloc = yNumPoint - j + 1;
        %currPoint = (i-1) * yNumPoint + j; 
          
        if currPoint == 11 && i==1
            break;
        end
        if currPoint > 100 %number of points in first run 
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
        
        % Get current data frame. First is the run in the bottom.
        currRun = dir([TIFs_PATH{runNumber}, '/Point*']);
        %if pointList(currPoint) <= numel(currRun) + pointsProcessed % number of point inthe botom run.
            currData = (imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
            currPoint = currPoint + 1;
        %else % the consecutive run
            %runNumber = runNumber + 1;
            %pointsProcessed = pointsProcessed + numel(currRun);
            %currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
        %end
        currData = currData(3 : dataSize + 2, 3 : dataSize + 2,:);
        %currData=img(:,:,[1 1 1]);
        
        % Get first position
        if (i == 1) && (j == 1) % first point. No coregistering
            currPos = startPosGlobal;
            %currData = insertText(currData,[20 40],['Set1-Pt',num2str(currPoint-1)],'FontSize',80,'TextColor','red','BoxColor','black','BoxOpacity',0); 
            allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1, :) = currData;
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
       %currData = MibiCalcOverlap(prevData, currData, weights);
       %currData = insertText(currData,[20 40],['Set1-Pt',num2str(currPoint-1)],'FontSize',80,'TextColor','red','BoxColor','black','BoxOpacity',0); 
       allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1, :) = currData;
       
   end
end


%for Bryan's pt numbers
pointList = pointsProcessed + 1:pointsProcessed + 105; 

yNumPoint = 15;
xNumPoint = 13;
currPoint = 0;
pointsProcessed = 0;
direction = 0;
ydRight = 0; % Shift this many pixels when moving right
xdRight = 0; % Vertical tilt of the image, shift this many pixels up each time when moving right
ydTop = 0; % Shift right by ydTop pixels when moving up one row. Horizontal tilt
xdTop = 0; % Should be negative or 0. Controls vertical coregistration when moving up one row. Positive value would yield blank space between rows
%% Calculating the rest of the offsets and starting the loop
ydRight = dataSize - ydRight; % Adjusting for the relative shift when going right
ydLeft = dataSize - ydRight; % do not change, similar to ydRight
xdLeft = -xdRight; % do not change, similar to xdRight

for i=1:xNumPoint
    if currPoint > 105 %number of points in second run 
        runNumber = runNumber + 1;
        pointsProcessed = pointsProcessed + numel(currRun);
        break;
    end
    xloc = xNumPoint - i - 1;
    for j=1:yNumPoint
        yloc = yNumPoint - j + 1;
        currPoint = (i-1) * yNumPoint + j;
       
        if currPoint > 105 %number of points in second run 
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
         %if pointList(currPoint) <= numel(currRun) + pointsProcessed % number of point inthe botom run.
            currData =(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
         %else % the consecutive run
             %runNumber = runNumber + 1;
             %pointsProcessed = pointsProcessed + numel(currRun);
             %currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
         %end
        currData = currData(3 : dataSize + 2, 3 : dataSize + 2,:);
        %currData=img(:,:,[1 1 1]);
        
%         % Get first position
        if (i == 1) && (j == 1) % first point. No coregistering
            currPos = ([(xloc - 1/2 + 2) * dataSize, (1/2 + yNumPoint -1) * dataSize]);
            %currData = insertText(currData,[20 40],['Set2-Pt',num2str(currPoint)],'FontSize',80,'TextColor','red','BoxColor','black','BoxOpacity',0); 
            allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1, :) = currData;
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
       %currData = MibiCalcOverlap(prevData, currData, weights);
       %currData = insertText(currData,[20 40],['Set2-Pt',num2str(currPoint)],'FontSize',80,'TextColor','red','BoxColor','black','BoxOpacity',0); 
       allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1, :) = currData;

   end
end

%for Bryan's pt numbers
pointList = 206:275;

yNumPoint = 15;
xNumPoint = 6;
currPoint = 1;
pointsProcessed = 0;
direction = 1;
for i=1:xNumPoint
    if currPoint > 70
            break;
    end
    xloc = xNumPoint - i + 1;
    for j=1:yNumPoint
        yloc = yNumPoint - j + 1;
        %currPoint = (i-1) * yNumPoint + j;
       
        if (currPoint == 15 && i==1) || (currPoint == 28 && i==2) || (currPoint == 40 && i==3)|| (currPoint == 51 && i==4) || (currPoint == 61 && i==5)
            break;
        end
        
        if currPoint > 70
            break;
        end
       
        % Set serpentine direction for even and odd rows
        if (mod(i,2) == 0)
            currentDirection = ~direction;
            if direction == 0
                yloc = j + i;
            end
        else
            currentDirection = direction; 
            if direction == 1
                yloc = j + i;
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
        %if pointList(currPoint) <= numel(currRun) + pointsProcessed % number of point inthe botom run.
            currData = (imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
            currPoint = currPoint + 1;
            %else % the consecutive run
%              runNumber = runNumber + 1;
%              pointsProcessed = pointsProcessed + numel(currRun);
%              currData = double(imread([TIFs_PATH{runNumber}, '/Point', num2str(pointList(currPoint) - pointsProcessed), '/TIFs/', channel]));
%          end
        currData = currData(3 : dataSize + 2, 3 : dataSize + 2,:);
        %currData=img(:,:,[1 1 1]);
        
%         % Get first position
        if (i == 1) && (j == 1) % first point. No coregistering
            currPos = ([(xloc - 1/2) * dataSize, (1/2 + 1) * dataSize]);
            %currData = insertText(currData,[20 40],['Set3-Pt',num2str(currPoint-1)],'FontSize',80,'TextColor','red','BoxColor','black','BoxOpacity',0); 
            allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1, :) = currData;
            continue;
        end

       % Stitching starts here.
       lastPos = currPos;

       if (currentDirection == 1) && ~(yloc == i + 1) % registering along right movement
           currPos = ([lastPos(1) + xdRight, lastPos(2) + ydRight]);          
       elseif ((currentDirection == 0) && (yloc == yNumPoint)) 
           currPos = ([lastPos(1) - dataSize - xdTop, lastPos(2) + ydTop]); % shift coordinates. if more pos moving up or increasing gap between frames. 
       elseif ((currentDirection == 1) && (yloc == i + 1)) % registration along top movement
           currPos = ([lastPos(1) - dataSize - xdTop, lastPos(2)+ dataSize + ydTop]); % shift coordinates. if more pos moving up or increasing gap between frames.
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
       %currData = MibiCalcOverlap(prevData, currData, weights);
       %currData = insertText(currData,[20 40],['Set3-Pt',num2str(currPoint-1)],'FontSize',80,'TextColor','red','BoxColor','black','BoxOpacity',0); 
       allDataStitch(currPos(1) : currPos(1) + dataSize - 1, currPos(2) : currPos(2) + dataSize - 1, :) = currData;
       
    end
end

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
imwrite((uint8(data)),[OutputFolder,'/',channel]);

%Scaling
capImage=max(max(max(data)));
percentile=capImage*.9;
%data(data>capImage) = capImage;
data = (data/percentile)*4;
%data2(data2(data/capImage)) = capImage;
%imwrite(data,[OutputFolder, '/StitchedScaled_',channel]);


close all;

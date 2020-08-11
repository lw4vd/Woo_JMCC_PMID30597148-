function [nucData,cellData] = analyzeWell(wellNum,channels,segThresholdLevel,minCellArea,maxCellArea)
% analyzeWell.m
% Called by 'Master.m' to analyze images in a particular well.
% History:
% Dec 2011, Jason Yang: original code
% Jan 2012, Jeff Saucerman: restructuring
% Sep 2013, Philip Tan: added error messages and ability to work with TIFFs
% 2017, Laura Woo: adapted code for binucleation and ploidy analysis

%% Load images and filter images for each channel
for chNum = 1:numel(channels)
    filename = sprintf('Well%0.2d_ch%0.2d.tif',wellNum,chNum);
    disp([filename]);
    
    images{chNum} = double(imread([filename])); % load tiff files
    imageInfo = imfinfo(filename); 
    
    if isempty(images{chNum})
        disp([' image ' [filename] ' not successfully loaded.']);
    end
    
    images{chNum} = wiener2(images{chNum}); % filter image
end

%% Segment objects in the first channel
disp('segmenting...'); tic
[segmentedImageNuc,segmentedImageCell] = segmentObjects(images{1},segThresholdLevel,minCellArea,maxCellArea); toc

% TEMP
assignin('base','segmentedImagesNuc',segmentedImageNuc);
assignin('base','segmentedImagesCell',segmentedImageCell);
assignin('base','images',images);

%% Measure objects 
disp('measuring...'); tic
measurementsNuc = measureObjects(images,segmentedImageNuc,imageInfo); 
measurementsCell = measureObjects(images,segmentedImageCell,imageInfo);toc

%% Arrange data
disp('organizing data...'); tic
assignin('base','segmentedImagesNuc',segmentedImageNuc);
assignin('base','segmentedImagesCell',segmentedImageCell);

numNuc = max(segmentedImageNuc(:));
nucLabels = {[1:numNuc]'};
numCells = max(segmentedImageCell(:));
cellLabels = {[1:numCells]'}; toc

if numNuc > 0
nucData = arrangeData(measurementsNuc,nucLabels);
cellData = arrangeData(measurementsCell,cellLabels);
else
    nucData = 0; 
    cellData = 0; 
end

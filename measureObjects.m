function [measurements] = measureObjects(images,segmentedImage,imageInfo)
% measureCells.m
% images: a cell array containing an image for each channel
% segmentedImage: output image from either bwlabel or bwconncomp
% Note: shape measurements are returned in units of pixels
%
% History:
% Dec 2011, Jason Yang: original code
% Jan 2012, Jeff Saucerman: restructuring
% 5/30/2012, JS: added background subtraction using image mode
% 2016, Laura Woo: modified to prevent background = saturated pixel intensity
% 2017, Laura Woo: adapted code for binucleation and ploidy analysis

%% Shape Measurements
stats = regionprops(segmentedImage,'Area','Perimeter');
nObjects = numel(stats);

area = cat(1,stats.Area);
perimeter = cat(1,stats.Perimeter);
formFactor = 4*pi.*area./perimeter./perimeter;

%% Pixel Intensity Measurements for each channel
numCh = numel(images);
meanIntensity = zeros(nObjects,numCh);
integratedIntensity = zeros(nObjects,numCh);

for chNum=1:numCh
    background = mode(reshape(images{chNum},[],1));% Image background calculated as mode of that channel
    if background == imageInfo.MaxSampleValue
        imageBS = images{chNum};
        [Hist,Bin] = histc(imageBS(:),0:10:(imageInfo.MaxSampleValue+5));
        background = imageInfo.MaxSampleValue/numel(Hist)*mode(Bin);
    end

    stats = regionprops(segmentedImage,images{chNum}-background,'PixelValues');
    X = struct2cell(stats);
    
    for i = 1:nObjects
        meanIntensity(i,chNum) = mean(X{i});
        integratedIntensity(i,chNum) = sum(X{i});
    end
end

measurements = struct('area',area,'perimeter',perimeter,'formFactor',formFactor,...
    'meanIntensity',meanIntensity,'integratedIntensity',integratedIntensity);
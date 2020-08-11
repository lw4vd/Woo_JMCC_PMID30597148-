function [segmentedImageNuc,segmentedImageCell] = segmentObjects(X,segThresholdLevel,minCellArea,maxCellArea)
% segmentObjects.m
% History:
% Dec 2011, Jason Yang: original code
% Jan 2012, Jeff Saucerman: restructuring
% Sep 2013, Philip Tan: prevent mode from being 1 in saturated images
% 2017, Laura Woo: adapted code for binucleation and ploidy analysis

dilateSize = 3; %pixels
splitThresh = 1; %pixels
BNOutlineDilateSize = 1; %normally 2 to visualize difference

% Rescale image from 0 to 1
Xoriginal = X;
X = X-min(X(:));
X = X./max(X(:)); 

modeX = mode(X(:));
if modeX == 1
    modeX = mode(X(X~=mode(X(:))));
end

threshold = modeX+segThresholdLevel*std(X(:));
maskBW = im2bw(X,threshold);
maskBW = imfill(maskBW,'holes');
maskBW = bwareaopen(maskBW, minCellArea); 

D = -bwdist(~maskBW);
minima = imextendedmin(D, splitThresh);  
D = imimposemin(D, minima); 
WS = watershed(D)>0;
maskBW2 = maskBW;
maskBW2(WS == 0) = 0; 
objects = maskBW2.*WS;
SE = strel('disk',dilateSize);
objects = imopen(objects,SE);
segmentedImageNuc = bwlabel(objects);

% Filter nuclei by size
stats = regionprops(segmentedImageNuc,'Area');
area = cat(1,stats.Area);
indBigCells = find(area>=minCellArea & area<=maxCellArea);

% Filter nuclei on border:
mask = imclearborder(objects);
mask(~ismember(segmentedImageNuc,indBigCells)) = 0;
segmentedImageNuc = bwlabel(mask);

%segment binucleated cells
SE = strel('disk',dilateSize);  
objects = imclose(segmentedImageNuc,SE);
objects = imfill(objects,'holes');
segmentedImageCell = bwlabel(objects);

% Save BN outline
BWOutlineMN = bwperim(segmentedImageNuc);
BWOutlineBN = bwperim(segmentedImageCell); 
SE = strel('disk',BNOutlineDilateSize);
BWOutlineBN = imdilate(BWOutlineBN,SE);

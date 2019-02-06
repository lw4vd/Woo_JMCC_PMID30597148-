function [segmentedImageNuc,segmentedImageCell] = segmentObjects(X,segThresholdLevel,localMaximaThreshold,minCellArea,maxCellArea)
% segmentObjects.m
% History:
% Dec 2011, Jason Yang: original code
% Jan 2012, Jeff Saucerman: restructuring
% Sep 2013, Philip Tan: prevent mode from being 1 in saturated images
% 2017, Laura Woo: adapted code for binucleation and ploidy analysis

% Rescale image from 0 to 1
Xoriginal = X;
X = X-min(X(:));
X = X./max(X(:));

SE = strel('disk',20);
XTH = imtophat(X,SE);

modeX = mode(XTH(:));
if modeX == 1
    modeX = mode(XTH(XTH~=mode(XTH(:))));
end
threshold = modeX+segThresholdLevel*std(XTH(:));
mask = im2bw(XTH,threshold);
mask = imfill(mask,'holes');

X1 = XTH;
maxima = imextendedmax(X1,localMaximaThreshold);
maxima = imfill(maxima,'holes');
overlay = imimposemin(1-X1,maxima);
WS = watershed(overlay)>0;

objects = mask.*WS;
SE = strel('disk',5);
objects = imopen(objects,SE);

segmentedImageNuc = bwlabel(objects);
BWOutline = bwperim(segmentedImageNuc);
OutlineOverlay = 3.*X + BWOutline;

FFList = regionprops(segmentedImageNuc, 'Area','Perimeter');
AR = cat(1,FFList.Area);
PR = cat(1,FFList.Perimeter);
FF = 4*pi.*AR./PR./PR;

multiNuc = find(FF<0.75 & AR< 2500);
mask2 = ismember(segmentedImageNuc,multiNuc);
maskSingleNuc = segmentedImageNuc.*~mask2;
X3 = XTH.*mask2;

bw = mask2.*X3;
D = bwdist(~bw);
maxima = imextendedmax(D,1);
maxima = imfill(maxima,'holes');

overlay = imimposemin(1-X3,maxima);
WS = watershed(overlay)>0;
mask2 = mask2.*WS;
SE = strel('disk',5);  
objects = imopen(mask2,SE);

SingleNuc = maskSingleNuc + objects;
segmentedImageNuc = bwlabel(SingleNuc);
BWOutline = bwperim(segmentedImageNuc);

%segment binucleated cells
SE = strel('disk',5);  
objects = imclose(SingleNuc,SE);
objects = imfill(objects,'holes');
segmentedImageCell = bwlabel(objects);

% Filter nuclei by size
stats = regionprops(segmentedImageNuc,'Area');
area = cat(1,stats.Area);
indBigCells = find(area>=minCellArea & area<=maxCellArea);

% Filter nuclei on border:
mask = imclearborder(SingleNuc);
mask(~ismember(segmentedImageNuc,indBigCells)) = 0;
segmentedImageNuc = bwlabel(mask);

% Filter cells by size
statsBN = regionprops(segmentedImageCell,'Area');
areaBN = cat(1,statsBN.Area);
indBigCellsBN = find(areaBN>=minCellArea & areaBN<=(maxCellArea+10000));

% Filter cells on border
maskBN = imclearborder(objects);
maskBN(~ismember(segmentedImageCell,indBigCellsBN)) = 0;
segmentedImageCell = bwlabel(maskBN);
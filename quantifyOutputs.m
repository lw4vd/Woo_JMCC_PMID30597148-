% quantify # of nuclei, # of cells, % binucleated cells, % 2N nuclei, % 4N nuclei, % >4N nuclei

clear all; clc; close all;

%inputs: 
load('nucOutput_well01.mat');
load('cellOutput_well01.mat');

% quantify # of nuclei:
meanInt = zeros(numel(nucData),1);
intInt = meanInt;
areaInt = meanInt;

for a = 1:numel(nucData)
    meanInt(a) = nucData(a).meanIntensity;
    intInt(a) = nucData(a).integratedIntensity;
    areaInt(a) = nucData(a).area;
end

pdArea = makedist('Normal',mean(areaInt(:)),std(areaInt(:))); %probability distribution
ThreshArea = icdf(pdArea,0.75)+3*(icdf(pdArea,0.75)-icdf(pdArea,0.25));
pdInt = makedist('Normal',mean(intInt(:)),std(intInt(:)));
ThreshInt = icdf(pdInt,0.75)+3*(icdf(pdInt,0.75)-icdf(pdInt,0.25));

% nuclear ploidy analysis:
intIntList = intInt(find(intInt(:) < ThreshInt & areaInt < ThreshArea));
areaList = areaInt(find(intInt(:) < ThreshInt & areaInt < ThreshArea));
A(:,1) = intIntList;
A(:,2) = areaList;

minc = min(A(:,1));
midc = abs(max(A(:,1))-abs(min(A(:,1))))/2;
maxc = max(A(:,1));
initialCent = [minc; midc; maxc];
[idx,Cent] = kmeans(A(:,1),3,'Start',initialCent);
G1 = intIntList(find(idx == 1));
G2 = intIntList(find(idx == 2));
G3 = intIntList(find(idx == 3));
[SortG,SortGind] = sort(Cent);

G1Thresh = abs(SortG(1)-SortG(2))/2+SortG(1);
G2Thresh = abs(SortG(1)-SortG(2))/2+SortG(2);

nNuclei = numel(find(areaInt(:) < ThreshArea & intInt(:) < ThreshInt)); % # of nuclei
p2N = numel(find(areaInt(:) < ThreshArea & intInt(:) < ThreshInt & intInt(:) < G1Thresh))/nNuclei*100; % perecent of nuclei with 2N ploidy
p4N = numel(find(areaInt(:) < ThreshArea & intInt(:) < ThreshInt & intInt(:) >= G1Thresh & intInt(:) < G2Thresh))/nNuclei*100; % perecent of nuclei with 4N ploidy
pMultiN = numel(find(areaInt(:) < ThreshArea & intInt(:) < ThreshInt & intInt(:) >= G2Thresh))/nNuclei*100; % perecent of nuclei with >4N ploidy

% quantify # of cells :
meanInt = zeros(numel(cellData),1);
intInt = meanInt;
areaInt = meanInt;

for a = 1:numel(cellData)
    meanInt(a) = cellData(a).meanIntensity;
    intInt(a) = cellData(a).integratedIntensity;
    areaInt(a) = cellData(a).area;
end

pdArea = makedist('Normal',mean(areaInt(:)),std(areaInt(:))); %probability distribution
ThreshArea = icdf(pdArea,0.75)+3*(icdf(pdArea,0.75)-icdf(pdArea,0.25));
pdInt = makedist('Normal',mean(intInt(:)),std(intInt(:)));
ThreshInt = icdf(pdInt,0.75)+3*(icdf(pdInt,0.75)-icdf(pdInt,0.25));

nCells = numel(find(areaInt(:) < ThreshArea & intInt(:) < ThreshInt)); % # of cells
pBinucleated = (nNuclei-nCells)/nCells*100; % percent binucleated cells

quantOutputs = struct('nNuclei',nNuclei,'nCells',nCells,'pBinucleated',pBinucleated,'p2N',p2N,'p4N',p4N,'pMultiN',pMultiN)

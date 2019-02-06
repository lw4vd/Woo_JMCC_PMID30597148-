% Master.m
% Top level function used to run batch analysis of cells from multiple wells and channels.
% History:
% Dec 2011, Jason Yang: original code
% Jan 2012, Jeff Saucerman: restructuring
% 5/30/2012, JS: updated segmentCells to use imextendedmax to provide seeds
% for watershed segmentation. 
% Sep 2013, Philip Tan: added error messages and ability to work with TIFFs
% 2017, Laura Woo: adapted code for binucleation and ploidy analysis

clear all; clc; close all;

%% Specify parameters for analysis

wells = [1:1]; % range of wells to analyze
channels = [1:1]; % range of channels to analyze

% segmentation parameters
segThresholdLevel = 1.5; % # of std deviations for object/background thresholding (1.5 is 93.32% of noise)
localMaximaThreshold = 0.1; % size of local maxima for dividing clumped objects (range 0-1, 0.2 is typical, 0.1 more aggressive)
minCellArea = 100;
maxCellArea = 2500;


%% loop over each well
for i=1:numel(wells) 
    wellNum = wells(i);
    disp(['processing well ' num2str(wells(i))]);
    
    [nucData,cellData] = analyzeWell(wellNum,channels,segThresholdLevel,localMaximaThreshold,minCellArea,maxCellArea);

    savefile = sprintf('nucOutput_well%0.2d',wellNum); % save data for nuclei analysis
    save([savefile],'nucData');
    
    savefile = sprintf('cellOutput_well%0.2d',wellNum); % save data for cell analysis
    save([savefile],'cellData');
end
function objectData = arrangeData(measurements,objectLabels)
% arrangeCellData.m
% Creates a structure containing cell data of the format:
%   objectData(objectNum).area or
%   objectData(objectNum).meanIntensity(channel)
% History:
% Jan 2012, Jeff Saucerman: original code
% 2017, Laura Woo: adapted code for binucleation and ploidy analysis

numObjects = max(objectLabels{end}); % total number of objects
numCh = size(measurements(1).meanIntensity,2); % number of channels for intensity measurements

% loop over objects

for objectNum=1:numObjects
    
    objectPos = find(objectLabels{:}==objectNum); % find the position of the object in objectLabels
    if ~isempty(objectPos)
        
        % assign cell shape features
        objectData(objectNum).area = measurements.area(objectPos);
        objectData(objectNum).perimeter = measurements.perimeter(objectPos);
        objectData(objectNum).formFactor = measurements.formFactor(objectPos);
        
        % assign cell intensity features
        % note: "measurements" stores intensity channels like meanIntensity(cellPos,channel),
        % so we need to collect from all channels.
        objectData(objectNum).meanIntensity(1:numCh) = measurements.meanIntensity(objectPos,1:numCh);
        objectData(objectNum).integratedIntensity(1:numCh) = measurements.integratedIntensity(objectPos,1:numCh);
    else
        % assign cell shape features
        objectData(objectNum).area = NaN;
        objectData(objectNum).perimeter = NaN;
        objectData(objectNum).formFactor = NaN;
        
        % assign cell intensity features
        objectData(objectNum).meanIntensity(1:numCh) = NaN;
        objectData(objectNum).integratedIntensity(1:numCh) = NaN;
    end
end % loop over objects
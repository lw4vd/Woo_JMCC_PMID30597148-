Instructions: 
- modify parameters in masterRun.m and run masterRun.m followed by quantifyOutput.m
- example image "Well01_ch01.tif" 
- segmentation results of example image "Well01_ch01_segmented.tif" (white outline = nuclei segmentation results; red outline = binculeation cell segmentation results)  

1) masterRun.m
define paramters: 
- wells = range of wells to analyze
- channels = range of channels to analyze
- segThresholdLevel = # of std for obj/background thresholding
- localMaximaThreshold = size of local maxima for dividing clumped objects
- minCellArea = minimum area of object to segment
- maxCellArea = maximum area of object to segment

*image file name format = "Well##_ch##.tif" (example image included in the zip folder: Well01_ch01.tif); Hoechst/DAPI image = ch01 

outputs: .mat files containing measurements for all segmented nuclei ("nucOutput_well##.mat") and cells ("cellOutput_well##.mat")

2) analyzeWell.m (function called by masterRun.m)
inputs: all inputs defined in masterRun
- wellNum = well # to analyze 
- channels = range of channels to analyze 
- segThresholdLevel = # of std for obj/background thresholding
- localMaximaThreshold = size of local maxima for dividing clumped objects
- minCellArea = minimum area of object to segment
- maxCellArea = maximum area of object to segment

outputs: measurements for all segmented nuclei (nucData) and cells (cellData) 

3) segmentObjects.m (function called by analyzeWell.m to segment nuclei or cells)
inputs: all inputs defined in analyzeWell.m or masterRun.m
- X = DAPI/Hoechst (channel 1) image loaded into MATLAB (from analyzeWell.m)
- segThresholdLevel = # of std for obj/background thresholding (from masterRun.m)
- localMaximaThreshold = size of local maxima for dividing clumped objects (from masterRun.m)
- minCellArea = minimum area of object to segment (from masterRun.m)
- maxCellArea = maximum area of object to segment (from masterRun.m) 

outputs: segmented nuclei and cell images stored in "segmentedImageNuc" and "segmentedImageCell" variables

4) measureObjects.m (function called by analyzeWell.m that makes shape/intensity measurements of segmented object) 
inputs: all inputs defined in analyzeWell.m or segmentObjects.m
- images = all images in well loaded into MATLAB (from analyzeWell.m)
- segmentedImage = segmented image (output from segmentObjects.m in analyzeWell.m)
- imageInfo = information about image file (from analyzeWell.m)

outputs: structure array ("measurements) containing various shape/intensity measurements of each segmented object 
- shape measurements = area, perimeter, formFactor
- pixel intensity measurements = meanIntensity, integratedIntensity

5) arrangeData.m (function called by analyzeWell.m to arrange measurements in a specific format) 
inputs: all inputs defined in analyzeWell.m or measureObjects.m
- measurements = structure array of measurements (output from measureObjects.m in analyzeWell.m)
- objectLabels = object labels (from analyzeWell.m)

outputs: structure containing nuclei/cell data in a specific format
	objectData(objectNum).area(channel)
	objectData(objectNum).meanIntensity(channel)
	objectData(objectNum).integratedIntensity(channel)

6) quantifyOutputs.m (quantifies nuclei/cells from measurements output in masterRun.m; including ploidy analysis)
inputs: 
- .mat files containing measurements for all segmented nuclei ("nucOutput_well##.mat") created in masterRun.m
- .mat file containing measurements for all segmented cells ("cellOutput_well##.mat") createdin masterRun.m

outputs: structure array (quantOutputs) containing multiple quantified outputs: 
# of nuclei (nNuclei), # of cells (nCells), % of binucleated cells (pBinucleated), 
% of 2N nuclei (p2N), % of 4N nuclei (p4N), % of >4N nuclei (pMultiN)

*ploidy analysis performs k-means clustering with 3 clusters (assumes population contains 2N, 4N, >4N nuclei)
*binucleated cells quantification assumes all cells are binucleated (although some may be multinucleated)
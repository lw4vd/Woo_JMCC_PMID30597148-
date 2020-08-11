function [ploidyThreshAll] = findPloidyThresholds(dataDirectory,imageFileNames,ch,sFit,kernelWidthBins)


%find DNA content thresholds by fitting 2 gaussain distributions
%DISCLAIMER: pipeline is robust to clean images of hiPSC-CMs treated with
%inhibitors in AZ omics paper for 6 days and stained with DAPI (on Operetta) and Hoechst
%on Olympus microscope

%INPUT VARS _______
%dataDirectory = directory with processed .mat files (containing nucData from masterRun.m)
%imageFileNames = file names for processed .mat files
%ch =  DAPI channel position number in nucData
%sFit = 0 or 1; if using hoechst data set sFit to 0 (don't fit "S-phase" uniform
%distribution) 
%kernelWidthBins = 1 (DAPI); 2 (Hoechst)

%OUTPUT TABLE ______
%output = table of ploidy thresholds for each image
%threshLG1 = left threshold for 2c population (G1) 
%threshRG1 = right threshold for 2c population (G1) & left threshold for
%S-phase (if sFit = 1)
%threshLG2 = left threshold for 4c population (G2/M) & right threhold for
%S-phase (if sFit = 1) 
%threshRG2 = right trheshold for 4c population (G2/M) & left threshold for
%>4c population
%muG1 = mean of 2c gaussian dist
%muG2 = mean of 4c gaussian dist
%muRatio = ratio of muG2/muG1 (should be ~2) 
%cvG1 = coefficient of variation for 2c dist 
%cvG2 = coefficient of variation for 4c dist

clear ThreshList threshLG1 threshRG1 threshLG2 threshRG2 muG1 muG2 muRatio cvG1 cvG2

stdevG1 = 2; %# of stdev right of G1 peak to measure uniform distribution S amplitude
stdevG2 = 2; %# of stdev left of G2 peak to measure uniform distribution S amplitude
percentSigma = 0.5; %percent of peak amplitude to measure sigma estimate (in fraction not %)
maxY = 50;

for i = 1:numel(imageFileNames)
    clear nucData meanInt intInt areaObj centroidObj 
    
    %load file data
    load([dataDirectory '\' imageFileNames{i}]);
    
    meanInt = zeros(numel(nucData),1);
    intInt = meanInt;
    areaObj = meanInt;
    centroidObj = zeros(numel(nucData),2);
    
    for j = 1:numel(nucData)
        meanInt(j,:) = nucData(j).meanIntensity(ch);
        intInt(j,:) = nucData(j).integratedIntensity(ch);
        areaObj(j,:) = nucData(j).area;
        centroidObj(j,:) = nucData(j).centroid;
    end
    
    nbins = 100;
    
    %% DNA content analysis:
    
    [nCounts,binEdges] = histcounts(intInt, nbins);
    binWidth = binEdges(2)-binEdges(1);
    binCenters = binEdges(1:end-1)+binWidth/2;
    kernelWidth = binWidth*kernelWidthBins;
    kernelFit = fitdist(intInt, 'Kernel', 'Kernel', 'normal', 'Bandwidth', kernelWidth);
    
    xKernel = binCenters;
    yKernel = pdf(kernelFit, xKernel)*kernelWidth*length(intInt)/kernelWidthBins;

    peakAll = sort(findpeaks(yKernel), 'descend');
    
    if ~isempty(peakAll) && sum(nCounts) > 50 %if no peaks then not worth doing ploidy analysis
        if length(peakAll) >= 2
            peak1 = peakAll(1);
            peak2 = peakAll(2);
        else
            peak1 = peakAll(1);
            peak2 = peak1;
        end
        
        peak1_x = xKernel(find(yKernel == peak1));
        peak1_x = peak1_x(1);
        peak2_x = xKernel(find(yKernel == peak2));
        peak2_x = peak2_x(1);
        
        if peak1_x > peak2_x && abs(log2(peak1/peak2)) < 0.5 
            temp = peak2_x;
            peak2_x = peak1_x;
            peak1_x = temp;
            
            peak1 = peakAll(2);
            peak2 = peakAll(1);
        end
        
        ratioList = peak2_x/peak1_x;
        
        %fit G1
        peak1_ind = find(binCenters == peak1_x);
        leftPeak1_x = binCenters(1:peak1_ind(1));
        leftPeak1_y = yKernel(1:peak1_ind(1));
        sigma1_est_ind = find(leftPeak1_y <= peak1*percentSigma); %find ind of 25% max
        if isempty(sigma1_est_ind)
            rightPeak1_x = binCenters(peak1_ind(1):end);
            rightPeak1_y = yKernel(peak1_ind(1):end);
            sigma1_est_ind = find(rightPeak1_y <= peak1*percentSigma); %find ind of 25% max
            sigma1_est = abs(peak1_x-rightPeak1_x(sigma1_est_ind(1)));
        else
            sigma1_est = peak1_x-leftPeak1_x(sigma1_est_ind(end));
        end
        
        peak1_xRange = binCenters(find(binCenters >= peak1_x-3*sigma1_est & binCenters <= peak1_x + 1*sigma1_est));
        peak1_yRange = nCounts(find(binCenters >= peak1_x-3*sigma1_est & binCenters <= peak1_x + 1*sigma1_est));
        
        f1= @(param, xdata)param(3).*exp((-1/2)*((xdata-param(1))./param(2)).^2);
        param0 = [peak1_x, sigma1_est, peak1];
        lb = [peak1_xRange(1), sigma1_est*0.5, peak1*0.5];
        ub = [peak1_xRange(end), sigma1_est*1.1, peak1*1.1];
        [paramFit, resnorm, ~, exitflag, output] = lsqcurvefit(f1, param0, peak1_xRange, peak1_yRange, lb, ub);
        param1 = [paramFit(1), paramFit(2), paramFit(3)];
        
        fitG1_mu = paramFit(1);
        fitG1_sigma = paramFit(2);
        
        %fit second peak
        
        if ratioList >= 1.8 && ratioList <= 2.2
            
            peak2_ind = find(binCenters == peak2_x);
            rightPeak2_x = binCenters(peak2_ind(1):end);
            rightPeak2_y = yKernel(peak2_ind(1):end);
            sigma2_est_ind = find(rightPeak2_y <= peak2*percentSigma); %find ind of 25% max
            sigma2_est = rightPeak2_x(sigma2_est_ind(1)) - peak2_x;
            
            peak2_xRange = binCenters(find(binCenters <= peak2_x+3*sigma2_est & binCenters >= peak2_x - 1*sigma2_est));
            peak2_yRange = nCounts(find(binCenters <= peak2_x+3*sigma2_est & binCenters >= peak2_x - 1*sigma2_est));
            
            param0 = [peak2_x, sigma2_est, peak2];
            lb = [fitG1_mu*1.9, sigma2_est*0.5, peak2*0.5];
            ub = [fitG1_mu*2.1, sigma2_est*1.1, peak2*1.1];
            
        else
            
            peak2_x_est = fitG1_mu*2;
            sigma2_est = fitG1_sigma;
            peak2_est = nCounts(find(binCenters >= peak2_x_est - sigma2_est & binCenters <= peak2_x_est + sigma2_est));
            peak2_est = max(peak2_est);
            
            peak2_xRange = binCenters(find(binCenters <= peak2_x_est+3*sigma2_est & binCenters >= peak2_x_est - 1*sigma2_est));
            peak2_yRange = nCounts(find(binCenters <= peak2_x_est+3*sigma2_est & binCenters >= peak2_x_est - 1*sigma2_est));
            
            param0 = [peak2_x_est, sigma2_est, peak2_est];
            lb = [fitG1_mu*1.9, sigma2_est*0.5, peak2_est*0.5];
            ub = [fitG1_mu*2.1, sigma2_est*1.1, peak2_est*1.1];
            
        end
        
        [paramFit, resnorm, ~, exitflag, output] = lsqcurvefit(f1, param0, peak2_xRange, peak2_yRange, lb, ub);
        param2 = [paramFit(1), paramFit(2), paramFit(3)];
        
        fitG2_mu = paramFit(1);
        fitG2_sigma = paramFit(2);
        
        %% find thresholds for G1 & G2
        if sFit == 1
            % fit S uniform distribution:
            Sparam_est = mean(nCounts(find(binCenters > fitG1_mu + stdevG1*fitG1_sigma & binCenters < fitG2_mu - stdevG2*fitG2_sigma)));
            
            intersectRange = fitG1_mu:1:fitG2_mu;
            intersectRange_y = nCounts(find(binCenters > fitG1_mu & binCenters < fitG2_mu));
            
            intersectInd = find(f1(param1, intersectRange)-Sparam_est <= 0); %first instance it's negative
            if ~isempty(intersectInd)
                threshRightG1 = intersectRange(intersectInd(1));
            else
                threshRightG1 = fitG1_mu + 1*fitG1_sigma;
            end
            
            
            intersectInd = find(f1(param2,intersectRange)-Sparam_est >= 0); %first instance it's positive
            if ~isempty(intersectInd)
                threshLeftG2 = intersectRange(intersectInd(1));
            else
                threshLeftG2 = fitG2_mu - 1*fitG2_sigma;
            end
            
            threshLeftG1 = fitG1_mu-(threshRightG1-fitG1_mu);
            threshRightG2 = fitG2_mu+(fitG2_mu-threshLeftG2);
        else
            
            threshRightG1 = fitG1_mu + (fitG2_mu-fitG1_mu)/2;
            threshLeftG2 = threshRightG1;
            
            threshLeftG1 = fitG1_mu-(threshRightG1-fitG1_mu);
            threshRightG2 = fitG2_mu+(fitG2_mu-threshLeftG2);
            
        end
        
       
    else
        threshLeftG1 = binCenters(find(nCounts == max(nCounts)));
        threshLeftG1 = threshLeftG1(1);
        threshRightG1 = threshLeftG1;
        threshLeftG2 = threshLeftG1;
        threshRightG2 = threshLeftG1;
    end
    
    if ~isempty(peakAll) && sum(nCounts) > 50
        muG1(i) = fitG1_mu;
        muG2(i) = fitG2_mu;
        muRatio(i) = fitG2_mu/fitG1_mu;
        cvG1(i) = fitG1_sigma/fitG1_mu*100; 
        cvG2(i) = fitG2_sigma/fitG2_mu*100; 
    else
        muG1(i) = 0;
        muG2(i) = 0;
        muRatio(i) = 0;
        cvG1(i) = 0; 
        cvG2(i) = 0;
    end
    
    threshLG1(i) = threshLeftG1;
    threshRG1(i) = threshRightG1;
    threshLG2(i) = threshLeftG2;
    threshRG2(i) = threshRightG2;
    
end %by plate

ploidyThreshAll = array2table([threshLG1', threshRG1', threshLG2', threshRG2', muG1', muG2', muRatio', cvG1', cvG2'], ...
    'VariableNames', {'threshLG1', 'threshRG1', 'threshLG2', 'threshRG2', 'muG1', 'muG2', 'muRatio','cvG1', 'cvG2'});
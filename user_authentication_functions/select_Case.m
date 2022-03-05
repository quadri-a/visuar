function [ peaks, numChirps, GRTruth, sigEnergy, frPeaks, unwrap_phase_mm, fileName, framePeriod, timeunit, unwrp_phase_max_peaks, ph_diff, ranges] = select_Case(caseNum, viewRanges, TargetMinIdx, TargetMaxIdx)%, sampMin, sampMax)

dataType = caseNum;

switch dataType
% 2 Loops 50 ms Dataset =++++++++++++++++===============================

    case 1 % ADC 100 Best Range - Best Range -  HR -Mean/MEdian GR Truth 90
        fileName ='cases/user_1.bin'; fs = 2000e3; S = 70e6; adcsamples = 100; GRTruth=90; framePeriod = 50e-3;numChirps =2; TargetMinIdx= TargetMinIdx; TargetMaxIdx = TargetMaxIdx; 
        [ adcData] = readDCA1000_BR_holdBreathing(fileName, adcsamples); plotEnable = viewRanges;
%         adcData = adcData(1:12000,:);
         adcData = adcData(12000:24000,:);
        [ sigEnergy, ranges,timeunit,peakVals, peaks, frPeaks] = rangeBinSelection(adcData, adcsamples, fs, S, framePeriod,numChirps, TargetMinIdx, TargetMaxIdx,  fileName, plotEnable);

    case 2 % ADC 100 Best Range - Best Range -  HR -Mean/MEdian GR Truth 70
        fileName ='cases/user_2.bin'; fs = 2000e3; S = 70e6; adcsamples = 100; GRTruth=70; framePeriod = 50e-3;numChirps =2; TargetMinIdx= TargetMinIdx; TargetMaxIdx = TargetMaxIdx; 
        [ adcData] = readDCA1000_BR_holdBreathing(fileName, adcsamples); plotEnable = viewRanges;
%         adcData = adcData(1:12000,:);
         adcData = adcData(12000:24000,:);
        [ sigEnergy, ranges,timeunit,peakVals, peaks, frPeaks] = rangeBinSelection(adcData, adcsamples, fs, S, framePeriod,numChirps, TargetMinIdx, TargetMaxIdx,  fileName, plotEnable);

    case 3 % ADC 100 Best Range - Best Range -  HR -Mean/MEdian GR Truth 88/84
        fileName ='cases/user_3.bin'; fs = 2000e3; S = 70e6; adcsamples = 100; GRTruth=84; framePeriod = 50e-3;numChirps =2; TargetMinIdx= TargetMinIdx; TargetMaxIdx = TargetMaxIdx; 
        [ adcData] = readDCA1000_BR_holdBreathing(fileName, adcsamples); plotEnable = viewRanges;
%         adcData = adcData(1:12000,:);
         adcData = adcData(12000:24000,:);
        [ sigEnergy, ranges,timeunit,peakVals, peaks, frPeaks] = rangeBinSelection(adcData, adcsamples, fs, S, framePeriod,numChirps, TargetMinIdx, TargetMaxIdx,  fileName, plotEnable);
%% ML dataset ends here--------------------------------------------------------------------------------------------------
    otherwise
        warning('No bin file selected');
            

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Investigate Phase samples for right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% range selection %%%%%%%%%%%%%%%%%%%
%%
framePeriod = framePeriod;
frPeaks = frPeaks;
phase_max_peaks = angle(peakVals);
unwrp_phase_max_peaks = unwraphase(phase_max_peaks');
% phase_max_peaksC = angle(peakValsC);
% unwrp_phase_max_peaksC = unwraphase(phase_max_peaksC');
sigEnergy = sigEnergy;
%% Scaling phasesamples to represent chest displacement
wavelength_mm = 3.9;
scale_factor = wavelength_mm/(4*pi);
phaseDiff = unwrp_phase_max_peaks;
% phaseDiffC = unwrp_phase_max_peaksC;

up =2:length(phaseDiff);
length(up);

for up =1:length(phaseDiff) 
    unwrap_phase_mm(up) = unwrp_phase_max_peaks(up)*scale_factor;
end

for up =2:length(phaseDiff) 
        ph_diff(up) = phaseDiff(up) - phaseDiff(up-1);    
end

end
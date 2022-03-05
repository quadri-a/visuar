function [pPeakLocs, pPeakVals, num_peaks] = find_peaks(filter_HeartOut)
numPeaks = 0;
for dd = 2:length(filter_HeartOut)-1
   
    if (filter_HeartOut(dd) > filter_HeartOut(dd-1)) && (filter_HeartOut(dd) > filter_HeartOut(dd+1))
        numPeaks = numPeaks +1 ;%numPeaks = numPeaks +1;
        pPeakLocs(numPeaks) = dd;
        pPeakVals(numPeaks) = filter_HeartOut(dd);%filter_HeartOut(numPeaks);
        
    end
    
end

num_peaks = numPeaks;
pPeakLocs = pPeakLocs;
pPeakVals = pPeakVals;
end
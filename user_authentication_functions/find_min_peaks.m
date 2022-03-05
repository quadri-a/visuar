function [minPeakLocs, minPeakVals, num_min_peaks] = find_min_peaks(filter_HeartOut)
numPeaks = 1;
for dd = 2:length(filter_HeartOut)-1
   
    if (filter_HeartOut(dd) < filter_HeartOut(dd-1)) && (filter_HeartOut(dd) < filter_HeartOut(dd+1))
        numPeaks = numPeaks +1 ;%numPeaks = numPeaks +1;
        pPeakLocs(numPeaks) = dd;
        pPeakVals(numPeaks) = filter_HeartOut(dd);%filter_HeartOut(numPeaks);
        
    end
    
end

num_min_peaks = numPeaks;
minPeakLocs = pPeakLocs;
minPeakVals = pPeakVals;
end
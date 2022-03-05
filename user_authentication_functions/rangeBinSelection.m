function [sigEnergy, ranges, timeunit, peakVals, peaks, frPeaks] = rangeBinSelection(adcData, adcsamples, fs, S, framePeriod,numChirps, TargetMinIdx, TargetMaxIdx, fileName, plotEnable)
fc_c = [];
[adcrw adccol] = size(adcData);
    for y = 1:adccol
        
        % Compute Range for N ADC Samples
        fc = (y/adccol)*fs;
        fc_c = [fc_c fc];
        ranges(1,y) = (3e8 * fc* 1e-6)/(2 * S);
    end  
%     figure;
    for k = 1:adcrw
       timeunit(k) = (k*framePeriod)/numChirps;
        rangeFFT(k,:) = fft(adcData(k,:)); 
        
        sigEnergy(k,:) = log(real(rangeFFT(k,:)).^2 + imag(rangeFFT(k,:)).^2);%abs(rangeFFT(k,:)).^2;
        
        TargetMinIdx = TargetMinIdx;
        TargetMaxIdx = TargetMaxIdx;
        
        % % DBFS convertion for range - power plot------------
        for gg = 1:adccol
          Pout1(k,gg) = (20*log10((abs(rangeFFT(k,gg))^2)/adccol) - ( 20*log10(2^(16-1)/1 ) + 20*log10(sum(abs(rangeFFT(k,gg)))) - 20*log10(sqrt(2))));
          if Pout1(k,gg)>-40
                pout(k,gg) = Pout1(k,gg);
          end
        end
        [rgmaxVal rgMax] = max(ranges);
        if plotEnable == 1 
          plot(ranges, Pout1(k,:));
          axis([0 rgmaxVal -120 0]);%axis([0 5 -120 0]);
          xlabel('Distance (m)');
          grid on;
          title(['Chirp, ' , num2str(k), ' of file ', fileName], 'Interpreter', 'none');
          pause(0.01);
        end
         
        [maxVal maxPeak] = max(Pout1(k,TargetMinIdx:TargetMaxIdx));
        peaks(k) = maxPeak+(TargetMinIdx-1);
        peakVals(k) = rangeFFT(k,(maxPeak+(TargetMinIdx-1)));%adcData(k,(maxPeak+(TargetMinIdx-1)));
        
%         [maxValC maxPeakC] = max(pout(k,TargetMinIdx:TargetMaxIdx));
%         peaksC(k) = maxPeakC+(TargetMinIdx-1);
%         peakValsC(k) = rangeFFT(k,(maxPeakC+(TargetMinIdx-1)));
        
    end
     
    frPeaks = mode(peaks);
    peakVals = peakVals;
    peaks = peaks;
    
%     frPeaksC = mode(peaksC);
%     peakValsC = peakValsC;
%     peaksC = peaksC;
    
end
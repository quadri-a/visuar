function further_processing(caseNum, SelectedRange)


caseNum = caseNum; 
viewRanges = 0;
TargetMinIdx = SelectedRange;
TargetMaxIdx = SelectedRange;
HeartFreq = 1.5;
[peaks, numChirps, GRTruth, sigEnergy, frPeaks, unwrap_phase_mm, fileName, framePeriod, timeunit, unwrp_phase_max_peaks, ph_diff, ranges] = select_Case(caseNum, viewRanges, TargetMinIdx, TargetMaxIdx);%select_sondadmomCase(caseNum, viewRanges, TargetMinIdx, TargetMaxIdx);

figure;
subplot(2,1,1);
plot(timeunit, unwrp_phase_max_peaks', '.-');xlabel('Time (s)');ylabel('Phase');grid on;title('Unwrapped Phase ');title(['Unwrapped Phase of range bin - ', num2str(frPeaks), ' at distance ', num2str(ranges(frPeaks)), ' m']);

subplot(2,1,2);plot( ph_diff, '.-');xlabel('Samples (Chirps for 5 min)');ylabel('Phase differences');grid on;title('Phase Difference between successive samples after unwrapping');

figure;
pspectrum( unwrp_phase_max_peaks',(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('Unwrapped phase before filtering');
figure;
pspectrum( ph_diff,(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('Ph Diff before filtering');


% EMD based filtering on unwrapped phase
[imfs, residual, info] = emd( unwrp_phase_max_peaks','MaxNumIMF',5,'Display',10); %emd( ph_diff,'MaxNumIMF',5,'Display',10);

unwrp_filt_phase = unwrp_phase_max_peaks' - imfs(:,1)'  - residual'; %- imfs(:,4)' - imfs(:,5)'

imf_fft1 = abs(fft(imfs(:,1)));
imf_fft2 = abs(fft(imfs(:,2)));
imf_fft3 = abs(fft(imfs(:,3)));
imf_fft4 = abs(fft(imfs(:,4)));
imf_fft5 = abs(fft(imfs(:,5)));
unwrp_filt_ph_fft = abs(fft(unwrp_filt_phase));

figure;
subplot(3,2,1);
plot(imf_fft1); title('IMF 1 Unwrapped Phase');
subplot(3,2,2);
plot(imf_fft2); title('IMF 2');
subplot(3,2,3);
plot(imf_fft3); title('IMF 3');
subplot(3,2,4);
plot(imf_fft4); title('IMF 4');
subplot(3,2,5);
plot(imf_fft5); title('IMF 5');
subplot(3,2,6);
plot(unwrp_filt_ph_fft); title('FFT Unwrapped Phase Removed IMF 1 2 Res');

figure;
plot(unwrp_filt_phase);

figure;
pspectrum( unwrp_filt_phase,(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('Unwrapped phase after filtering');

% Phase diff on the EMD filtered unwrapped phase
uwFFT = abs(fft(unwrp_filt_phase));
figure;plot(uwFFT);title('FFT - Unwrapped phase Filt(removed IMF 1 & 2 and residual)')

for up =2:length(unwrp_filt_phase) 
        ph_diffU(up) = unwrp_filt_phase(up) - unwrp_filt_phase(up-1);    
end
figure;plot(unwrp_filt_phase, 'r-'); hold on; plot(ph_diffU, 'b-');

fft_phDiff = abs(fft(ph_diffU));
figure;plot(fft_phDiff);title('FFT on PH DIFF after filtering unwrpd phase with EMD');

figure;
pspectrum(  ph_diffU,(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('Ph Diff from filtered unwrpd phase');

% Perform EMD on phDiffU
[imfs1, residual1, info1] = emd( ph_diffU,'MaxNumIMF',5,'Display',10);

figure;
subplot(3,2,1);
plot(abs(fft(imfs1(:,1))));title('PhDiff IMF 1 FFT');

subplot(3,2,2);
plot(abs(fft(imfs1(:,2))));title('PhDiff IMF 2 FFT');

subplot(3,2,3);
plot(abs(fft(imfs1(:,3))));title('PhDiff IMF 3 FFT');

subplot(3,2,4);
plot(abs(fft(imfs1(:,4))));title('PhDiff IMF 4 FFT');

subplot(3,2,5);
plot(abs(fft(imfs1(:,5))));title('PhDiff IMF 5 FFT');


%% Heart Rate Filter
 ph_diff_hr = imfs1(:,2);%ph_diffU;
% filter_HeartOut = sig;
IIR_FILTER_HEART_NUM_STAGES = 4;
IIR_FILTER_COEFS_SECOND_ORDER = 6;
HEART_WFM_IIR_FILTER_TAPS_LENGTH = (IIR_FILTER_HEART_NUM_STAGES *IIR_FILTER_COEFS_SECOND_ORDER) +2;

pDelayHeart = zeros(1,HEART_WFM_IIR_FILTER_TAPS_LENGTH);

% For 0.8 to 2 Hz
pFilterCoefsHeart = [1.0000, 0, -1.0000, 1.0000, -1.4907, 0.8216, 1.0000, 0, -1.0000, 1.0000, -1.8530, 0.9166,1.0000, 0, -1.0000, 1.0000, -1.6588, 0.7512, 1.0000, 0, -1.0000, 1.0000, -1.4616, 0.6555];
pScaleValsHeart = [0.1754, 0.1754, 0.1619, 0.1619, 1.0000];

for hh = 1:length(ph_diff_hr)
    filter_HeartOut(hh) = filter_IIR_biquad_heart(ph_diff_hr(hh),pFilterCoefsHeart, pScaleValsHeart, pDelayHeart, IIR_FILTER_HEART_NUM_STAGES);
end
figure;
pspectrum( filter_HeartOut,(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('After iir');

figure;plot(filter_HeartOut, 'r-'); hold on; plot(ph_diff_hr, 'b-');
figure;plot(abs(fft(filter_HeartOut)));title('EMD IIR Filtered Ph DiffU');
% % %% Extract cardiac cycles
% % % Extract cardiac cycles
heartPeaks = imfs1(:,2);%filter_HeartOut;% sig;% FOR 2 LOOP FRAMES imfs(:,2)';
plotcyc = 0;

[cycle_features, cardiacPeaks, diastolic_p1, diastolic_p2, systolic_p1, Maxloc_1, Maxloc_2 ,Minloc_2, Minloc_1, Minloc_3]= cycles_emd_iir_percent(filter_HeartOut', framePeriod, numChirps, plotcyc, HeartFreq);
%cycles_observation_2(heartPeaks, framePeriod, numChirps, plotcyc);
% sel_features = cycle_features(1:40,:);
% dlmwrite('features_24_emdU_iirPhDiff_phDiff_sondad_rb2220.csv',sel_features,'-append');
% dlmwrite('unwrappedEMD_signal_48_rb18.csv',unwrp_filt_phase,'-append');
% % % % % figure;plot(ph_diffU, 'r-');hold on; plot(filter_HeartOut, 'b-');

figure;plot(filter_HeartOut, 'b-');




% Usual EMD on phdiff -> then IIR filter then doppler and then cycles
% extraction with cycle_observation_2 - Best case scenario was 46 sig - 19
% HeartFreq - 1.38
% RB 47 - 20 RB  HeartFreq - 1.15  and 48 - 18 RB heartFreq 1.52
% % caseNum = 48; 
% % viewRanges = 0;
% % TargetMinIdx = 17;
% % TargetMaxIdx = 17;
% % 
% % [peaks, numChirps, GRTruth, sigEnergy, frPeaks, unwrap_phase_mm, fileName, framePeriod, timeunit, unwrp_phase_max_peaks, ph_diff, ranges] = select_sondadmomCase(caseNum, viewRanges, TargetMinIdx, TargetMaxIdx);
% % 
% % figure;
% % subplot(2,1,1);
% % plot(timeunit, unwrp_phase_max_peaks', '.-');xlabel('Time (s)');ylabel('Phase');grid on;title('Unwrapped Phase ');title(['Unwrapped Phase of range bin - ', num2str(frPeaks), ' at distance ', num2str(ranges(frPeaks)), ' m']);
% % 
% % subplot(2,1,2);plot( ph_diff, '.-');xlabel('Samples (Chirps for 5 min)');ylabel('Phase differences');grid on;title('Phase Difference between successive samples after unwrapping');
% % 
% % figure;
% % pspectrum( ph_diff,(1/framePeriod),'spectrogram','TimeResolution',(1/framePeriod));
% % title('Unwrapped phase before filtering');
% % 
% % 
% % % EMD based filtering on unwrapped phase
% % [imfs, residual, info] = emd( ph_diff,'MaxNumIMF',5,'Display',10); %emd( ph_diff,'MaxNumIMF',5,'Display',10);
% % 
% % 
% % sig = imfs(:,3)';%FOR 46 47 48 it is imfs(:,3)' for 45 43 41 it is imfs(:,2)%ph_diff - imfs(:,1)' ;% - imfs(:,3)' - imfs(:,4)' - imfs(:,5)' - residual';
% % figure;
% % pspectrum(sig,(1/framePeriod),'spectrogram','TimeResolution',(1/framePeriod));
% % title('EMD Filtered sig');
% % 
% % figure;subplot(6,1,1);plot(imfs(:,1));title('IMF 1');
% % subplot(6,1,2);plot(imfs(:,2));title('IMF 2');
% % subplot(6,1,3);plot(imfs(:,3));title('IMF 3');
% % subplot(6,1,4);plot(imfs(:,4));title('IMF 4');
% % subplot(6,1,5);plot(imfs(:,5));title('IMF 5');
% % subplot(6,1,6);plot(residual);title('IMF 6');
% % 
% % x1 = abs(fft(imfs(:,1)));
% % figure;plot(x1);title('FFT of IMF1');
% % 
% % x2 = abs(fft(imfs(:,2)));
% % figure;plot(x2);title('FFT of IMF2');
% % 
% % x3 = abs(fft(imfs(:,3)));
% % figure;plot(x3);title('FFT of IMF3');
% % 
% % x4 = abs(fft(imfs(:,4)));
% % figure;plot(x4);title('FFT of IMF4');
% % 
% % x5 = abs(fft(imfs(:,5)));
% % figure;plot(x5);title('FFT of IMF5');
% % 
% % x6 = abs(fft(sig));
% % figure;plot(x6);title('FFT of SIG');
% % 
% % % x7 = abs(fft(dopFilter_heartOut));
% % % figure;plot(x7);title('FFT of Dopfiltered');
% % 
% % %% Heart Rate Filter
% %  ph_diff_hr = imfs(:,3)';
% % % filter_HeartOut = sig;
% % IIR_FILTER_HEART_NUM_STAGES = 4;
% % IIR_FILTER_COEFS_SECOND_ORDER = 6;
% % HEART_WFM_IIR_FILTER_TAPS_LENGTH = (IIR_FILTER_HEART_NUM_STAGES *IIR_FILTER_COEFS_SECOND_ORDER) +2;
% % 
% % pDelayHeart = zeros(1,HEART_WFM_IIR_FILTER_TAPS_LENGTH);
% % 
% % % For 0.8 to 2 Hz
% % pFilterCoefsHeart = [1.0000, 0, -1.0000, 1.0000, -1.4907, 0.8216, 1.0000, 0, -1.0000, 1.0000, -1.8530, 0.9166,1.0000, 0, -1.0000, 1.0000, -1.6588, 0.7512, 1.0000, 0, -1.0000, 1.0000, -1.4616, 0.6555];
% % pScaleValsHeart = [0.1754, 0.1754, 0.1619, 0.1619, 1.0000];
% % 
% % for hh = 1:length(ph_diff_hr)
% %     filter_HeartOut(hh) = filter_IIR_biquad_heart(ph_diff_hr(hh),pFilterCoefsHeart, pScaleValsHeart, pDelayHeart, IIR_FILTER_HEART_NUM_STAGES);
% % end
% % figure;
% % pspectrum( filter_HeartOut,(1/framePeriod),'spectrogram','TimeResolution',(1/framePeriod));
% % title('After iir');
% % 
% % 
% % dopFilter_heartOut = [];
% % dopplerWinSize = 4;%16;
% % dopplerWin = floor(length(filter_HeartOut)/dopplerWinSize);
% % dopplerCoefs = [0.0800, 0.0894, 0.1173, 0.1624];%,0.2231, 0.2967, 0.3802, 0.4703, 0.5633, 0.6553, 0.7426, 0.8216, 0.8890, 0.9422, 0.9789, 0.9976];
% % for dp = 1:dopplerWin
% %         dopFilter_heartOut = [ dopFilter_heartOut (filter_HeartOut((dp-1)*dopplerWinSize+1:dopplerWinSize*dp).*dopplerCoefs) ];
% % end
% % figure;
% % pspectrum( dopFilter_heartOut,(1/framePeriod),'spectrogram','TimeResolution',(1/framePeriod));
% % title('dopfilter');
% % 
% % 
% % %% Extract cardiac cycles
% % % Extract cardiac cycles
% % heartPeaks = dopFilter_heartOut; %filter_HeartOut;% sig;% FOR 2 LOOP FRAMES imfs(:,2)';
% % plotcyc = 0;
% % 
% % [cycle_features, cardiacPeaks, diastolic_p1, diastolic_p2, systolic_p1, Maxloc_1, Maxloc_2 ,Minloc_2, Minloc_1, Minloc_3]= cycles_observation_2(heartPeaks, framePeriod, numChirps, plotcyc);
% % sel_features = cycle_features(1:200,:);
% % dlmwrite('features_24_emd_iir_phdiff_sondadmom_rb221817.csv',sel_features,'-append');
% % % heartPeaks1 = imfs(:,3)';% FOR 2 LOOP FRAMES imfs(:,2)';
% % % plotcyc = 0;
% % % 
% % % [cycle_features1, cardiacPeaks1, diastolic_p11, diastolic_p21, systolic_p11, Maxloc_11, Maxloc_21 ,Minloc_21, Minloc_11, Minloc_31]= cycles_observation_2(heartPeaks1, framePeriod, numChirps, plotcyc);
% % 
% % % all_features = [];
% % % all_features = [cycle_features cycle_features1];
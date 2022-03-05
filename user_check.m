clear all;
clc;
close all;

dirc = addpath('user_authentication_functions');


disp('03 Users in Database: Start Signal Processing and Generation of Cardiac Profile for - ');
prompt = input('User No. = ');
caseNum = prompt; 
disp(['Started processing data of User no. ', num2str(prompt), '. . . . . . . . . . . . . . .']);
viewRanges = 0;
framePeriod = 50e-3;
showL = 8000;
phase_fft_size = 12000;
sampling_Freq = 2/framePeriod;
freqIncrement = sampling_Freq/phase_fft_size;
heart_startFreq = 0.8; % CHANGE HERE ______________________________
heart_endFreq =2;
br_startFreq = 0.1;
br_endFreq = 0.5;
heart_startFreq_Idx = floor(heart_startFreq/freqIncrement);
heart_endFreq_Idx = ceil(heart_endFreq/freqIncrement);
br_startFreq_Idx = floor(br_startFreq/freqIncrement); 
br_endFreq_Idx = ceil(br_endFreq/freqIncrement);

if prompt == 1
    TargetMinIdx = 17;
    TargetMaxIdx = 17;
    HeartFreq = 1.51;
end
if prompt == 2
    TargetMinIdx = 18;
    TargetMaxIdx = 18;
    HeartFreq = 1.2;
end
if prompt == 3
    TargetMinIdx = 18;
    TargetMaxIdx = 18;
    HeartFreq = 1.51;
end
disp('Processing ./.\./.\./.\./.\./.\./.\./.\');
[peaks, numChirps, GRTruth, sigEnergy, frPeaks, unwrap_phase_mm, fileName, framePeriod, timeunit, unwrp_phase_max_peaks, ph_diff, ranges] = select_Case(caseNum, viewRanges, TargetMinIdx, TargetMaxIdx);%select_sondadmomCase(caseNum, viewRanges, TargetMinIdx, TargetMaxIdx);

uw_ph_fft(1,:) = abs(fft(ph_diff));
[brfft brIdx] = max(uw_ph_fft(1,br_startFreq_Idx:br_endFreq_Idx));
[mxfft fftidx] = max(uw_ph_fft(1,heart_startFreq_Idx:heart_endFreq_Idx));
bridx = (brIdx + br_startFreq_Idx) - 1;
hridx = (fftidx + heart_startFreq_Idx) -1;
br = ((bridx/phase_fft_size)*sampling_Freq)*60;
hr = ((hridx/phase_fft_size)*sampling_Freq)*60;

figure(1);
subplot(2,1,1);
plot(timeunit, unwrp_phase_max_peaks', '.-');xlabel('Time (s)');ylabel('Phase');grid on;title('Signal Phase ');title(['1. Signal from radar range bin - ', num2str(frPeaks), ' at distance ', num2str(ranges(frPeaks)), ' m']);

subplot(2,1,2);plot( ph_diff, '.-');xlabel('Samples');ylabel('Phase differences');grid on;title({['2. Body Movement Related Motions Supressed - Cardiac Motion Retained '],['Estimated Heart Beat Rate: ',num2str(floor(hr)), ' bpm']});

%% Filtering to remove body motion entirely
disp('Starting the process to filter out any signal other than cardiac signal of user. . . . . . . . . . . . . . . ');
disp('Processing ./.\./.\./.\./.\./.\./.\./.\');
% EMD based filtering on unwrapped phase
[imfs, residual, info] = emd( unwrp_phase_max_peaks','MaxNumIMF',5, 'Display',0);

unwrp_filt_phase = unwrp_phase_max_peaks'  - imfs(:,1)' - residual' ; 
%  dlmwrite('unwrappedphEMD_signal_48_rb18_2ndHalf_imf1RESremov.csv',unwrp_filt_phase,'-append');
figure(2);
subplot(4,1,1);plot(unwrp_filt_phase(:,showL:end));title('3. Filtering Original Signal: Before Decomposition');
subplot(4,1,2);plot(imfs(showL:end,1));title('4. After Empirical Mode Decomposition: Filtered Signal 1');
subplot(4,1,3);plot(imfs(showL:end,2));title('5. Filtered Signal 2');
subplot(4,1,4);plot(imfs(showL:end,3));title('6. Filtered Signal 3');

f = figure(3);
f.Position = [200, 400, 1200, 400];
subplot(1,2,1);
pspectrum( unwrp_phase_max_peaks(showL:end,:)',(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('7. Before Filtering: Corrupted by Body Motion');
subplot(1,2,2);
pspectrum( unwrp_filt_phase(:,showL:end),(2/framePeriod),'spectrogram','TimeResolution',(2/framePeriod));
title('8. Filtered Signal 1: Without any Body Motion');

%% NEXT PART STARTS WITH SELECTED DATASETS

if prompt == 1
    user_data = csvread('filtered_sigs/user_1_sig.csv');
end
if prompt == 2
    user_data = csvread('filtered_sigs/user_2_sig.csv');
end
if prompt == 3
    user_data = csvread('filtered_sigs/user_3_sig.csv');
end

framePeriod= 50e-3;
numChirps = 2;
%% Heart Rate Filter
 ph_diff_hr = user_data;
%% Phase difference of the EMD filtered signal
for up =2:length(ph_diff_hr) 
        ph_diffU_uF(up) = ph_diff_hr(up) - ph_diff_hr(up-1);    
end

IIR_FILTER_HEART_NUM_STAGES = 4;
IIR_FILTER_COEFS_SECOND_ORDER = 6;
HEART_WFM_IIR_FILTER_TAPS_LENGTH = (IIR_FILTER_HEART_NUM_STAGES *IIR_FILTER_COEFS_SECOND_ORDER) +2;

pDelayHeart = zeros(1,HEART_WFM_IIR_FILTER_TAPS_LENGTH);

% For 0.8 to 2 Hz
pFilterCoefsHeart = [1.0000, 0, -1.0000, 1.0000, -1.4907, 0.8216, 1.0000, 0, -1.0000, 1.0000, -1.8530, 0.9166,1.0000, 0, -1.0000, 1.0000, -1.6588, 0.7512, 1.0000, 0, -1.0000, 1.0000, -1.4616, 0.6555];
pScaleValsHeart = [0.1754, 0.1754, 0.1619, 0.1619, 1.0000];

for hh = 1:length(ph_diff_hr)
%     filter_HeartOut(hh) = filter_IIR_biquad_heart(ph_diff_hr(hh),pFilterCoefsHeart, pScaleValsHeart, pDelayHeart, IIR_FILTER_HEART_NUM_STAGES);
    filter_HeartOut(hh) = filter_IIR_biquad_heart(ph_diffU_uF(hh),pFilterCoefsHeart, pScaleValsHeart, pDelayHeart, IIR_FILTER_HEART_NUM_STAGES);
end

figure(4);plot(filter_HeartOut(:,showL:end));grid on;xlabel('Samples');ylabel('Amplitude');title(['Cardiac Motion of User ', num2str(prompt)]);
disp('Filtering process ended ++++++++++++++++++++++++++++++++++++++++++++');
%% Doppler Coefficients 

dopFilter_heartOut = [];
dopplerWinSize = 4;%16;
dopplerWin = floor(length(filter_HeartOut)/dopplerWinSize);
dopplerCoefs = [0.0800, 0.0894, 0.1173, 0.1624];%,0.2231, 0.2967, 0.3802, 0.4703, 0.5633, 0.6553, 0.7426, 0.8216, 0.8890, 0.9422, 0.9789, 0.9976];
for dp = 1:dopplerWin
    dopFilter_heartOut = [ dopFilter_heartOut (filter_HeartOut((dp-1)*dopplerWinSize+1:dopplerWinSize*dp).*dopplerCoefs) ];
%         dopFilter_heartOut = [ dopFilter_heartOut (filter_HeartOut((dp-1)*dopplerWinSize+1:dopplerWinSize*dp).*dopplerCoefs) ];
end
% % figure;plot(filter_HeartOut, 'r-'); hold on; plot(dopFilter_heartOut, 'b-');title('Doopler Filtered');%

%% Phase difference of the iir filtered signal
for up =2:length(filter_HeartOut) 
        ph_diffU(up) = filter_HeartOut(up) - filter_HeartOut(up-1);    
end
% % figure;plot(filter_HeartOut, 'r-'); hold on; plot(ph_diffU, 'b-');

disp('Starting the process to detect cardiac cycles and extract unique features to create a cardiac profile. . . . . . . . . . . . . . . ');
disp('Processing ./.\./.\./.\./.\./.\./.\./.\');
%% Cycle Extraction from cleaned unwrapped phase
heartPeaks = dopFilter_heartOut;
plotcyc = 0;
HeartFreq = 1.51;

phase_fft_size = length(heartPeaks);

heartPeaks_fft = abs(fft(heartPeaks));
% % figure;plot(heartPeaks_fft);
sampling_Freq = 2/framePeriod;
freqIncrement = sampling_Freq/phase_fft_size;
heart_startFreq = 0.8; % CHANGE HERE ______________________________
heart_endFreq =2;
br_startFreq = 0.1;
br_endFreq = 0.5;
heart_startFreq_Idx = floor(heart_startFreq/freqIncrement);
heart_endFreq_Idx = ceil(heart_endFreq/freqIncrement);
br_startFreq_Idx = floor(br_startFreq/freqIncrement); 
br_endFreq_Idx = ceil(br_endFreq/freqIncrement);

[blockL blockIdx] = max(heartPeaks_fft(:,br_startFreq_Idx:br_endFreq_Idx ));
blockWidth = blockIdx + br_startFreq_Idx - 1;

Hrfreq = (blockWidth/phase_fft_size)*40;
period = 1/Hrfreq;
chirpBlock = floor((period*numChirps)/framePeriod);

blocks = chirpBlock; 
window = floor(length(heartPeaks)/blocks) - 1;
min_locs = [];
for t = 1:window
    [min_vals min_idx] = min(heartPeaks(((t -1)*blocks)+1 : (t*blocks)+3));
        min_locs(:,t) = min_idx + ((t -1)*blocks);


end
showHP = showL+3000;
min_ps = min_locs(min_locs> showHP);
hs = heartPeaks(:,showHP:end);
figure(5);plot(hs);hold on;plot((min_ps- showHP+1), hs((min_ps- showHP+1)), 'ro');grid on;xlabel('Samples');ylabel('Amplitude');title('Detected Consecutive Cardiac Cycles');

nfid = 10;
mxlocs = [];
mnlocs = [];
features = [];
height_pk_dip_pk = [];
norm_height_min1_PkdipPk = [];
countNFIDmax = 0;
for h =1:length(min_locs)-1
    if (h == 1) && ((min_locs(h) - 1) > blocks/2 )
        
        [mxploc1, pval1, num_p1] = find_peaks(heartPeaks(1:min_locs(h)) );
        [mnploc1, pval1, num_p1] = find_min_peaks(heartPeaks(1:(min_locs(h)-1)) );
        mnploc1 = mnploc1(mnploc1>0);
        
        cyc_strt = 1;
        % Find the max peak in the whole cycle and pick all the dips and
        % peaks after max peak
        full_cyc = (min_locs(h));
        new_mx = mxploc1(mxploc1< full_cyc );
        
        [mxlocV,mxlocSidx] = max(heartPeaks(new_mx));
        sel_mxlocs = mxploc1(mxploc1>=new_mx(mxlocSidx));
        sel_mnplocs = mnploc1(mnploc1>=sel_mxlocs(1));
            
        % Choose only selected amount of peaks and dips as defined by
        % nfid
        peaks_mix = [sel_mnplocs sel_mxlocs];
        peaks_sort = sort(peaks_mix, 'ascend');
        % If there aren't enough peaks and dips in the selected cycle
        if (length(peaks_sort) < nfid)

            countNFIDmax = countNFIDmax + 1;
% %             if (countNFIDmax < 2)
% %                 figure;plot(heartPeaks(1:min_locs(h))); hold on; plot(sel_mxlocs, heartPeaks(sel_mxlocs), 'r*');plot(sel_mnplocs, heartPeaks(sel_mnplocs), 'b*');
% %             end 
        end
        
        if (nfid < length(peaks_sort))
            peaks = peaks_sort(1:nfid);%peaks_sort(end-nfid:end);
            pulse_width = ((min_locs(h) - cyc_strt)*framePeriod)/numChirps;
            height_min1 = heartPeaks(peaks(1)) - heartPeaks(cyc_strt);
            % Time features - each peaks comapared to cyc start location
            for ft = 1:length(peaks)
                time_cycS(:,ft) = ((peaks(ft) - cyc_strt)*framePeriod)/numChirps;
                height_cycS(:,ft) = abs(heartPeaks(peaks(ft)) - heartPeaks(cyc_strt));            
                norm_timeCYCS(:,ft) = time_cycS(:,ft)/pulse_width;
                norm_height_cycS(:,ft) = height_cycS(:,ft)/height_min1;
            end
            for ft = 1:length(peaks)-1    
                time_pk_pk(:, ft) = ((peaks(ft + 1) - peaks(ft))*framePeriod)/numChirps;
                height_pk_pk(:, ft) = abs(heartPeaks(peaks(ft+1)) - heartPeaks(peaks(ft)));
                norm_timePKPK(:,ft) = time_pk_pk(:,ft)/pulse_width;
            end 
                        
            pk_size =(length(peaks)/2);
            for fts=1:pk_size
                PKdip(:,fts) = heartPeaks(peaks((2*fts) -1)) - heartPeaks(peaks(2*fts));  
                norm_PKdip(:,fts) = PKdip(:,fts)/height_min1;
            end
            for fts=1:pk_size -1
                pekloc = peaks(end - ((2*fts)-1));
                diploc = peaks(end - (2*fts));
                pek = heartPeaks(peaks(end - ((2*fts)-1)));
                dip = heartPeaks(peaks(end - (2*fts)));
                dipPK(:, fts) = pek - dip;
                norm_dipPK(:,fts) = dipPK(:,fts)/height_min1;
            end
            
            height_min2 = heartPeaks(peaks(end)) - heartPeaks(min_locs(h));
            time_min1 = ((peaks(1) - cyc_strt)*framePeriod)/numChirps;
            time_min2 = ((min_locs(h) - peaks(end))*framePeriod)/numChirps;
            norm_timemin1 = time_min1/pulse_width;
            norm_timemin2 = time_min2/pulse_width;
                
            height_pk_dip_pk = [height_min1  PKdip dipPK height_min2];
            time_min_pk_min = [time_min1 time_pk_pk time_min2];
            norm_timemin_pk_min = [norm_timemin1 norm_timePKPK norm_timemin2];
            norm_height = height_min2/height_min1;
            norm_height_min1_PkdipPk = [norm_PKdip norm_dipPK norm_height];   

            features = [features; time_cycS height_cycS  norm_timeCYCS time_min_pk_min height_pk_dip_pk norm_timemin_pk_min 0 height_pk_pk 0 norm_height_min1_PkdipPk 0 norm_height_cycS];
        end
    end    
        if (h ~=1)
            minlocDiffs(h) = min_locs(h) - min_locs(h-1);
            if h==60
                [mxv1, pKval12, num_pK12] = find_peaks(heartPeaks( min_locs(h-4): min_locs(h) ) );
                 [mnv1, pKva21, num_pK21] = find_min_peaks(heartPeaks( (min_locs(h-4)): (min_locs(h)) ) );
                f2 = figure(6);
                f2.Position = [200, 400, 800, 400];
                plot(heartPeaks(min_locs(h-4):min_locs(h)));hold on;plot(mxv1, heartPeaks(mxv1+min_locs(h-4)-1), 'ro');plot(mnv1, heartPeaks(mnv1+min_locs(h-4)-1), 'bo');grid on;xlabel('Samples');ylabel('Signal Amplitude');title(['Cardiac Profile of User ', num2str(prompt) ,' created with features extracted from cardiac cycles']);
%                 plot(mnv1, heartPeaks(mnv1), 'bo');
            end
        end
        if (h ~=1) && ( ( min_locs(h) - min_locs(h-1) > (blocks - 20)  ) && (min_locs(h) - min_locs(h-1) < (blocks + 20)  ) ) %&&((min_locs(h) - min_locs(h-1)) > blocks/2 )%
            
            [mxploc1, pval1, num_p1] = find_peaks(heartPeaks( (min_locs(h-1)+1): (min_locs(h)-1) ) );
            mxploc1 = mxploc1 + (min_locs(h-1));
            [mnploc1, pval1, num_p1] = find_min_peaks(heartPeaks( (min_locs(h-1)+1): (min_locs(h)-1) ) );
            mnploc1 = mnploc1 + (min_locs(h-1));
            mnploc1 = mnploc1(mnploc1>0);
           
            cyc_strt = min_locs(h-1);
            % Find the max peak in the whole cycle and pick all the dips and
            % peaks after max peak
            full_cyc = (min_locs(h));
            new_mx = mxploc1(mxploc1< full_cyc );
            
            [mxlocV,mxlocSidx] = max(heartPeaks(new_mx));
            sel_mxlocs = mxploc1(mxploc1>=new_mx(mxlocSidx));
            sel_mnplocs = mnploc1(mnploc1>=sel_mxlocs(1));
            
            % Choose only selected amount of peaks and dips as defined by
            % nfid
            peaks_mix = [sel_mnplocs sel_mxlocs];
            peaks_sort = sort(peaks_mix, 'ascend');
            if (length(peaks_sort) < nfid)
%                 min_locs(h)
%                 mxploc1
%                 new_mx(mxlocSidx)
                countNFIDmax = countNFIDmax + 1;
% %                 if (countNFIDmax < 2)
% %                     figure;plot(heartPeaks(1:min_locs(h))); hold on; plot(sel_mxlocs, heartPeaks(sel_mxlocs), 'r*');plot(sel_mnplocs, heartPeaks(sel_mnplocs), 'b*');
% %                 end 
            end
            
            if (nfid < length(peaks_sort))
                peaks = peaks_sort(1:nfid);%peaks_sort((end-nfid):end);
    %             figure;plot(heartPeaks(1:min_locs(h))); hold on;plot(min_locs(h-1), heartPeaks(min_locs(h-1)), 'bo');plot(min_locs(h), heartPeaks(min_locs(h)), 'bo'); plot(peaks, heartPeaks(peaks), 'r*');

                pulse_width = ((min_locs(h) - cyc_strt)*framePeriod)/numChirps;
                height_min1 = heartPeaks(peaks(1)) - heartPeaks(cyc_strt);
                %Time features - each peaks compared to cyc start location 
                for ft = 1:length(peaks)

                    time_cycS(:,ft) = ((peaks(ft) - cyc_strt)*framePeriod)/numChirps;
                    height_cycS(:,ft) = abs(heartPeaks(peaks(ft)) - heartPeaks(cyc_strt));
                    norm_height_cycS(:,ft) = height_cycS(:,ft)/height_min1;
                    norm_timeCYCS(:,ft) = time_cycS(:,ft)/pulse_width;
                end
                for ft = 1:length(peaks) - 1
                    time_pk_pk(:, ft) = ((peaks(ft + 1) - peaks(ft))*framePeriod)/numChirps;
                    height_pk_pk(:, ft) = abs(heartPeaks(peaks(ft+1)) - heartPeaks(peaks(ft)));
                    norm_timePKPK(:,ft) = time_pk_pk(:,ft)/pulse_width;           
                end
                              
                pk_size =(length(peaks)/2);
                for fts=1:pk_size
                    PKdip(:,fts) = heartPeaks(peaks((2*fts) -1)) - heartPeaks(peaks(2*fts)); 
                    norm_PKdip(:,fts) = PKdip(:,fts)/height_min1;
                end
                for fts=1:pk_size -1
                    pekloc = peaks(end - ((2*fts)-1));
                    diploc = peaks(end - (2*fts));
                    pek = heartPeaks(peaks(end - ((2*fts)-1)));
                    dip = heartPeaks(peaks(end - (2*fts)));
                    dipPK(:, fts) = pek - dip;
                    norm_dipPK(:,fts) = dipPK(:,fts)/height_min1;
                end
                
                height_min2 = heartPeaks(peaks(end)) - heartPeaks(min_locs(h));
                time_min1 = ((peaks(1) - cyc_strt)*framePeriod)/numChirps;
                time_min2 = ((min_locs(h) - peaks(end))*framePeriod)/numChirps;
                norm_timemin1 = time_min1/pulse_width;
                norm_timemin2 = time_min2/pulse_width;
                
                height_pk_dip_pk = [height_min1  PKdip dipPK height_min2];
                time_min_pk_min = [time_min1 time_pk_pk time_min2];
                norm_timemin_pk_min = [norm_timemin1 norm_timePKPK norm_timemin2];
                norm_height = height_min2/height_min1;
                norm_height_min1_PkdipPk = [norm_PKdip norm_dipPK norm_height];

                features = [features; time_cycS height_cycS  norm_timeCYCS time_min_pk_min height_pk_dip_pk norm_timemin_pk_min 0 height_pk_pk 0 norm_height_min1_PkdipPk 0 norm_height_cycS];
        
            end
        end
% %         if h == length(min_locs)-1
% %               
% %                 figure(6);plot(heartPeaks(:,sel_mnplocs(1):sel_mnplocs(10)));hold on;
% %                 plot((sel_mnplocs(1)-sel_mnplocs(1)+1):(sel_mnplocs(10)-sel_mnplocs(1)+1), heartPeaks((sel_mnplocs(1):sel_mnplocs(10))), 'rx');grid on;xlabel('Cardiac Samples');ylabel('Signal Amplitude');
% %                 title(['9. Cardiac Profile of User ', num2str(prompt) ,' : Features extracted from cardiac cycles']);
% %         end
                    
%     end
end
disp('Starting the process to perform Principal Component Analysis on the extracted features . . . . . . . . . . . . . . . ');
disp('Processing ./.\./.\./.\./.\./.\./.\./.\');
%% NEXT IS THE PCA PART TO EXTRACT FEATURES ADN PREP FOR ML SVM TRAINING
%% Mixed norm_time_pk_pk and norm_Height_min1_pk_pk
c1 = 57; c2 = 60; c3 = 81; c4 = 84; % Time norm Height feats PK-PK (SVM ACCURACY 90)
p1 = 1; p2 = 8;

mom2 = csvread('pca_analysis/user_3m_pca_transform.csv');
mom = [mom2(:, c1:c2) mom2(:, c3:c4)];% mom2(:, c1:c2) mom2(:, c3:c4)];

dad2 = csvread('pca_analysis/user_2d_pca_transform.csv');
dad = [dad2(:, c1:c2) dad2(:, c3:c4)]; %dad2(:, c1:c2) dad2(:, c3:c4)];

son2 = csvread('pca_analysis/user_1s_pca_transform.csv');
son = [son2(:, c1:c2) son2(:, c3:c4)];% son2(:, c1:c2) son2(:, c3:c4)];


mom_stack = [];
for tr = 1:length(mom)
    mom_stack = [mom_stack; 1 mom(tr,:)];
end

dad_stack = [];
for tr = 1:length(dad)
    dad_stack = [dad_stack; 2 dad(tr,:)];
end

son_stack = [];
for tr = 1:length(son)
    son_stack = [son_stack; 3 son(tr,:)];
end

%% Imbalanced CLass %%%% PCA Over entire data set then using that coeff to transform
feats = [mom_stack(1:(end-5),p1:p2);dad_stack(1:(end-5),p1:p2);son_stack(1:(end-5),p1:p2)];%[mom_stack(1:30,p1:p2);dad_stack(1:30,p1:p2);son_stack(1:30,p1:p2)];%
[coeff, score, latent] = pca(feats);

mom_stackTr = mom_stack(1:(end-5),p1:p2)*coeff(:,1:8);%mom_stack(1:(end-5),p1:p2);%mom_stack(1:30,p1:p2)*coeff(:,1:8);
dad_stackTr = dad_stack(1:(end-5),p1:p2)*coeff(:,1:8);%dad_stack(1:(end-5),p1:p2);%dad_stack(1:30,p1:p2)*coeff(:,1:8);%%%
son_stackTr = son_stack(1:(end-5),p1:p2)*coeff(:,1:8);%son_stack(1:30,p1:p2)*coeff(:,1:8);%son_stack(1:(end-5),p1:p2);%%

mom_stackP = mom_stack((end-4):end,p1:p2)*coeff(:,1:8);%mom_stack(31:34,p1:p2)*coeff(:,1:8);%mom_stack((end-4):end,p1:p2);%%
dad_stackP = dad_stack((end-4):end,p1:p2)*coeff(:,1:8);%dad_stack(31:34,p1:p2)*coeff(:,1:8);%dad_stack((end-4):end,p1:p2);%%
son_stackP = son_stack((end-4):end,p1:p2)*coeff(:,1:8);%son_stack(31:34,p1:p2)*coeff(:,1:8);%son_stack((end-4):end,p1:p2);%%

f1 = 1; f2 =2; f3 =3;
f_pic = figure(7);
f_pic.Position = [200, 400, 1200, 400];
subplot(1,2,1);
% figure(7);
% % p1 = plot3( mom_stackP(:,f1), mom_stackP(:,f2), mom_stackP(:,f3), 'kx');
% % hold on;
% % p2 = plot3(dad_stackP(:,f1),dad_stackP(:,f2),dad_stackP(:,f3), 'rx');
% % p3 = plot3(son_stackP(:,f1),son_stackP(:,f2),son_stackP(:,f3),'bx');

p4 = plot3( mom_stackTr(:,f1), mom_stackTr(:,f2), mom_stackTr(:,f3), 'ko');
hold on;
p5 = plot3(dad_stackTr(:,f1),dad_stackTr(:,f2),dad_stackTr(:,f3), 'ro');
p6 = plot3(son_stackTr(:,f1),son_stackTr(:,f2),son_stackTr(:,f3),'bo');
xlabel('1. PCA 4 transformed features');ylabel('PCA 5 transformed features');zlabel('PCA 6 transformed features');title('10. First 3 Features of All Users in Database');legend('User 1', 'User 2', 'User 3');


f1 = 4; f2 =5; f3 =6;
% f1 = 10; f2 =11; f3 =12;
subplot(1,2,2);
% % p1 = plot3( mom_stackP(:,f1), mom_stackP(:,f2), mom_stackP(:,f3), 'kx');
% % hold on;
% % p2 = plot3(dad_stackP(:,f1),dad_stackP(:,f2),dad_stackP(:,f3), 'rx');
% % p3 = plot3(son_stackP(:,f1),son_stackP(:,f2),son_stackP(:,f3),'bx');

p4 = plot3( mom_stackTr(:,f1), mom_stackTr(:,f2), mom_stackTr(:,f3), 'ko');
hold on;
p5 = plot3(dad_stackTr(:,f1),dad_stackTr(:,f2),dad_stackTr(:,f3), 'ro');
p6 = plot3(son_stackTr(:,f1),son_stackTr(:,f2),son_stackTr(:,f3),'bo');
xlabel('2. PCA transformed feature 4');ylabel('PCA transformed feature 5');zlabel('PCA transformed feature 6');title(' Second 3 Features of All Users in Database');legend('User 1', 'User 2', 'User 3');


mom_data = [];
for tr = 1:length(mom_stackTr)
    mom_data = [mom_data; 1 mom_stackTr(tr,:)];
end

dad_data = [];
for tr = 1:length(dad_stackTr)
    dad_data = [dad_data; 2 dad_stackTr(tr,:)];
end

son_data = [];
for tr = 1:length(son_stackTr)
    son_data = [son_data; 3 son_stackTr(tr,:)];
end

mom_dataP = [];
[rw cl] = size(mom_stackP);
for te = 1:rw
    mom_dataP = [mom_dataP; 1 mom_stackP(te,:)];
end

dad_dataP = [];
[rw cl] = size(dad_stackP);
for te = 1:rw
    dad_dataP = [dad_dataP; 2 dad_stackP(te,:)];
end

son_dataP = [];
[rw cl] = size(son_stackP);
for te = 1:rw
    son_dataP = [son_dataP; 3 son_stackP(te,:)];
end
%% For oneVsALL data
train_data = [mom_data; dad_data; son_data];
test_data = [mom_dataP; dad_dataP; son_dataP];

[rw col] = size(train_data);
tr_rep = randperm(rw);
shuff_tr_d = train_data(tr_rep,:);
[rw col] = size(test_data);
te_rep = randperm(rw);
shuff_te_d = test_data(te_rep,:);

disp('PCA transformation complete. Packing data for the training process of machine learning algorithm ++++++++++++++++++++++');
disp('Processing ./.\./.\./.\./.\./.\./.\./.\');

%% FINALLT ML TRAINING WITH SVM AND PREDICTION ON TEST DATA

train_data = csvread('dataset_for_ML/train_data.csv');%dopp_momdadson_2ndHalf_train_data_imbalanced_withPCA_BlockConst_feat21to24_91to94.csv');%csvread('mom_dad_son_train_data_all.csv');%csvread('mom_dad_son_1stHalf_train_data_balanced_withPCA_feat21to29_took8frmPCA.csv');%csvread('mom_dad_son_2ndHalf_train_data_balanced_withPCA_feat21to29_took8frmPCA.csv');%csvread('mom_dad_son_2ndHalf_train_data_balanced_withPCA_feat1to10_0meanNormSD.csv');%csvread('mom_dad_son_2ndhalf_train_data_balanced_withPCA.csv');%%csvread('mom_dad_son_train_data_all_feat1_10PCA_took8.csv');% csvread('mom_dad_son_2ndhalf_train_data_balanced_withPCA.csv');%csvread('dad_v_son_train_data_1.csv');
test_data = csvread('dataset_for_ML/test_data.csv');%dopp_momdadson_2ndHalf_test_data_imbalanced_withPCA_BlockConst_feat21to24_91to94.csv');%csvread('mom_dad_son_test_data_all.csv');%csvread('mom_dad_son_1stHalf_test_data_balanced_withPCA_feat21to29__took8frmPCA.csv');%csvread('mom_dad_son_2ndHalf_test_data_balanced_withPCA_feat21to29__took8frmPCA.csv');%csvread('mom_dad_son_2ndHalf_test_data_balanced_withPCA_feat1to10_0meanNormSD.csv');%csvread('mom_dad_son_2ndhalf_test_data_balanced_withPCA.csv');%csvread('mom_dad_son_test_data_all.csv');%csvread('mom_dad_son_test_data_all_feat1_10PCA_took8.csv');%csvread('mom_dad_son_2ndhalf_test_data_balanced_withPCA.csv');%csvread('dad_v_son_test_data_1.csv');

%% preparing dataset
X_train = train_data(:,2:end); %
y_train = train_data(:,1);
X_test = test_data(:,2:end);
y_test = test_data(:,1);
% % %% CV partition
% % c = cvpartition(y_train,'k',10);

%% feature selection
disp('Running ML algorithm to classify and detect authorized users for the system . . . . . . . . . . . . . . . ');
% % % opts = statset('display','iter');
% % % % SVM Templates -----------------------------------------------------
% % % t = templateSVM('Standardize',false,'KernelFunction','rbf','IterationLimit', 1e3,'SaveSupportVectors',true);
% % % % % % Fitcecoc training
% % % classf = @(train_data, train_labels, test_data, test_labels)...
% % %     sum(predict(fitcecoc(train_data, train_labels,'Learners',t), test_data) ~= test_labels);
[rwX colX] = size(X_train);
fs = colX;
X_train_w_best_feature = X_train(:,fs);
% % % % % % % % %% Model GENRATOR---------------------------------------------------
% % % % % % % % Md1  = fitcecoc(X_train,y_train,'Learners', t ,'Coding','onevsone','OptimizeHyperparameters','all',...
% % % % % % % %     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
% % % % % % % %     'expected-improvement-plus','MaxObjectiveEvaluations',45, 'ShowPlots',false,  'Verbose', 0)); % 'Verbose',0 [Verbose 2 to show display of SVM processing
% % % % % % % % 
% % % % % % % % 
% % % % % % % % %% Final test with test set - if model is trained and instantly acceses for prediction
% % % % % % % % X_test_w_best_feature = X_test(:,fs);
% % % % % % % % predD = predict(Md1,X_test);
% % % % % % % % test_accuracy_for_iter = sum((predict(Md1,X_test) == y_test))/length(y_test)*100;
% % % % % % % % % Save trained Model ---------------------------------------------
% % % % % % % % fileName =strcat('svm_model/model_accuracy_',num2str(test_accuracy_for_iter),'_perc.mat');
% % % % % % % % save(fileName, 'Md1');
%% Load Model For Prediction =================================================
fileNAME = strcat('svm_model/model_accuracy_100_perc.mat');
Md1 = load(fileNAME, 'Md1');

%% If model is loaded then Md1 struct needs to be accessed like this ----
predD = predict(Md1.Md1,X_test);
test_accuracy_for_iter = sum((predict(Md1.Md1,X_test) == y_test))/length(y_test)*100;


% Confusion MAtrix Plot
figure(8);
cm = confusionchart(y_test,predD,'RowSummary','total-normalized');%
cm.NormalizedValues;
cm.XLabel = 'ML Inferred User No.';
cm.YLabel = 'Actual User No.';
cm.Title = strcat('12. ML algorithm recognized users with a accuracy of:  ' , num2str(test_accuracy_for_iter), ' %');
disp(['++++++++++++++++++++++++++  Accuracy :', num2str(test_accuracy_for_iter), ' % +++++++++++++++++++++++++++++++++++++']);

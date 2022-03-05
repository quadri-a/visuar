function [T1_norm, T2_norm, T3_norm, T4_norm, H1_norm, H2_norm, H3_norm, H4_norm, S1_norm, S2_norm, S3_norm, S4_norm, ...
    T1, T2, T3 ,T4, H1, H2, H3, H4, S1, S2, S3 ,S4] = fid_markers(Minval_1, Peakval_1, Minval_2, Peakval_2, Minval_3, Minloc_1, Peakloc_1, Minloc_2, Peakloc_2, Minloc_3, framePeriod,minWin,numChirps )

    % FIDUCIAL MARKERS _________________________
            %H1 is ST to AFP(systolic first)/ASP to VFP(diastolic first) height 
            H1 = Peakval_1 - Minval_1; % refers to the 1st peak val and min 1 val
%             H1_mm = H1*(3.9e-3/4*pi);
            %H2 is AFP to ASP(systolic first)/VFP to ED(diastolic first)
            %heigth
            H2 = Peakval_1 - Minval_2;%refers to the 1st peak val and min 2 val
%             H2_mm = H2*(3.9e-3/4*pi);
            %H3 is ASP to VFP(systolic first)/ED to AFP (diastolic
            %first)height
            H3 = Peakval_2 - Minval_2;% refers to 2nd peak val and min 2 val
%             H3_mm = H3*(3.9e-3/4*pi);
            %H4 is VFP to ED(systolic first)/AFP to ASP (diastolic first)
            H4 = Peakval_2 - Minval_3; % refers to 2nd peak val and min 3 val 
%             H4_mm = H4*(3.9e-3/4*pi);
            % Height differences
            Heights = [H1 H2 H3 H4];
            [mxval mxidx] = max(Heights);
            [mnval mnidx] = min(Heights);
            maxHeight = mxval;
            minHeight = mnval;
            disp(['Max height is H', num2str(mxidx)]);
            
            % Normalize Heights
            H1_norm = (H1)/(maxHeight);%(H1 - minHeight)/(maxHeight - minHeight );
            H2_norm = (H2)/(maxHeight);
            H3_norm = (H3)/(maxHeight);
            H4_norm = (H4)/(maxHeight );
            
            Pi = Minloc_3 - Minloc_1;
            cycle_duration = (Pi*(framePeriod))/numChirps;%minWin
            
            % Time differences
            % T1 is difference between 1st Peak and 1st min
            D1 = Peakloc_1 - Minloc_1;
            T1 = (D1*(framePeriod))/numChirps;
            T1_norm = T1/cycle_duration;
            
            % T2 is difference between 2nd Min and 1st Peak
            D2 = Minloc_2 - Peakloc_1;
            T2 = (D2*(framePeriod))/numChirps;
            T2_norm = T2/cycle_duration;
            
            % T3 is difference between 2nd Peak and 2nd min
            D3 = Peakloc_2 - Minloc_2;
            T3 = (D3*(framePeriod))/numChirps;
            T3_norm = T3/cycle_duration;
            
            % T4 is difference between 3rd/last Min and 2nd Peak
            D4 = Minloc_3 - Peakloc_2;
            T4 = (D4*(framePeriod))/numChirps;
            T4_norm = T4/cycle_duration;
            
            % Cardiac Output
            S1 = H1/T1;
            S2 = H2/T2;
            S3 = H3/T3;
            S4 = H4/T4;
            
            Slopes = [S1 S2 S3 S4];
            [maxSv maxSid] = max(Slopes);
%             S1_norm = S1/maxSv;
%             S2_norm = S2/maxSv;
%             S3_norm = S3/maxSv;
%             S4_norm = S4/maxSv;

            S1_norm = H1_norm/T1_norm;
            S2_norm = H2_norm/T2_norm;
            S3_norm = H3_norm/T3_norm;
            S4_norm = H4_norm/T4_norm;
            



end % end of funtion
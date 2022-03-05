% The input must contain phases of the ranges where each column is for a
% specific range . In other words , treating each column as an independent
% channel .
%
% UPDATEs :
% September 25, 2018: after adding or subtracting 2pi , all the phases to
% the end of the phase vector will be shifted by that amount . Because
%each
% shift corresponds to the unit circle shift up or down to a new circle .
% Thus , all later phases should be assumed to be in a new circle .
% Otherwise , the previous procedure did not allow to have a phase beyond
%of
% [-3pi ,3 pi ].
function UnwrappedPh = unwraphase(PhaseRange)
% unwrapping the phase
for i = 1: length ( PhaseRange ) -1
cols = PhaseRange(i+1 ,:) -PhaseRange(i ,:) > pi;
if (sum ( cols ) >0)
    colsIdx = find( cols ~=0) ;
for jj = 1: length( colsIdx )
PhaseRange(i+1: end , colsIdx (jj)) = PhaseRange(i+1: end,colsIdx (jj)) -2* pi;
end
end
cols1 = PhaseRange(i+1 ,:) -PhaseRange(i ,:) < -pi;
if (sum( cols1 ) >0)
colsIdx = find( cols1 ~=0) ;
for jj = 1: length( colsIdx )
PhaseRange(i+1: end , colsIdx (jj)) = PhaseRange(i+1: end , colsIdx (jj))+2* pi;
end
end
if cols1 == cols & sum ( cols1 ) > 0
display ('At the same time PhaseRange (i+1) -PhaseRange (i) <> \pi ')
end
end
UnwrappedPh = PhaseRange ;
end
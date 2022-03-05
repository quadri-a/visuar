function [output] = filter_IIR_biquad_heart(DataIn, pFilterCoefs, pScaleVals, pDelay, numStages)

numCoefStage = 6;
input = DataIn;

for indexStage = 1:numStages
    
    indexTemp = numCoefStage*(indexStage-1);
    
    b0 = pFilterCoefs((indexTemp+1) );
    b1 = pFilterCoefs((indexTemp+1) + 1);
    b2 = pFilterCoefs((indexTemp+1) + 2);
    a1 = pFilterCoefs((indexTemp+1) + 4);
    a2 = pFilterCoefs((indexTemp+1) + 5);

    scaleVal = pScaleVals(indexStage);
    
    pDelay(indexTemp+1) = scaleVal*input - a1*pDelay((indexTemp+1)+1 ) - a2*pDelay((indexTemp +1) +2 );
    y = b0*pDelay(indexTemp+1) + b1*pDelay((indexTemp+1)+1 ) + b2*pDelay((indexTemp+1) +2);
        
    pDelay((indexTemp+1)+2) = pDelay((indexTemp+1)+1);
    pDelay((indexTemp+1)+1) = pDelay((indexTemp+1));
   
    
    
    input = y;
end

output = y;

end
function [vSignal,vRand,vState,mP] = FiniteAutomate(iLength,varargin)
    vRand   = rand(iLength,1);
    vSignal = zeros(iLength,1);
    vState  = vSignal;
    vOutput = [0,1,2,3];
    mP = [0.4, 0.3, 0.0, 0.3; 
          0.3, 0.4, 0.3, 0.0; 
          0.0, 0.2, 0.7, 0.1;
          0.0, 0.0, 0.1, 0.9];
    iStateOffset = 0;
    if(~isempty(varargin))
        mP = [0.5, 0.5, 0.0; 
              0.0, 0.5, 0.5; 
              0.5, 0.25, 0.25];    
        %mP = mP';
        %vOutput = [4,5,6,7];
        vOutput = [0,5,2];
        iStateOffset = 4;
    end
    % ensure valid probabilities
    mP = bsxfun(@times,mP,1./sum(mP,2));
    iState = 1;
    vState(1) = iState;
    for(k=1:iLength)
        vSignal(k) = vOutput(iState);
        dRand = vRand(k);
        vCS = cumsum(mP(iState,:)); vTemp = ((vCS-dRand) > 0) .* (vCS-dRand); vTemp(vTemp==0) = Inf; [Y,iState] = min(vTemp);
        vState(k) = iState+iStateOffset;
    end
end
% test like
% [vSignal,vRand,vState] = FiniteAutomate(2000);
% vCounter = zeros(4,4); for(k=2:1999), vCounter(vState(k),vState(k+1)) = vCounter(vState(k),vState(k+1))+1;end
% repmat(1.0./sum(vCounter,2),1,4) .* vCounter
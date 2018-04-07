function [vSignal,vRand,vState,mP] = FiniteAutomate(iLength,varargin)
    vRand   = rand(iLength,1);
    vSignal = zeros(iLength,1);
    vState  = vSignal;
    vOutput = [4,5];
    mP = [0.5, 0.5; 
          0.2, 0.8];
    iStateOffset = 0;
    if(~isempty(varargin))
        mP = [0.3, 0.7; 
              0.8, 0.2];    
        vOutput = [4,5];
        iStateOffset = 2;
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
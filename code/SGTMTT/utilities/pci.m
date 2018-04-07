function [Y,beta] = pci(X,U,FI,oMetric)

% initialize using PRINCIPAL COMPONENTS.
%
%
% Options:	[Y, beta] = pci(X,U)
%           [W, beta] = pci(X,U,FI)
%
% assumes data already normalized

if nargin <3 || nargin >4
    error('The number of input arguments must be 3 or 4.');
end
    
[K,L] = size(U);

% Calculate the L first principal components and scale them by their respective eigenvalues
if(~isempty(oMetric.Par)) % thats for the case you call pci with metric scaling
    X=bsxfun(@times,X,oMetric.Par');
end
[eVts, eVls] = pca(X);
A = eVts(:,1:L)*diag(sqrt(eVls(1:L)));

% Normalise X to ensure 1:1 mapping of variances and calculate W
% as the solution of the equation: FI*W = normX*A'

normU = (U - ones(size(U))*diag(mean(U)))*diag(1./std(U)); % (X - mean)/std

if(sum(isnan(normU(:))))
    error('Latent grid contains nans after normalization - cant proceed')
end

if nargin==2
    Y = normU*A';
    interDistBeta = betai(Y);
else
    Y = FI \ (normU*A');
    Mplus1 = size(FI,2);
    Y(Mplus1,:) = mean(X);
    interDistBeta = betai(FI*Y); % keep regular square distance in betai - Y is already scaled above if needed
    %interDistBeta = betai(bsxfun(@times,(FI * Y),oMetric.Par'));
end

if (L < length(X(1,:))) 
    beta = min(interDistBeta,(1/eVls(L+1)));
else
    beta = interDistBeta;
end


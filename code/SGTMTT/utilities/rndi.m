function [Y,beta] = rndi(K, X, FI,oMetric)

% initialize using PRINCIPAL COMPONENTS.
%
%
% Options:	[Y, beta] = pci(X,U)
%           [W, beta] = pci(X,U,FI)
%

if nargin <3 || nargin >4
    error('The number of input arguments must be 3 or 4.');
end
[N D] = size(X);
varX=mean(std(X).^2);

Y = mvnrnd(zeros(D,1),varX*eye(D),K);
dist = oMetric.distance(Y,X);
%dist = sdist(Y,X);
beta = N*D/sum(sum((1/K)*dist));

if nargin==4
    Y = FI'*Y;
    Mplus1 = size(FI,2);
    Y(Mplus1,:) = mean(X);
end


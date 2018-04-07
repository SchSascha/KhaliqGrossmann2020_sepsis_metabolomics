function Pkl = trprobi(X, zeta,varargin)
% Evaluates prior transition probabilities (Kaban and Girolami's paper)
%
% Format:   Pkl = gtm_trprob(X, zeta)
%
% where:    Pkl is the prior transition probabilities matrix
%           X : hidden state vector
%           zeta: standard deviation
%
%Author:  Ivan Olier
%2005

if nargin==1
    zeta=0.3;
elseif nargin>3
    error('too many input arguments!');
end
if(isempty(varargin))
    dist = sdist(X,X);
else
    oMetric = varargin{1};
    dist = oMetric.distance(X, X);
end

Pkl = exp(-0.5*dist/zeta^2);
sumexp = sum(Pkl,2);
Pkl = Pkl./repmat(sumexp,[1 length(sumexp)]);


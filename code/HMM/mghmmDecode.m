function [pStates, loglik, fs, bs, s] = mghmmDecode(seq,tr,mu,sigma,varargin)
%MGHMMDECODE performs forward and backward algorithm (posterior decoding)
%on a multivariate Gaussian-emission HMM.

[covTypeStr, arg] = varArgRemove('cov', 1, varargin);
COV_DIAG = 0; COV_FULL = 1; covType = COV_DIAG; 
if strcmpi('full',covTypeStr), covType = COV_FULL; end

[numStates, checkTr] = size(tr);
if checkTr ~= numStates, error('TR matrix must be square.'); end
if numStates==0, error('Number of states must be positive.'); end
[nDim, statesMu] = size(mu);
if (statesMu ~= numStates-1)
    error('MU must have the same number of columns as TRANSITION.');
end
if covType == COV_FULL
    if (ndims(sigma) ~= 3), error('Full SIGMA must be 3D'); end
    [dimSigma, dimSigma2, statesSigma] = size(sigma);
    if (dimSigma ~= dimSigma2), error('Full SIGMA must be square'); end
else
    [dimSigma, statesSigma] = size(sigma);
end
if (statesSigma ~= numStates-1)
    error('SIGMA matrix must have the same number of columns as TRANSITION.');
elseif (nDim ~= dimSigma)
    error('SIGMA matrix must have the same number of rows (dimensions).');
end

emilik = mghmmEmiLik(seq, mu, sigma, varargin{:});
if all(max(emilik) > 0)
    [pStates, loglik, fs, bs, s] = hmmDecodeLik(tr, emilik);
else
    L = size(seq, 2); pStates = zeros(numStates, L); loglik = -1e100; 
    fs = zeros(numStates, L+1); bs = zeros(numStates, L+1); s = zeros(1,L);
end

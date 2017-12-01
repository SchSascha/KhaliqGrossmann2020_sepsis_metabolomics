function [estStates, loglikSeq] = mghmmViterbi(seq,tr,mu,sigma,varargin)
%MGHMMVITERBI performs Viterbi algorithm on a multivariate
%Gaussian-emission HMM. 

[numStates, checkTr] = size(tr);
if checkTr ~= numStates, error('TR matrix must be square.'); end
[nDim, statesMu] = size(mu);
[dimSigma, statesSigma] = size(sigma);
if (statesMu ~= numStates-1 || statesSigma ~= numStates-1)
    error('MU/SIGMA matrix must have the same number of columns as TR.');
elseif (nDim ~= dimSigma)
    error('MU/SIGMA matrix must have the same number of rows (dimensions).');
end
if numStates==0, error('Number of states must be positive.'); end

[estStates, loglikSeq] = hmmViterbiLik(tr, mghmmEmiLik(seq, mu, sigma, varargin{:}));

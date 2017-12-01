function [estStates, logP] = hmmViterbiLik(tr, emilik)
%HMMVITERBILIK performs the Viterbi algorithm on a HMM according to the
%specified emission likelihood. 
%   [ESTSTATES, LOGP] = HMMVITERBILIK(TR,EMILIK)
%   calculates the most likely path of hidden states of a given sequence. 
%   The sequence and the emission probability is
%   incorporated in the emission likelihood matrix, EMILIK. The sequence
%   can be discrete or continuous, as long as the emission probability is
%   set up accordingly (e.g. multinomial or normal).
%   TR(i,j) is the transition probability from state i to j, and
%   EMILIK(i,t) is the likelihood of observed value at t given its state i. 
%   LOGP is the joint probability of observing the optimal path and the
%   observation.
%   A hidden state and corresponding emission (observation) is added before
%   the first observation, at time zero (t=0). The state at time zero is
%   assumed to be the first state, so the emission probability from the
%   first to all other states is equivalent to the initial probability. 

%   This script is modified from hmmDecodeLik.m and MathWorks hmmviterbi.m.

[numStates, checkTr] = size(tr);
if checkTr ~= numStates, error('TR matrix must be square.'); end
if size(emilik,1) ~= numStates, 
    error('EMILIK matrix size must agree with TR matrix.'); 
end

L = size(emilik, 2);
w = warning('off'); % get log of zero warnings
logTR = log(tr);% work in log space to avoid numerical issues
logE = log(emilik');
warning(w);

pTR = zeros(numStates, L);
% assumption is that model is in state 1 at step 0
v = repmat(-inf, 1, numStates);
v(1,1) = 0;
for count = 1 : L
    [maxValue, maxPtr] = max(repmat(v', 1, numStates) + logTR);
    v = maxValue + logE(count,:);
    ptr(:, count) = maxPtr';
    if any(isnan(v)), 
        warning('hmmViterbiLik:vNaN','pr-max-path include NaN'); 
    end
end
[logP, maxPtr] = max(v);
estStates(L) = maxPtr;
for count = L - 1 : -1 : 1
    estStates(count) = ptr(estStates(count + 1), count + 1);
end

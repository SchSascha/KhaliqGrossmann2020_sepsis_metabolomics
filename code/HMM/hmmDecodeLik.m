function [pStates, loglikSeq, fs, bs, s] = hmmDecodeLik(tr, emilik)
%HMMDECODELIK performs forward and backward algorithm (posterior decoding)
%on a HMM according to emission likelihood of each state.
%   [PSTATES, LOGLIKSEQ, FORWARD, BACKWARD, S] = HMMDECODELIK(TR,EMILIK)
%   calculates the posterior probability of a given sequence (time series),
%   the log likelihood of the sequence, and the forward and backward 
%   probabilities of the sequence scaled by S. In Rabiner's notations, 
%   PSTATES is gamma matrix, and FORWARD and BACKWARD are alpha and beta
%   matrix, respectively. The sequence and the emission probability is
%   incorporated in the emission likelihood matrix, EMILIK. The sequence
%   can be discrete or continuous, as long as the emission probability is
%   set up accordingly (e.g. multinomial or normal).
%   TR(i,j) is the transition probability from state i to j, and
%   EMILIK(i,t) is the likelihood of observed value at t given its state i. 
%   LOGLIKSEQ is the marginal probability of the observation given the model.
%   A hidden state and corresponding emission (observation) is added before
%   the first observation, at time zero (t=0). The state at time zero is
%   assumed to be the first state, so the emission probability from the
%   first to all other states is equivalent to the initial probability. 
%
%   The actual forward probabilities can be recovered by using:
%        f = FORWARD.*repmat(cumprod(s),size(FORWARD,1),1);
%   The actual backward probabilities can be recovered by using:
%       bscale = fliplr(cumprod(fliplr(S)));
%       b = BACKWARD.*repmat([bscale(2:end), 1],size(BACKWARD,1),1);

%   This script is modified from MathWorks hmmdecode.m.

[numStates, checkTr] = size(tr);
if checkTr ~= numStates, error('TR matrix must be square.'); end
if size(emilik,1) ~= numStates, 
    error('EMILIK matrix size must agree with TR matrix.'); 
end

% shift the sequence to make algorithm cleaner at f0 and b0
emilik = [NaN(numStates,1), emilik];
L = size(emilik, 2);

fs = zeros(numStates,L);
fs(1,1) = 1;  % assume that we start in state 1.
s = zeros(1,L); % log is numerically instable, so use a scaling factor
s(1) = 1;
for count = 2:L
    fs(:,count) = (tr' * fs(:,count-1)) .* emilik(:,count);
    % scale factor normalizes sum(fs,count) to be 1. 
    s(count) =  sum(fs(:,count));
    if s(count) < 0
        % if seq likelihood is 0, return (and have caller skip this seq)
        % following HMMER implementation (not documented)
        %turn off warning so that dbstop at caller
        %warning('hmmDecodeLik:Likeli0', 'Sequence likelihood close to zero');
        pStates = []; loglikSeq = -Inf; fs = []; bs = []; return; 
    end
    fs(:,count) =  fs(:,count)./s(count);
%     if ~all(isfinite(fs(:,count)))
%         warning('hmmDecodeLik:FsNaN', 'FS include NaN'); 
%     end
end

%  The  actual forward and  probabilities can be recovered by using
%   f = fs.*repmat(cumprod(s),size(fs,1),1);

bs = ones(numStates,L);
for count = L-1:-1:1
    bs(:,count) = tr * (bs(:,count+1) .* emilik(:,count+1)) / s(count+1);
%     if ~all(isfinite(bs(:,count))) 
%         warning('hmmDecodeLik:BsNaN', 'BS include NaN'); 
%     end
end

%  The  actual backward and  probabilities can be recovered by using
%  scales = fliplr(cumprod(fliplr(s)));
%  b = bs.*repmat([scales(2:end), 1],size(bs,1),1);

loglikSeq = sum(log(s));
pStates = fs.*bs;

% get rid of the column that we stuck in to deal with the f0 and b0 
pStates(:,1) = [];

if size(fs,1)==0 || size(fs,2)==0, error('fs empty'); end
if size(bs,1)==0 || size(bs,2)==0, error('bs empty'); end

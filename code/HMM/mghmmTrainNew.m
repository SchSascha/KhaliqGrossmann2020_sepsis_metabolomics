function [tr mu sigma loglikArray mghmmArr] = mghmmTrainNew(seqs, varargin)
%MGHMMTRAINNEW trains a multivariate Gaussian-emission HMM from the sequences.
%   [TR1 MU1 SIGMA1 LOGLIKS] = MGHMMTRAINNEW(SEQS)
%   Multiple restarts with random initialization is employed to overcome
%   local optimum. SEQS{n}(g,t) or SEQS(n,g,t) denoting the observed value
%   at time t of sample n.

% Process arguments
[nState, algArg] = varArgRemove('states', 3, varargin);
nState = nState + 1;
[nRepeat, algArg] = varArgRemove('replicates', 1, algArg);
[modelStr, algArg] = varArgRemove('model', 'loophmm', algArg);
isLoopHmm = false; isJumpHmm = false; isEqLenHmm = false; isLoopJumpHmm = false;
if     strcmpi('jumphmm', modelStr), isJumpHmm  = true; 
elseif strcmpi('eqlenhmm',modelStr), isEqLenHmm = true;
elseif strcmpi('loopjumphmm',modelStr), isLoopJumpHmm = true;
else                                 isLoopHmm  = true; 
end

covTypeStr = varArgRemove('cov', 'diag', algArg);
COV_DIAG = 0; COV_FULL = 1; covType = COV_DIAG; 
if strcmpi('full',covTypeStr), covType = COV_FULL; end

pseudoStay = .1;

if isnumeric(seqs) && ndims(seqs)==3
    [numSeqs, nDim, maxL] = size(seqs);
elseif iscell(seqs)
    numSeqs = numel(seqs);
    nDim = size(seqs{1}, 1);
    maxL = 0;
    for count = 1:numSeqs
        seq = seqs{count};
        if size(seq,2) > maxL, maxL = size(seq,2); end
    end
else
    error('SEQS must be a numeric 3D matrix or a cell array.');
end
if isLoopHmm
    [pars, cumPars, nRepeat] = loopHmmPrepare(nState, maxL, nRepeat);
elseif isLoopJumpHmm
    [pars, cumPars, nRepeat] = loopHmmPrepare(nState, maxL, nRepeat);
elseif isJumpHmm
    [tr0, sigma0, muArr, nRepeat] = jumpHmmPrepare(seqs, numSeqs, nDim, ...
                                     nState, maxL, covType, varargin{:});
%elseif isEqLenHmm, do nothing
end

if nargout > 4, mghmmArr = cell(nRepeat, 1); end
loglikArray = cell(nRepeat, 1);
topLogLiks = -Inf;
if nRepeat == 0, error('replicates cannot be zero'); end
for iRepeat = 1 : nRepeat
    if isLoopHmm
        [tr1, mu1, sigma1] = loopHmmInit(seqs, numSeqs, nDim, nState, maxL, ...
                pseudoStay, pars, cumPars, iRepeat, covType, varargin{:});
        %if nDim > 1, disp('sigma1:'); disp(sigma1); end
    elseif isJumpHmm
        tr1 = tr0; sigma1 = sigma0; mu1 = muArr{iRepeat};
    elseif isLoopJumpHmm
        [tr1, mu1, sigma1] = loopJumpHmmInit(seqs, numSeqs, nDim, nState, maxL, ...
                pseudoStay, pars, cumPars, iRepeat, covType, varargin{:});
    elseif isEqLenHmm
        [tr1, mu1, sigma1] = eqLenHmmInit(seqs, numSeqs, nDim, maxL, ...
                                           iRepeat, covType, varargin{:});
    end
    %sigma1(sigma1 <= 0) = 1;
    [tr1,mu1,sigma1,loglik1] = mghmmTrain(seqs,tr1,mu1,sigma1,varargin{:});
    if nargout > 4
        mghmmArr{iRepeat}.tr = tr1; 
        mghmmArr{iRepeat}.mu = mu1; 
        mghmmArr{iRepeat}.sigma = sigma1;
    end
    loglikArray{iRepeat} = loglik1;
    %if ~isfinite(max(loglik1)), error('nan', 'loglik1 nan'); end
    if max(loglik1) > max(topLogLiks)
        topLogLiks = loglik1; tr = tr1; mu = mu1; sigma = sigma1;
    end
end
if ~(max(topLogLiks,2) > -Inf), error('LOGLIKS never assigned'); end

function [pars, cumPars, nRepeat] = loopHmmPrepare(nState, maxL, nRepeat)
    if nRepeat == 1
        pars = listAvgPartition(maxL, nState-1);
    else
        pars = listPartitions(maxL, nState-1);
        parsIdx = randperm(size(pars,1));
        nRepeat = min(size(pars,1), nRepeat);
        pars = pars(parsIdx(1 : nRepeat), :);
    end
    cumPars = [zeros(size(pars,1),1), cumsum(pars, 2)];

function pars = listAvgPartition(L, N)
    pars = floor(L/N:L/N:L) - floor(0:L/N:L-L/N); % N x 1

function [tr1, mu1, sigma1] = loopHmmInit(seqs, numSeqs, nDim, nState, maxL, ...
                              pseudoStay, pars, cumPars, iRepeat, covType, varargin)
    [glbCovRt, algArg] = varArgRemove('rtglbcov', 0.5, varargin);
    glbCov = varArgRemove('glbcov', [], algArg);
    escape = ((nState-1)*pseudoStay+maxL)./maxL./(pseudoStay+pars(iRepeat,:))';    
    tr1 = diag([0; 1 - escape; 1]) + diag([1; escape], 1);
    %no terminal state: tr1 = diag([0; 1 - escape(1:nState-2); 1]) + diag([1; escape(1:nState-2)], 1);
    if covType == 1 %COV_FULL
        totalPost = zeros(nDim, nState-1);
        meanSeq = zeros(nDim, nState-1); 
        nCov = zeros(nState-1, 1);
        covSeq = zeros(nDim, nDim, nState-1); 
        for count = 1 : numSeqs
            seq = seqs{count};
            for j = 1 : nState-1
                seq1 = seq(:, cumPars(iRepeat,j)+1 : ...
                                 min(size(seq, 2), cumPars(iRepeat,j+1)));
                isFin = isfinite(seq1); seq1(~isFin) = 0;
                meanSeq(:,j) = meanSeq(:,j) + sum(seq1,2);
                totalPost(:,j) = totalPost(:,j) + sum(isFin, 2);
            end
        end
        mu1 = [zeros(nDim,1), meanSeq ./ totalPost];
        for count = 1 : numSeqs
            seq = seqs{count}; 
            for j = 1 : nState-1
                seq1 = seq(:, cumPars(iRepeat,j)+1 : ...
                                 min(size(seq, 2), cumPars(iRepeat,j+1)));
                ut = ~isfinite(seq1); 
                for t=1:size(seq1,2), seq1(ut(:,t),t) = mu1(ut(:,t), j+1); end
                covSeq(:,:,j) = covSeq(:,:,j) + seq1 * seq1';
                nCov(j) = nCov(j) + size(seq1,2);%FIXME
            end
        end
        sigma1 = zeros(nDim, nDim, nState);
        sigma1(:, :, 1) = eye(nDim);
        for j = 1 : nState-1
            sigma1(:,:,j+1) = glbCovRt * glbCov + (1-glbCovRt) ...
                * (covSeq(:,:,j)/nCov(j)-mu1(:,j+1)*mu1(:,j+1)');
            [upper,part] = chol(sigma1(:,:,j+1)); %check positive definite
            if part~=0, sigma1(:,:,j+1) = sigma1(:,:,j+1) + 0.01*eye(nDim); end
        end
    else
        sumPars = zeros(nDim, nState-1);
        sum2Pars = zeros(nDim, nState-1);
        for count = 1 : numSeqs
            seq = seqs{count};
            for j = 1 : nState-1
                for iVar = 1 : nDim
                    seq1 = seq(iVar, cumPars(iRepeat,j)+1 : ...
                                     min(size(seq, 2), cumPars(iRepeat,j+1)));
                    seq1 = seq1(isfinite(seq1));
                    if ~isempty(seq1), 
                        sumPars(iVar,j) = sumPars(iVar,j)+mean(seq1);
                        sum2Pars(iVar,j) = sum2Pars(iVar,j)+mean(seq1.^2);
                    end
                end
            end
        end
        mu1 = [zeros(nDim,1), sumPars / numSeqs];
        sigma1 = [ones(nDim,1), sum2Pars / numSeqs - mu1(:,2:nState).^2]; 
    end

function [tr1, mu1, sigma1] = loopJumpHmmInit(seqs, numSeqs, nDim, nState, maxL, ...
                              pseudoStay, pars, cumPars, iRepeat, covType, varargin)
    [glbCovRt, algArg] = varArgRemove('rtglbcov', 0.5, varargin);
    glbCov = varArgRemove('glbcov', [], algArg);
    %esc_seq = (1./(nState:-1:1))';
    %tr1 = diag([0; esc_seq]);
    %for k = 1:nState
    %    tr1 = tr1 + diag([esc_seq(1); esc_seq(1:(nState-k))], k);
    %end
    
    escape = ((nState-1)*pseudoStay+maxL)./maxL./(pseudoStay+pars(iRepeat,:))';    
    tr1 = diag([0; 1 - escape; 1]) + diag([1; escape], 1);
    tr1(1, 2:nState) = 1/(nState-1);
    
    %no terminal state: tr1 = diag([0; 1 - escape(1:nState-2); 1]) + diag([1; escape(1:nState-2)], 1);
    if covType == 1 %COV_FULL
        totalPost = zeros(nDim, nState-1);
        meanSeq = zeros(nDim, nState-1); 
        nCov = zeros(nState-1, 1);
        covSeq = zeros(nDim, nDim, nState-1); 
        for count = 1 : numSeqs
            seq = seqs{count};
            for j = 1 : nState-1
                seq1 = seq(:, cumPars(iRepeat,j)+1 : ...
                                 min(size(seq, 2), cumPars(iRepeat,j+1)));
                isFin = isfinite(seq1); seq1(~isFin) = 0;
                meanSeq(:,j) = meanSeq(:,j) + sum(seq1,2);
                totalPost(:,j) = totalPost(:,j) + sum(isFin, 2);
            end
        end
        mu1 = [zeros(nDim,1), meanSeq ./ totalPost];
        for count = 1 : numSeqs
            seq = seqs{count}; 
            for j = 1 : nState-1
                seq1 = seq(:, cumPars(iRepeat,j)+1 : ...
                                 min(size(seq, 2), cumPars(iRepeat,j+1)));
                ut = ~isfinite(seq1); 
                for t=1:size(seq1,2), seq1(ut(:,t),t) = mu1(ut(:,t), j+1); end
                covSeq(:,:,j) = covSeq(:,:,j) + seq1 * seq1';
                nCov(j) = nCov(j) + size(seq1,2);%FIXME
            end
        end
        sigma1 = zeros(nDim, nDim, nState);
        sigma1(:, :, 1) = eye(nDim);
        for j = 1 : nState-1
            sigma1(:,:,j+1) = glbCovRt * glbCov + (1-glbCovRt) ...
                * (covSeq(:,:,j)/nCov(j)-mu1(:,j+1)*mu1(:,j+1)');
            [upper,part] = chol(sigma1(:,:,j+1)); %check positive definite
            if part~=0, sigma1(:,:,j+1) = sigma1(:,:,j+1) + 0.01*eye(nDim); end
        end
    else
        sumPars = zeros(nDim, nState-1);
        sum2Pars = zeros(nDim, nState-1);
        for count = 1 : numSeqs
            seq = seqs{count};
            for j = 1 : nState-1
                for iVar = 1 : nDim
                    seq1 = seq(iVar, cumPars(iRepeat,j)+1 : ...
                                     min(size(seq, 2), cumPars(iRepeat,j+1)));
                    seq1 = seq1(isfinite(seq1));
                    if ~isempty(seq1), 
                        sumPars(iVar,j) = sumPars(iVar,j)+mean(seq1);
                        sum2Pars(iVar,j) = sum2Pars(iVar,j)+mean(seq1.^2);
                    end
                end
            end
        end
        mu1 = [zeros(nDim,1), sumPars / numSeqs];
        sigma1 = [ones(nDim,1), sum2Pars / numSeqs - mu1(:,2:nState).^2]; 
    end
    
function [tr0, sigma0, muArr, nRepeat] = jumpHmmPrepare(seqs, numSeqs, nDim, ...
                                         nState, maxL, covType, varargin)
    [isTermState, algArg] = varArgRemove('termstate', 0, varargin);
    [glbCovRt, algArg] = varArgRemove('rtglbcov', 0.5, algArg);
    glbCov = varArgRemove('glbcov', [], algArg);
    
    nState = nState + 1;
    tr0 = diag([.1*ones(1,nState-2) 1], 1) ...
          + diag([.7*ones(1,nState-3) .9], 2) + diag(.2*ones(1,nState-3), 3);
    tr0(nState, nState) = 1;
    
    total = zeros(nDim, maxL); y1 = zeros(nDim, maxL); y2 = zeros(nDim, maxL); 
    for q = 1 : numSeqs
        seq = seqs{q}; L = size(seq, 2);
        finiteSeq = seq; isFinSeq = isfinite(seq); finiteSeq(~isFinSeq) = 0;
        total(:,1:L) = total(:,1:L) + isFinSeq; 
        y1(:,1:L) = y1(:,1:L) + finiteSeq; 
        if covType ~= 1, y2(:,1:L) = y2(:,1:L) + finiteSeq.^2; end
    end
    mu = y1 ./ total;
    if any(any(~isfinite(mu))), error('jumpHmmPrepare, mu has NaN'); end
    if covType == 1 %COV_FULL, global covariance already calculated
        sigma0 = [eye(nDim), repmat(glbCov, [1 1 nState-2])];
    else % calculate global variance; set as Sigma
        sigma0 = sum(y2,2) ./ sum(total,2) - (sum(y1,2) ./ sum(total,2)).^2;
        sigma0 = [ones(nDim,1), repmat(sigma0, 1, nState-2)];
    end
    
    % fill missing values by avg of same time; pick a random q; set as mu's
    for q=1:numSeqs, isFin=isfinite(seqs{q}); seqs{q}(~isFin)=mu(~isFin); end
    emilik = [zeros(1, maxL); ones(nState-2, maxL); zeros(1,maxL)];
    if isTermState, emilik = [emilik, [zeros(nState-1,1); 1]]; end
    expectedGamma = hmmDecodeLik(tr0, emilik);
    seqs = seqs(randperm(numSeqs));
    muArr = cell(numSeqs, 1);
    for q = 1 : numSeqs
        finiteP = expectedGamma(2 : nState-1, 1 : size(seqs{q},2))';
        finiteP(end, all(finiteP==0)) = 1; %if cannot reach, use last
        mu1 = (seqs{q} * finiteP) ./ repmat(sum(finiteP), nDim, 1);
        muArr{q} = [zeros(nDim,1), mu1];
        if any(any(~isfinite(muArr{q}))) 
            fprintf('q is %2.0f', q);
            error('jumpHmmPrepare, muArr NaN');
        end
    end
    if any(any(~isfinite(tr0))), error('jumpHmmPrepare, tr0 has NaN'); end
    if any(any(~isfinite(sigma0))), error('jumpHmmPrepare, sigma0 has NaN'); end
    nRepeat = length(muArr);

function [tr1, mu1, sigma1] = eqLenHmmInit(seqs, numSeqs, nDim, maxL, ...
                                           iRepeat, covType, varargin)
%EQLENHMMINIT create equal-length one-state-per-time-point HMM
%(in fact not a HMM but a Gaussian Naive Bayes model)
    [glbCovRt, algArg] = varArgRemove('rtglbcov', 0.5, varargin);
    glbCov = varArgRemove('glbcov', [], algArg);
    if iRepeat>1, error('EqLenHMM does not support random initialization'); end
    nState = maxL + 2;
    tr1 = diag(ones(nState - 1, 1), 1); tr1(nState, nState) = 1;
    totalPost = zeros(nDim, nState-2); 
    Mu = zeros(nDim, nState-2);
    if covType == 1 %COV_FULL
        covY = zeros(nDim, nDim, nState-2);
        coTotalP = zeros(nDim, nDim, nState-2);
    else
        y2 = zeros(nDim, nState-2);
    end
    for count = 1 : numSeqs
        seq = seqs{count};
        finiteSeq = seq; 
        isFinSeq = isfinite(seq);
        finiteSeq(~isFinSeq) = 0;
        finiteP = eye(nState - 2, size(seq, 2));
        Mu = Mu + finiteSeq * finiteP';
        totalPost = totalPost + isFinSeq * finiteP';
        if covType == 1 %COV_FULL
            for j = 1 : nState-2
                covY(:,:,j) = covY(:,:,j) + finiteSeq(:,j) * finiteSeq(:,j)';
                coTotalP(:,:,j) = coTotalP(:,:,j)+isFinSeq(:,j)*isFinSeq(:,j)';
            end
        else
            y2 = y2 + (finiteSeq.^2) * finiteP';
        end
    end
    mu1 = [zeros(nDim,1), Mu ./ totalPost]; 
    if covType == 1 %COV_FULL
        sigma1(:,:,1) = eye(nDim);
        for j = 2 : nState-1
            sigma1(:,:,j) = glbCovRt * glbCov + (1 - glbCovRt) ...
                * (covY(:,:,j-1)./coTotalP(:,:,j-1) - mu1(:,j-1)*mu1(:,j-1)');
            sigma1(:,:,j) = (sigma1(:,:,j) + sigma1(:,:,j)') / 2;
        end
        if any(any(any(~isfinite(sigma1)))), error('eqLenHmmInit, init sigma include NaN'); end
        if any(any(any(sigma1 < 0))), error('eqLenHmmInit, init sigma include negative'); end
    else
        sigma1 = [ones(nDim,1), y2 ./ totalPost - mu1(:,2:nState-1).^2];
        if any(any(~isfinite(sigma1))), error('eqLenHmmInit, init sigma include NaN'); end
        if any(any(sigma1 < 0)), error('eqLenHmmInit, init sigma include negative'); end
    end

function [mghmms, logCondLiksArr, trainAcc] = mmieGenTrain(trainX, ...
         mghmms, varargin)

[tol, maxiter, verbose, isCorrTr, learnRtType, learnRtTrType, isJumpHmm, ...
    covType] = processArguments(varargin{:});
verifyInput(mghmms, covType);

logCondLiks = 1; % make sure not to converge after first step
logCondLiksArr = nan(maxiter,1);
trainAcc = nan(maxiter,1);
for iter = 1 : maxiter
    oldMghmms = mghmms; oldLogCondLiks = logCondLiks;
    [numGamma,numO1,numO2,numXi,genGamma,genO1,genO2,genXi,trainErr,logCondLiks]...
        = calExptCount(trainX,mghmms,isCorrTr,isJumpHmm,varargin{:});
    learnRt = calLearnRate(learnRtType, iter, trainErr);
    learnRtTr = calLearnRate(learnRtTrType, iter, trainErr);
    logCondLiksArr(iter) = logCondLiks; trainAcc(iter) = 1 - trainErr;
    if verbose, verboseOutput(iter, logCondLiks); end
    if learnRt < 1e-5, 
        if verbose, fprintf(['Algorithm converged after %d iterations'...
                ' (no training error).'], iter); end
        break;
    end
    if covType == 1
        mghmms = updateParameterCovFull(mghmms, numGamma, numO1, numXi, ...
                 genGamma, genO1, genXi, learnRt, learnRtTr);
    else
        mghmms = updateParameter(mghmms, numGamma, numO1, numO2, numXi, ...
                 genGamma, genO1, genO2, genXi, learnRt, learnRtTr);
    end
    if isConvergent(mghmms, oldMghmms, logCondLiks, oldLogCondLiks, tol)
        if verbose
            fprintf('Algorithm converged after %d iterations.',iter);
        end
        break; 
    end
end
logCondLiksArr(isnan(logCondLiksArr)) = [];
trainAcc(isnan(trainAcc)) = [];

%---------------------------------------------------------------------

function verifyInput(mghmms, covType)
    for c = 1 : length(mghmms)
        [numStates, checkTr] = size(mghmms{c}.tr);
        if checkTr ~= numStates, error('TRANSITION matrix must be square.'); end
        [nDim, statesMu] = size(mghmms{c}.mu);
        if covType == 1
            if (ndims(mghmms{c}.sigma)~=3), error('Full SIGMA must be 3D'); end
            [dimSigma, dimSigma2, statesSigma] = size(mghmms{c}.sigma);
            if (dimSigma ~= dimSigma2), error('Full SIGMA must be square'); end
        else
            [dimSigma, statesSigma] = size(mghmms{c}.sigma);
        end
        if (statesMu ~= numStates-1 || statesSigma ~= numStates-1)
            error('MU/SIGMA matrix must have the same number of columns as TRANSITION.');
        elseif (nDim ~= dimSigma)
            error('MU/SIGMA matrix must have the same number of rows (dimensions).');
        end
        if numStates==0, error('Number of states must be positive.'); end
    end

function [tol, maxiter, verbose, isCorrTr, learnRtType, learnRtTypeTr, ...
    isJumpHmm, covType] = processArguments(varargin)
    % default values
    tol.lcl = 1e-6; tol.mu = tol.lcl; tol.sigma = tol.lcl; 
    maxiter = 100; %500;
    verbose = false;
    isCorrTr = false;
    learnRtType = 'counted';
    learnRtTypeTr = 'counted';
    covType = 0;

    modelStr = varArgRemove('model', 'loophmm', varargin);
    isJumpHmm = strcmp(modelStr, 'jumphmm');

    if nargin > 0
        if rem(nargin,2) == 1, 
            error('Incorrect number of arguments to %s.',mfilename);
        end
        okargs = {'mmietol','mmieiterations','mmieverbose','mmietolmu',...
            'mmietolsigma','mmiecorr','mmielearnrate', 'mmietrlearnrate', 'cov'};
        for j = 1: 2 : nargin
            pname = varargin{j};
            pval = varargin{j+1};
            k = strmatch(lower(pname), okargs);
            if isempty(k), continue;%error('Unknown parameter name:%s.',pname);
            elseif length(k)>1, error('Ambiguous parameter name:  %s.',pname);
            end
            switch(k)
            case 1  % tolerance
                tolLcl = tol.lcl;
                tol.lcl = pval;
                if tol.mu == tolLcl, tol.mu = tol.lcl; end
                if tol.sigma == tolLcl, tol.sigma = tol.lcl; end
            case 2 % max iterations
                maxiter = pval;
            case 3 % verbose
                if islogical(pval) || isnumeric(pval), verbose = pval; end
            case 4 % mu tolerance
                tol.mu = pval;
            case 5 % sigma tolerance
                tol.sigma = pval;
            case 6 % corrective training
                if islogical(pval) || isnumeric(pval), isCorrTr = pval; end
            case 7 % learning rate for emission
                if ischar(pval) && ~strcmp('counted', pval) ...
                        && ~strcmp('inverse', pval)
                    error('Invalid value for LearnRate: %s', pval);
                end
                learnRtType = pval; %character or numeric
            case 8 % learning rate for transition
                if ischar(pval) && ~strcmp('counted', pval) ...
                        && ~strcmp('inverse', pval)
                    error('Invalid value for LearnRate: %s', pval);
                end
                learnRtTypeTr = pval; %character or numeric
            case 9 % cov
                if     strcmpi('diag',pval), covType = 0; 
                elseif strcmpi('full',pval), covType = 1; 
                end
            end
        end
    end
    
function [numGamma, numO1, numO2, numXi, genGamma, genO1, genO2, genXi, ...
    trainErr,logCondLiks] = calExptCount(seqs,mghmms,isCorrTr,isJumpHmm,varargin)
    nCate = length(mghmms);
    logCondLiks = 0;
    genHMM = buildGeneralModel(mghmms, isJumpHmm);
    genNState = size(genHMM.tr, 1); genNDim = size(genHMM.mu, 1);
    genGamma = zeros(genNDim, genNState-1); % general model
    genO1 = zeros(genNDim, genNState-1); genO2 = zeros(genNDim, genNState-1);
    genXi = zeros(genNState, genNState);
    numGamma = cell(size(mghmms)); % numerator models
    numO1 = cell(size(mghmms)); numO2 = cell(size(mghmms));
    numXi = cell(size(mghmms));
    trainErr = 0; nSeqs = 0; 
    logJoint = zeros(nCate, 1); pState = NaN;
    for c = 1 : nCate
        nSeqs = nSeqs + length(seqs{c});
        numNState = size(mghmms{c}.tr, 1); numNDim = size(mghmms{c}.mu, 1);
        numGamma{c} = zeros(numNDim, numNState-2);
        numO1{c}    = zeros(numNDim, numNState-2); 
        numO2{c}    = zeros(numNDim, numNState-2);
        numXi{c}    = zeros(numNState-1, numNState-1);
        for iSeq = 1 : length(seqs{c})
            seq = seqs{c}{iSeq}; 
            if nCate ~= 2
                for m = 1 : nCate
                    emilik1 = mghmmEmiLik(seq, mghmms{m}.mu, mghmms{m}.sigma, varargin{:});
                    [pState1,loglik1,fs1,bs1,scale1] = hmmDecodeLik(...
                        mghmms{m}.tr,emilik1);
                    if m == c, 
                        pState=pState1; fs=fs1; bs=bs1; scale=scale1; 
                        emilik = emilik1;
                    end
                    logJoint(m) = loglik1 + log(mghmms{m}.prior);
                end
                genLogLik1 = log(sum(exp(logJoint)));
            else
                emilik = mghmmEmiLik(seq, mghmms{c}.mu, mghmms{c}.sigma, varargin{:});
                genEmiLik = mghmmEmiLik(seq, genHMM.mu, genHMM.sigma, varargin{:});
                [pState,loglik1,fs,bs,scale]=hmmDecodeLik(mghmms{c}.tr,emilik);
                % if seq likelihood is 0, skip this seq, following HMMER
                if size(pState,1) == 0, continue; end %seq likelihood zero
                [genPState,genLogLik1,genFs,genBs,genScale] = hmmDecodeLik(...
                    genHMM.tr, genEmiLik);
                if size(genPState,1) == 0, continue; end %seq likelihood zero
                logJoint(c) = loglik1 + log(mghmms{c}.prior);
                diffJoint = exp(genLogLik1) - exp(logJoint(c));
                if diffJoint >= eps, logJoint(3 - c) = log(diffJoint);
                else logJoint(3 - c) = -Inf;
                end                
            end
            logCondLiks = logCondLiks + logJoint(c) - genLogLik1;
            if ~isCorrTr || any(logJoint > logJoint(c))
                if nCate ~= 2
                    genEmiLik = mghmmEmiLik(seq, genHMM.mu, genHMM.sigma, varargin{:});
                    [genPState, genLogLik1, genFs, genBs, genScale] ...
                        = hmmDecodeLik(genHMM.tr, genEmiLik);
                end
                if any(logJoint > logJoint(c)), trainErr = trainErr + 1; end
                [totGamma, eO1, eO2, totXi] = calSuffStat(seq, pState, fs, ...
                                              bs, scale, emilik, mghmms{c}.tr);
                numGamma{c} = numGamma{c} + totGamma(:, 2:numNState-1);
                numO1{c} = numO1{c} + eO1(:, 2:numNState-1); 
                numO2{c} = numO2{c} + eO2(:, 2:numNState-1);
                numXi{c} = numXi{c} + totXi(2:numNState, 2:numNState);
                [totGamma,eO1,eO2,totXi] = calSuffStat(seq,genPState,genFs,...
                    genBs, genScale, genEmiLik, genHMM.tr);
                genGamma = genGamma + totGamma;
                genO1 = genO1 + eO1; genO2 = genO2 + eO2; 
                genXi = genXi + totXi;
            end
        end
    end
    %disp(sprintf('nWrong = %g', nWrong))
    trainErr = trainErr / nSeqs;

function learnRt = calLearnRate(learnRtType, iter, trainErr)
    if isnumeric(learnRtType)
        learnRt = learnRtType;
    else
        switch(learnRtType)
        case 'inverse'
            learnRt = 1 / iter;
        case 'counted'
            learnRt = trainErr;
        end
    end

function genHMM = buildGeneralModel(mghmms, isJumpHmm)
    isCovFull = (ndims(mghmms{1}.sigma) == 3);
    if ~isJumpHmm
        genNState = 2; % init and terminal state
        for c = 1 : length(mghmms)
            genNState = genNState + size(mghmms{c}.tr, 1) - 2;
        end
        nDim = size(mghmms{1}.mu, 1);
        genHMM.tr = zeros(genNState, genNState);
        genHMM.mu = zeros(nDim, genNState - 1); 
        if isCovFull
            genHMM.sigma = zeros(nDim, nDim, genNState - 1);
            genHMM.sigma(:, :, 1) = eye(nDim); %mu(:,1) already zero
        else
            genHMM.sigma = zeros(nDim, genNState - 1);
            genHMM.sigma(:, 1) = ones(nDim, 1); %mu(:,1) already zero
        end
        n = 2;
        for c = 1 : length(mghmms)
            nState = size(mghmms{c}.tr, 1);
            n1 = n + nState - 3;
            genHMM.tr(1, n) = mghmms{c}.prior;
            genHMM.tr(n : n1, n : n1) = mghmms{c}.tr(2 : nState-1, 2 : nState-1);
            genHMM.tr(n1, genNState) = mghmms{c}.tr(nState - 1, nState);
            genHMM.mu(:, n : n1) = mghmms{c}.mu(:, 2 : nState-1);
            if isCovFull, genHMM.sigma(:,:,n:n1)=mghmms{c}.sigma(:,:,2:nState-1);
            else genHMM.sigma(:, n : n1) = mghmms{c}.sigma(:, 2 : nState-1);
            end
            n = n1 + 1;
        end
        genHMM.tr(genNState, genNState) = 1;
    else %JumpHmm
        genNState = 2; % init and terminal state
        for c = 1 : length(mghmms)
            genNState = genNState + size(mghmms{c}.tr, 1) - 2;
        end
        nDim = size(mghmms{1}.mu, 1);
        genHMM.tr = zeros(genNState, genNState);
        genHMM.mu = zeros(nDim, genNState - 1); 
        if isCovFull
            genHMM.sigma = zeros(nDim, nDim, genNState - 1);
            genHMM.sigma(:, :, 1) = eye(nDim); %mu(:,1) already zero
        else
            genHMM.sigma = zeros(nDim, genNState - 1);
            genHMM.sigma(:, 1) = ones(nDim, 1); %mu(:,1) already zero
        end
        n = 2;
        for c = 1 : length(mghmms)
            nState = size(mghmms{c}.tr, 1);
            n1 = n + nState - 3;
            genHMM.tr(1, n:n1) = mghmms{c}.prior * mghmms{c}.tr(1,2:nState-1);
            genHMM.tr(n:n1, n:n1) = mghmms{c}.tr(2 : nState-1, 2 : nState-1);
            genHMM.tr(n:n1, genNState) = mghmms{c}.tr(2 : nState-1, nState);
            genHMM.mu(:, n : n1) = mghmms{c}.mu(:, 2 : nState-1);
            if isCovFull, genHMM.sigma(:,:,n:n1)=mghmms{c}.sigma(:,:,2:nState-1);
            else genHMM.sigma(:, n : n1) = mghmms{c}.sigma(:, 2 : nState-1);
            end
            n = n1 + 1;
        end
        genHMM.tr(genNState, genNState) = 1;
    end

function [totGamma, eO1, eO2, totXi] = calSuffStat(seq, pState, fs, bs,scale,...
                                                   emilik, tr)
    L1 = size(emilik, 2);% length after terminal state %L = size(seq, 2);
    numStates = size(pState,1);
    if (size(fs, 2)~=L1+1), error('mmieGenTrain:FsSizeWrong', 'fs size not L+1'); end
    if (size(bs, 2)~=L1+1), error('mmieGenTrain:BsSizeWrong', 'bs size not L+1'); end
    totXi = (fs(:,1:L1) * (bs(:,2:L1+1).*emilik ...
                            ./repmat(scale(2:L1+1),numStates,1))') .* tr;
    finiteP = pState(1 : (numStates-1), 1 : size(seq, 2)); 
    finiteSeq = seq; 
    finiteSeq(~isfinite(seq)) = 0;
    totGamma = isfinite(seq) * finiteP';
    eO1 = finiteSeq * finiteP';
    eO2 = (finiteSeq.^2) * finiteP';
    if any(any(~isfinite(pState))) || any(any(~isfinite(eO1))) ...
            || any(any(~isfinite(eO2))),
        warning('mmieGenTrain:SuffStatNaN', 'Gamma/Mu/Sigma not finite');
    end
    
function mghmms = updateParameterCovFull(mghmms, numGamma, numO1, numXi, ...
    genGamma, genO1, genXi, learnRt, learnRtTr)
    smoothT = 0; smoothE = 0;
    n = 2;
    dGamma = cell(size(mghmms)); 
    dEO1 = cell(size(mghmms)); 
    dXi = cell(size(mghmms));
    genNState = size(genXi, 1);
    for c = 1 : length(mghmms)
        nState = size(mghmms{c}.tr, 1);
        n1 = n + nState - 3;
        dGamma{c} = numGamma{c} - genGamma(:,n:n1);
        dEO1{c} = numO1{c} - genO1(:,n:n1);
        %smoothE = max(smoothE, 2*max(max(genGamma(:,n:n1))));
        smoothE = max(smoothE, 2*max(genGamma(:,n:n1)));
        if (learnRtTr > eps) % Transition
            dXi{c} = numXi{c} - genXi([n:n1, genNState], [n:n1, genNState]);
            oldTr = max(mghmms{c}.tr(2:nState, 2:nState), eps);
            smoothT = max(smoothT, max(max(-dXi{c} ./ oldTr)));
            %smoothT = max(smoothT, max(-dXi{c} ./ oldTr, [], 2));
        end
        n = n1 + 1;
    end
    smoothT = smoothT * 2; %smoothE = smoothE * 2;
    %smoothT = 2 * repmat(smoothT, 1, size(mghmms{1}.tr,1) - 1); 
    smoothE = 2 * repmat(smoothE, size(genGamma,1), 1);
    w = warning('off','MATLAB:divideByZero');
    n = 2;
    for c = 1 : length(mghmms)
        nState = size(mghmms{c}.tr, 1);
        n1 = n + nState - 3;
        oldMu = mghmms{c}.mu(:, 2:nState-1); 
        %Mu = dEO1{c} ./ dGamma{c};
        Mu = (dEO1{c} + smoothE .* oldMu) ./ (dGamma{c} + smoothE);
        if any(any(~isfinite(Mu))), error('MU include NaN'); end
        mghmms{c}.mu(:, 2:nState-1) = (1-learnRt) * oldMu + learnRt * Mu; 
        if (learnRtTr > eps) % Transition
            oldTr = mghmms{c}.tr(2:nState,2:nState);
            Tr = dXi{c} + smoothT .* (oldTr > eps); 
            Tr = Tr ./ repmat(sum(Tr,2), 1, nState-1);
            if any(any(~isfinite(Tr))), 
                warning('mmieGenTrain:TrNaN', 'TR include NaN'); 
            else
                mghmms{c}.tr(2:nState, 2:nState) = (1-learnRtTr)*oldTr+learnRtTr*Tr; 
            end
        end
        n = n1 + 1;
    end
    warning(w);

function mghmms = updateParameter(mghmms, numGamma, numO1, numO2, numXi, ...
    genGamma, genO1, genO2, genXi, learnRt, learnRtTr)
    smoothE = 0; smoothT = 0;
    n = 2;
    mu2 = cell(size(mghmms)); dGamma = cell(size(mghmms)); 
    dEO1 = cell(size(mghmms)); dEO2 = cell(size(mghmms));
    dXi = cell(size(mghmms));
    genNState = size(genXi, 1);
    for c = 1 : length(mghmms)
        nState = size(mghmms{c}.tr, 1);
        n1 = n + nState - 3;
        % 2*genGamma Povey
        mu = mghmms{c}.mu(:, 2:nState-1);
        Sigma = mghmms{c}.sigma(:, 2:nState-1); %variance not std deviation
        if any(any(Sigma <= 0)), error('SIGMA negative'); end
        mu2{c} = mu .^ 2;
        dGamma{c} = numGamma{c} - genGamma(:,n:n1);
        dEO1{c} = numO1{c} - genO1(:,n:n1); dEO2{c} = numO2{c} - genO2(:,n:n1);
        tmp = 2 * mu .* dEO1{c} - dEO2{c} - (Sigma + mu2{c}) .* dGamma{c};
        % Twice smallest, so /2 canceled out
        smoothE1 = max(max((tmp + sqrt(tmp.^2 - 4 * Sigma.^2 ...
            .* (dGamma{c}.*dEO2{c}-dEO1{c}.^2))) ./ Sigma));
        smoothE = max(smoothE, max(2*max(max(genGamma(:,n:n1))), smoothE1));
        if (learnRtTr > eps) % Transition
            dXi{c} = numXi{c} - genXi([n:n1, genNState], [n:n1, genNState]);
            oldTr = max(mghmms{c}.tr(2:nState, 2:nState), eps);
            smoothT = max(smoothT, max(max(-dXi{c} ./ oldTr)));
            %smoothT = max(smoothT, 4 * max(max(dXi{c} ./ mghmms{c}.tr(2:nState, 2:nState))));
        end
        n = n1 + 1;
    end
    smoothE = smoothE * 2; smoothT = smoothT * 2;
    w = warning('off','MATLAB:divideByZero');
    n = 2;
    for c = 1 : length(mghmms)
        nState = size(mghmms{c}.tr, 1);
        n1 = n + nState - 3;
        denominator = dGamma{c} + smoothE;
        oldMu = mghmms{c}.mu(:, 2:nState-1); 
        Mu = (dEO1{c} + smoothE .* oldMu) ./ denominator;
        if any(any(~isfinite(Mu))), error('MU include NaN'); end
        mghmms{c}.mu(:, 2:nState-1) = (1-learnRt) * oldMu + learnRt * Mu; 
        oldSigma = mghmms{c}.sigma(:, 2:nState-1);
        Sigma = (dEO2{c} + smoothE.*(oldSigma+mu2{c})) ./ denominator - Mu.^2;
        if any(any(~isfinite(Sigma))), error('SIGMA include NaN'); end
        if any(any(Sigma < 0)), disp(Sigma); error('SIGMA negative'); end
        mghmms{c}.sigma(:,2:nState-1) = (1-learnRt)*oldSigma + learnRt*Sigma;
        if (learnRtTr > eps) % Transition
            oldTr = mghmms{c}.tr(2:nState,2:nState);
            Tr = dXi{c} + smoothT * (oldTr > eps); 
            Tr = Tr ./ repmat(sum(Tr,2), 1, nState-1);
            if any(any(~isfinite(Tr))), 
                warning('mmieGenTrain:TrNaN', 'TR include NaN'); 
            else
                mghmms{c}.tr(2:nState, 2:nState) = (1-learnRtTr)*oldTr+learnRtTr*Tr; 
            end
        end
        n = n1 + 1;
    end
    warning(w);
    
function convergent = isConvergent(m1, m0, logCondLiks, oldLogCondLiks, tol)
    % Check convergence on loglik, TR and E (options trtol and etol).
    isCovFull = (ndims(m0{1}.sigma) == 3);
    convergent = false;
    if (logCondLiks - oldLogCondLiks) / (abs(oldLogCondLiks) + 1) < tol.lcl
        convergent = true;
        for c = 1 : length(m1)
            if isCovFull
                if max(max(abs(m1{c}.mu - m0{c}.mu))) >= tol.mu || ...
                    max(max(max(abs(sqrt(m1{c}.sigma)-sqrt(m0{c}.sigma)))))...
                    >=tol.sigma,
                    convergent = false;
                    break; %for
                end
            else
                if max(max(abs(sqrt(m1{c}.sigma)-sqrt(m0{c}.sigma))))>=tol.sigma...
                    || max(max(abs(m1{c}.mu - m0{c}.mu))) >= tol.mu,
                    convergent = false;
                    break; %for
                end
            end
        end
    end

function verboseOutput(iter, logCondLiks)
    fprintf('It %d: logCondLik=%.3g\n', iter, logCondLiks);
%     if iter == 1
%         fprintf('Relative Changes in Log Likelihood\n');
%     else
%         fprintf('It %d: logCLik %.3f\n', iter, loglik);
%     end

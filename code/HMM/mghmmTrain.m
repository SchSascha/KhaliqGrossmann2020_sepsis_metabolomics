function [guessTR,guessMu,guessSigma,logliks] = mghmmTrain(seqs,...
    guessTR,guessMu,guessSigma,varargin)
%MGHMMTRAIN trains a multivariate Gaussian emissian HMM using EM algorithm.
%   [GUESSTR, GUESSMU, GUESSSIGMA, LOGLIKS] = MGHMMTRAIN(SEQS, GUESSTR,
%   GUESSMU, GUESSSIGMA) estimates the transition and emission
%   probabilities for a multivariate Gaussian-emission Hidden Markov Model
%   from the training sequences using the Baum-Welch algorithm. The training
%   sequences, SEQS, can be a single sequence or a cell array of sequences.
%   The input arguments GUESSTR, GUESSMU, GUESSSIGMA are initial estimates
%   of the transition and emission probability, and the corresponding
%   output arguments are the optimised parameters with maximal likelihood.
%   SEQS are sequences (time series) of real numbers (NaN denotes missing
%   values), SEQS{n}(g,t) or SEQS(n,g,t) denoting the observed value at time t
%   of sample n. GUESSTR(i,j) is the transition probability from state i to j, 
%   and GUESSMU(g,i) and GUESSSIGMA(g,i) are mean and standard deviation of
%   gene g and state i. 
%
%   A hidden state and corresponding emission (observation) is added before
%   the first observation, at time zero (t=0). The state at time zero is
%   assumed to be the first state, so the emission probability from the
%   first to all other states is equivalent to the initial probability. 

%   This script is modified from MathWorks hmmtrain.m.

% default values
tol = 1e-6;
trtol = tol;
mutol = tol;
sigmatol = tol;
maxiter = 100; %500;
pseudoTRcounts = false;
verbose = false;
COV_DIAG = 0; COV_FULL = 1; covType = COV_DIAG; 
glbCovRt = 0; glbCov = [];
MISS_EM_INDEP = 0; MISS_EM_COND = 1;
missFillType = MISS_EM_INDEP;

% verify the input
[numStates, checkTr] = size(guessTR);
if checkTr ~= numStates, error('TRANSITION matrix must be square.'); end
if numStates==0, error('Number of states must be positive.'); end

% options
if nargin > 4
    if rem(nargin,2) == 1, 
        error('Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'tolerance','pseudotransitions','maxiterations','verbose',...
              'trtol','mutol','sigmatol','globalmeanvar','cov','rtglbcov',...
              'glbcov','miss'};
    for j=1:2:nargin-4
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k), continue; %error('Unknown parameter name:  %s.',pname);
        elseif length(k)>1, error('Ambiguous parameter name:  %s.',pname);
        end
        switch(k)
        case 1  % tolerance
            tol = pval;
            trtol = tol;
            mutol = tol; sigmatol = tol;
        case 2  % pseudocount transitions
            pseudoTR = pval;
            [rows, cols] = size(pseudoTR);
            if rows ~= cols, error('PSEUDOTR matrix must be square.'); end
            if  rows ~= numStates
                error('PSEUDOTR must have the same size with TRANSITION.');
            end
            pseudoTRcounts = true;
        case 3 % max iterations
            maxiter = pval;
        case 4 % verbose
            if islogical(pval) || isnumeric(pval), verbose = pval; end
        case 5 % transtion tolerance
            trtol = pval;
        case 6 % mu tolerance
            mutol = pval;
        case 7 % sigma tolerance
            sigmatol = pval;
        case 8 % global mean & var
%             if isnumeric(pval) && numel(pval)==2
%                 globalVar = pval(2); %globalMean = pval(1); 
%             else error('GLOBALMEANVAR must have two elements'); 
%             end
        case 9 % cov
            if     strcmpi('diag',pval), covType = COV_DIAG; 
            elseif strcmpi('full',pval), covType = COV_FULL; 
            end
        case 10 %rtglbcov
            glbCovRt = pval;
        case 11 %glbcov
            glbCov = pval;
        case 12 %miss
            if     strcmpi('em_indep',pval), missFillType = MISS_EM_INDEP; 
            elseif strcmpi('em_cond' ,pval), missFillType = MISS_EM_COND; 
            end
        end
    end
end

[nDim, statesMu] = size(guessMu);
if (statesMu ~= numStates-1)
    error('MU must have the same number of columns as TRANSITION.');
end
if covType == COV_FULL
    if (ndims(guessSigma) ~= 3), error('Full SIGMA must be 3D'); end
    [dimSigma, dimSigma2, statesSigma] = size(guessSigma);
    if (dimSigma ~= dimSigma2), error('Full SIGMA must be square'); end
else
    [dimSigma, statesSigma] = size(guessSigma);
end
if (statesSigma ~= numStates-1)
    error('SIGMA matrix must have the same number of columns as TRANSITION.');
elseif (nDim ~= dimSigma)
    error('SIGMA matrix must have the same number of rows (dimensions).');
end

if isnumeric(seqs)
    [numSeqs, L] = size(seqs);
    cellflag = false;
elseif iscell(seqs)
    numSeqs = numel(seqs);
    cellflag = true;
else
    error('SEQS must be cell array or numerical array.');
end

if ~pseudoTRcounts, pseudoTR = zeros(size(guessTR)); end

if covType == COV_FULL
    SMALL_DIAG_COV = repmat(1e-2 * eye(nDim), [1,1,numStates-1]); 
    NEG_DIAG_COV = -SMALL_DIAG_COV;
    SMALL_DIAG_COV(SMALL_DIAG_COV == 0) = -Inf; 
    NEG_DIAG_COV(NEG_DIAG_COV == 0) = -Inf; 
    
    if isempty(glbCov), glbCovRt = 0; end
end

%-----------------------------------------------------------
% Core algorithm implementation
%-----------------------------------------------------------

converged = false;
loglik = 1; % make sure not to converge after first step
logliks = zeros(1,maxiter);
for iteration = 1:maxiter    
    oldLL = loglik;
    loglik = 0;
    oldGuessMu = guessMu;
    oldGuessSigma = guessSigma;
    oldGuessTR = guessTR;
    TR = pseudoTR;
    Mu = zeros(nDim, numStates-1);
    if covType ~= COV_FULL
        totalPost = zeros(nDim, numStates-1);
        y2 = zeros(nDim, numStates-1);
    elseif glbCovRt < 1
        totalPost = zeros(numStates-1, 1);
        covY = zeros(nDim, nDim, numStates-1);
    end
    for count = 1:numSeqs
        if cellflag
            seq = seqs{count}; L = size(seq, 2);
        else
            seq = squeeze(seqs(count,:,:));
        end
        emilik = mghmmEmiLik(seq, guessMu, guessSigma, varargin{:});
        if any(any(~isfinite(emilik))), disp(emilik); error('emilik NaN'); end
        L1 = size(emilik, 2); %length after adding terminal state
        % E-step (get the scaled forward and backward probabilities)
        [p,logPseq,fs,bs,scale] = hmmDecodeLik(guessTR, emilik);
        % if seq likelihood is 0, skip this seq, following HMMER implementation
        if size(p, 1) == 0, continue; end %seq likelihood zero
        if ~all(all(isfinite(fs))) || ~all(all(isfinite(bs))), %here for dbstop
            warning('mghmmTrain:FsBsNaN', 'fwd/bck matrix has NaN'); 
        end
        loglik = loglik + logPseq;
        % M-step (preparations)
        % Rabiner: xi_t(i,j) = alpha_t(i) a(i,j) b_t(j) beta_t+1(j)
        isFin = isfinite(seq);
        scale1 = repmat(scale(2:L1+1),numStates,1);
        TR = TR + (fs(:,1:L1) * (bs(:,2:L1+1).*emilik./scale1)') .* guessTR;
        if covType ~= COV_FULL 
            finiteP = p(1:numStates-1,1:L); 
            finiteSeq = seq; 
            finiteSeq(~isFin) = 0;
            Mu = Mu + finiteSeq * finiteP';
            totalPost = totalPost + isFin * finiteP';
            y2 = y2 + (finiteSeq.^2) * finiteP';
        elseif glbCovRt < 1
            totalPost = totalPost + sum(p(1:numStates-1,1:L), 2);
            for t = 1 : L  %before terminal state
                % Fill in missing values
                ot = isFin(:,t); ut = ~ot; 
                Yut = zeros(sum(ut),1); Yot = seq(ot,t);
                if missFillType ~= MISS_EM_COND
                    for j=2:numStates-1, Yut = Yut + p(j,t) * guessMu(ut,j); end
                else
                    for j = 2 : numStates-1
                        Yut = Yut + p(j,t)*(guessMu(ut,j)+guessSigma(ut,ot,j)...
                           *inv(guessSigma(ot,ot,j))*(Yot-guessMu(ot,j)));
                    end
                end
                Yt = zeros(nDim, 1); Yt(ut) = Yut; Yt(ot) = Yot;
                % Count sufficient statistics
                for j = 2 : numStates-1
                    covY(:,:,j) = covY(:,:,j) + p(j,t) * Yt * Yt';
                    Mu(:,j) = Mu(:,j) + p(j,t) * Yt;
                end
            end
        end
        if any(any(~isfinite(p))) || any(any(~isfinite(TR))) ...
                || any(any(~isfinite(Mu))) || ~isfinite(loglik)
            warning('ghmmTrain:PostNaN', ...
                'Posterior/Tr/Mu/Loglik include NaN'); 
        end
    end
    w = warning('off','MATLAB:divideByZero');
    totalTR = sum(TR,2);
    guessTR  = TR ./ repmat(totalTR,1,numStates);
    if any(totalTR == 0)
        noTransitionRows = find(totalTR == 0);
        guessTR(noTransitionRows,:) = 0;
        guessTR(sub2ind(size(guessTR),noTransitionRows,noTransitionRows)) = 1;
    end
    guessTR(isnan(guessTR)) = 0; % clean up any remaining NaN
    if covType ~= COV_FULL 
        guessMu = Mu ./ totalPost; 
        if any(any(totalPost == 0)), guessMu(totalPost == 0) = 0; end
        if any(any(~isfinite(guessMu))), 
            warning('guessMu include NaN (iteration %d)', iteration); 
            logliks(logliks ==0) = []; return; 
        end
        guessSigma = y2 ./ totalPost - guessMu.^2 ;
        if any(guessSigma<-1e-3), disp(guessSigma);error('Sigma include negative'); end
        guessSigma(guessSigma <= 1e-3) = 1e-3; 
        if any(any(totalPost == 0)), guessSigma(totalPost == 0) = 1; end
        if any(any(~isfinite(guessSigma))), error('guessSigma include NaN'); end
    elseif glbCovRt < 1
        for j = 2 : numStates-1
            if totalPost(j) == 0, continue; end
            guessMu(:,j) = Mu(:,j) ./ totalPost(j); 
            S = glbCovRt * glbCov + (1 - glbCovRt) ...
                * (covY(:,:,j)./totalPost(j) - guessMu(:,j)*guessMu(:,j)');
            S = (S + S') / 2; 
            [upper, part] = chol(S); %check positive definite
            if part == 0, guessSigma(:,:,j) = S; end
        end
        if any(totalPost(2:numStates-1) < 0), disp(totalPost); error('sum_t gamma_t(j) < 0'); end
        if any(any(any(guessSigma<NEG_DIAG_COV))), error('Sigma diagonal include negative'); end
        guessSigma = max(guessSigma, SMALL_DIAG_COV);
        if any(any(any(~isfinite(guessSigma)))), error('guessSigma include NaN'); end
    end
    warning(w);

    if verbose
        if iteration == 1
            fprintf(['Relative Changes in Log Likelihood, Transition Matrix'...
                ' and Emission Matrix\n']);
        else
            fprintf('Iter %d: loglik %.3f  dTR: %.3f  dMU: %.3f  dSIGMA: %.3f\n',...
                iteration, loglik,...
                norm(guessTR-oldGuessTR,inf)/numStates,...
                max(max(abs(guessMu-oldGuessMu))),...
                max(max(abs(sqrt(guessSigma)-sqrt(oldGuessSigma)))));
                %(abs(loglik-oldLL)/(1+abs(oldLL))), 
        end
    end
    % Check convergence on loglik, TR and E (options trtol and etol).
    logliks(iteration) = loglik;
    if (abs(loglik-oldLL)/(1+abs(oldLL))) < tol
        if norm(guessTR - oldGuessTR,inf)/numStates < trtol
            if covType == COV_FULL
                if max(max(abs(guessMu-oldGuessMu))) < mutol && max(max(max(abs(...
                        sqrt(guessSigma)-sqrt(oldGuessSigma))))) < sigmatol,
                    if verbose
                        fprintf('Algorithm converged after %d iterations.',iteration);
                    end
                    converged = true;
                    break
                end
            else
                if max(max(abs(guessMu-oldGuessMu))) < mutol && max(max(abs(...
                        sqrt(guessSigma)-sqrt(oldGuessSigma)))) < sigmatol,
                    if verbose
                        fprintf('Algorithm converged after %d iterations.',iteration);
                    end
                    converged = true;
                    break
                end
            end
        end
    end
end
if ~converged
    %warning('stats:hmmtrain:NoConvergence',...
    if verbose
        fprintf('Algorithm did not converge with tolerance %f in %d iterations.\n',...
            tol,maxiter);
    end
end
logliks(logliks ==0) = [];

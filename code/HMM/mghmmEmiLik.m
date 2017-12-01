function [emilik, emilik3d] = mghmmEmiLik(seq, mu, sigma, varargin)
%MGHMMEMILIK calculates the emission likelihood for a multivariate
%Gaussian-emission HMM. 

[isTermState, arg] = varArgRemove('termstate', 0, varargin);
[covTypeStr, arg] = varArgRemove('cov', '', arg);
COV_DIAG = 0; COV_FULL = 1; covType = COV_DIAG; 
if strcmpi('full',covTypeStr), covType = COV_FULL; end

numStates = size(mu, 2) + 1;
nDim = size(mu, 1);
L = size(seq, 2);
% For 'Maximum number of users for Statistics_Toolbox reached'
% for k = 1 : 10 
%     try
%         emilik3d = normpdf(permute(repmat(seq,[1,1,numStates-1]),[1,3,2]), ...
%                    repmat(mu,[1,1,L]), repmat(sqrt(sigma),[1,1,L]));
%         break;
%     catch
%         [errMsg, errId] = lasterr;
%         if strfind(lower(errMsg), 'license'), 
%             if k<10, pause(60); else error(errId,[errMsg '\n10 times\n']); end
%         else error(errId, errMsg); break; end
%     end
% end
if covType == COV_FULL
    emilik3d = []; emilik = zeros(numStates-1, L);
    for t = 1 : L
        x = seq(:,t)'; isFinX = isfinite(x); x = x(isFinX);
        if isempty(x), 
            emilik(:,t) = 1; 
        else 
            for j = 2 : numStates-1
                emilik(j,t) = mymvnpdf(x, mu(isFinX,j)', sigma(isFinX,isFinX,j));
            end
        end
    end
else
    mu = mu(:, 2 : numStates-1); sigma = sigma(:, 2 : numStates-1);
    emilik3d = zeros(nDim, numStates-1, L);
    %emilik3d(:, 1, :) = zeros(nDim, L);
    emilik3d(:, 2 : numStates-1, :) = mynormpdf(...
               permute(repmat(seq,[1,1,numStates-2]),[1,3,2]), ...
               repmat(mu,[1,1,L]), repmat(sqrt(sigma),[1,1,L]));
    emilik3d(~isfinite(emilik3d)) = 1;
    emilik = squeeze(prod(emilik3d, 1));
end
if isTermState,
    emilik = [emilik, zeros(numStates-1,1); zeros(1,L), 1];
else
    emilik = [emilik; zeros(1,L)];
end

function y = mynormpdf(x,mu,sigma)
    y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

% function y = mymvnpdf(X, Mu, Sigma)
%     % Make sure Sigma is a valid covariance matrix
%     [spd,R] = isspd(Sigma);
%     if spd
%         % Create array of standardized data, vector of inverse det
%         xRinv = X0 / R;
%         sqrtInvDetSigma = 1 / prod(diag(R));
%     else
%         error('stats:mvnpdf:BadSigma',...
%               'SIGMA must be symmetric and positive definite.');
%     end
%     quadform = sum(xRinv.^2, 2);
%     y = sqrt((2*pi)^(-d)) * sqrtInvDetSigma .* exp(-0.5*quadform);

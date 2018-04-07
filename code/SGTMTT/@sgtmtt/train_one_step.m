function [g]=train_one_step(g,iter)
% Train the GTM-TT model
% modifications: Frank-Michael Schleif, March-June, 2011
if ~strcmp(g.Status,'initialized')
    if g.Verb
        disp('MSG: Initializing the model...');
    end
    g = init(g);
end

FI         = g.FI;
K          = g.K;
FI_T       = FI';
N          = g.Data.N;
NFullData  = sum(N);
XMatFormat = cell2mat(g.Data.X');
Deff       = g.Par.Deff;

%% Training procedure 

% Forward/backward procedure (E-step)        

% calculate emission probabilities
dBeta = g.Par.Beta;
BjnMat = ((dBeta/(2*pi))^(Deff/2)*exp(-0.5*dBeta*g.DistXY));
Bjn = mat2cell(BjnMat,K,N);
    
% HMM forward-backward algorithm
% Gamma - resp. matrix, gamma-init - for PI reestimate , Xi - expected number of transitions
[Gamma, GammaInit, Xi, LLhood, LLhoodV]=forwback(g.Par.PI, g.Par.A, Bjn, N);
g.RespMat = Gamma;    
if g.Verb
    fprintf('It: %g\tlogLH: %g\tBeta: %g\n', iter, LLhood, g.Par.Beta);
end
g.LLhood(iter) = LLhood;
g.LLhoodxObs(iter,:) = LLhoodV;
    
%M-step    
g.Par.A = Xi./repmat(sum(Xi,2),[1 K]);    % re-estimate of the transition prob matrix (just normalize Xi)
g.Par.A(isnan(g.Par.A))=0;                  % cleanup
g.Par.PI = GammaInit/sum(GammaInit);        % update priors (just normalize gammainit)   
    
%Resolve equation Phi'*Gs*Phi*W' = Phi' * Rs *X
% Rs - the responsibility matrix r_kn = p_nk' / sum p_nk'
% Gs - diag(sum(Rs))
% phis (M+1 x K) * Gs (K x T x N) * phis (K x M+1) * weights (D x M+1) == phis (M+1 x K) * Rs * Data (T x N x D)
%Update W
GammaMat = cell2mat(Gamma);

%  Phi'*G*Phi 
A1 = full(FI_T*spdiags(sum(GammaMat,2), 0, K, K)*FI);
% get upper triangle matrix R such that R * R' = A1
[cholDcmp singular] = chol(A1);

if singular
    if g.Verb
        fprintf('WARNING: M-Step matrix singular, using pinv.\n');
    end
    pInvAI  = pinv(A1);
    g.Par.W = pInvAI*(FI_T*(GammaMat*XMatFormat));
else
    % R' * R * W' = Phi' * Rs *X  
    g.Par.W = cholDcmp \ (cholDcmp' \ (FI_T*(GammaMat*XMatFormat)));
end

% Calculate Y new
g.Par.Y    = FI*g.Par.W;  % y(x_k, W^), latend points mapped to the orig D space using basis functions FI and updated weights    
Y          = g.Par.Y;
X          = XMatFormat;
[g.DistXY] = g.Par.oMetric.distance(Y, X); % includes relevances here
g.Par.Beta = NFullData*Deff/sum(sum(g.DistXY.*GammaMat)); 

function [Gamma, GammaInit, Xi, LLhood, LLhoodV,LhoodV]=forwback_stable(PI, A, Bjn, N)
% Forward/Backward procedure
% for standard GTMTT class
%
% Input:
%   PI  - priors
%   A   - transition matrix
%   Bjn - emission prob. for each item at each time point with respect to each state / latent point (phi) 
%   N   - vector of time points per sample
%
% Output:
%   Xi         - expected number of transitions (K x K) - for A reestimate
%   Gamma      - responsibility matrix K x N
%   Gamma_init - for PI reestimate
%   LLhood     - log likelihood overall
%   LLhoodV    - llh per sample
%   LhoodV     - maximum likelihood of observing the sequence (no log, no sum)
%
% Ivan Olier, May 18th 2007
% modifications: Frank-Michael Schleif, March, 2011 
% also to achieve better numerical stability

if ~iscell(Bjn)
    BjnMat = Bjn;
    clear Bjn;
    Bjn{1} = BjnMat;
end

K         = size(Bjn{1},1);
NrObs     = length(N);
Xi        = zeros(K);
LLhoodV   = zeros(1,NrObs);
GammaInit = zeros(K,1);
LhoodV    = zeros(NrObs,N(1));
Gamma     = cell(1,NrObs);


for Obs = 1:NrObs
    NObs = N(Obs);            % time points per sample
    Alpha = zeros(K,NObs);    % forward probabilities, prob of sequence j from time 1 to t in state Xi at time t
    Beta  = ones(K,NObs);     % backward probabilities,  
    % last time step in last observation using forward prob 1/P(Y|\pars) - to literature this is not scaled  - thats why I init to ones    
    c_factor = zeros(1,NObs);
    BjnCurrentObs = Bjn{Obs};   % B is the matrix of the probability of events/emissions - event/emission matrix B    
    % probability to observe an emission in state pi
    Alpha(:,1) = log(max(PI,realmin)) + log(max(BjnCurrentObs(:,1),realmin)); % first step of the inductive forward estimation for j in time point 1 with priors PI for all states Xi
    shift      = max(Alpha(:,1)); % avoid overflow of exp
    Alpha(:,1) = exp(Alpha(:,1)-shift);
    
    %scaling factor...
    c_factor(1) = sum(Alpha(:,1));  % for numerical stability of the log-likelihoods - to sum to 1
    Alpha(:,1) = Alpha(:,1)/c_factor(1); % renormalization
    c_factor(1) = log(max(c_factor(1), realmin)) + shift;
    
    for n=2:NObs % A is the transition probability matrix, Alpha the forward prob matrix
        Alpha(:,n)  = log(max(A'*Alpha(:,n-1), realmin)) + log(max(BjnCurrentObs(:,n),realmin)); % recursive estimating step - propagating through time
        shift       = max(Alpha(:,n));
        Alpha(:,n)  = exp(Alpha(:,n) - shift);
        c_factor(n) = sum(Alpha(:,n));
        Alpha(:,n) = Alpha(:,n)/c_factor(n); % renormalization
        c_factor(n) = log(max(c_factor(n), realmin)) + shift;
    end

    %Calculates beta_n -backward, prob of sequence j from time t+1 to T in state Xi at time t
    % probability to generate no further symbols at T from any state is 1
    for n=NObs-1:-1:1 % go back in time, 
        % A is the transition probability matrix, Beta the backward prob matrix
        tmp       = log(max(Beta(:,n+1),realmin)) + log(max(BjnCurrentObs(:,n+1),realmin));
        shift     = max(tmp);
        Beta(:,n) = log(max(A* exp(tmp - shift), realmin)); 
        Beta(:,n) = exp(Beta(:,n) - c_factor(n+1) + shift);
    end

    % Normalization factors in acc. to Strachan, multiply 4 probs (p.112)
    % 1) being in state Xi at t given obs_j till then (forward only)
    % 2) Transition prob Xi->Xk via A (transition prob. matrix)
    % 3) Prob. being in Xk given obs y at t+1 (using emission prob B)
    % 4) Prob. of partial seq from t+1 to T in Xk (via backward)
    %%%
    %% - this would be (theoretically) the same but is numerical less stable
    %% Xi    = A .* (Alpha(:,1:end-1) * (Beta(:,[2:end]).*BjnCurrentObs(:,[2:end]))');
    %%
    % should be stable enough now
    % Xi - expected number of transitions from Xi to Xk
    for n=1:NObs-1
            Xi=Xi+ (A.*(Alpha(:,n)*(Beta(:,n+1).*BjnCurrentObs(:,n+1))'));    
    end

    Rkn = Alpha.*Beta;
    Rkn = bsxfun(@rdivide,Rkn,max(sum(Rkn,1),realmin));
    LLhoodV(Obs)  = sum(c_factor);        % loglikelihood
    LhoodV(Obs,:) = exp(c_factor)';            % maximum likelihood of observing such a sequence by the model 
    Gamma{Obs} = Rkn;                    % assign resp. for each observation
    GammaInit = GammaInit + Rkn(:,1);    % acc. responsibilities for being in state Xi at t=1
end

LLhood = sum(LLhoodV);
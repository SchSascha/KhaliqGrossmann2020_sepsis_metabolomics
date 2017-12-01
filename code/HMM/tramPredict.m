function logOdds = tramPredict(data, model, varargin)
%TRAMPREDICT predicts the classes of testing data based on the hidden
%Markov models (HMM).
%   LOGODDS = TRAMPREDICT(DATA, MODEL) returns the log odds of each
%   instance in DATA being a positive example, based on the HMMs given by
%   MODEL. DATA{k} are the time series gene expresion of patient k.
%   MODEL{1} and MODEL{2} are the HMMs for negative and positive examples,
%   respectively. See TRAMGENTRAIN for the data structure. LOGODDS(k) are
%   the log odds of patient k being a positive example.

lpr = log(model{2}.prior) - log(model{1}.prior);
logOdds = zeros(length(data),1);
for q = 1 : length(data)
    [pStates, loglik2] = mghmmDecode(data{q}, model{2}.tr, model{2}.mu, ...
                                     model{2}.sigma, varargin{:});
    [pStates, loglik1] = mghmmDecode(data{q}, model{1}.tr, model{1}.mu, ...
                                     model{1}.sigma, varargin{:});
    logOdds(q) = loglik2 - loglik1 + lpr;
end

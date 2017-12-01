function discModel = tramDiscTrain(trainData, model, varargin)
%TRAMDISCTRAIN discriminatively trains hidden Markov models (HMM) for the
%classification of time series gene expression data.
%   DISCMODEL = TRAMDISCTRAIN(TRAINDATA, MODEL) optimize the HMMs given by
%   MODELS, in terms of the conditional likelihood of the training gene
%   expression data, TRAINDATA. TRAINDATA{1} are the negative training
%   examples and TRAINDATA{2} are the positive training examples. MODEL{1}
%   and MODEL{2} are the HMMs for negative and positive examples,
%   respectively. See TRAMGENTRAIN for the data structure. 
%
%   TRAMDISCTRAIN(...,'mmieIterations',MAXITER) performs the Extended
%   Baum-Welch algorithm for maximally MAXITER iterations.

nState = size(model{1}.mu, 2) - 1;
discModel = mmieGenTrain(trainData, model, 'nState', nState, varargin{:});

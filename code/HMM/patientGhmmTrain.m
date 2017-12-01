function [patientGhmm, logLikArr] = patientGhmmTrain(trainX, varargin)
%PATIENTGHMMTRAIN trains one gene-specific GHMM for each category
%   [PATIENTGHMM] = PATIENTGHMMTRAIN(TRAINX, ...)
%   trains one gene-specific GHMM for each category. TRAINX{c}{q}(g,t)
%   denotes the gene expression level at time t of patient q and gene g of
%   category c. 

nCate = length(trainX);
patientGhmm = cell(nCate, 1);
logLikArr = cell(nCate, 1);
nPatients = 0;
for c = 1 : nCate,
    [patientGhmm{c} logLikArr{c}] = trainCategory(trainX{c}, varargin{:});
    nPatients = nPatients + patientGhmm{c}.prior;
end
for c = 1 : nCate, patientGhmm{c}.prior = patientGhmm{c}.prior / nPatients; end

function [patientGhmm, loglikArr] = trainCategory(cateX, varargin)
    patientGhmm.prior = size(cateX,1);
    [tr1 mu1 sigma1 loglikArr] = mghmmTrainNew(cateX, varargin{:});
    patientGhmm.tr = tr1; patientGhmm.mu = mu1; patientGhmm.sigma = sigma1;
%     for iRepeat = 1 : length(loglikArr), plot(loglikArr{iRepeat}); end

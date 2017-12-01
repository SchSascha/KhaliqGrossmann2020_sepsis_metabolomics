function model = tramGenTrain(trainData, nState, nFeatArr, varargin)
%TRAMGENTRAIN generatively trains hidden Markov models (HMM) for the
%classification of time series gene expression data.
%   MODEL = TRAMGENTRAIN(TRAINDATA, NSTATE, NFEATARR) generatively trains
%   HMMs of NSTATE states, one for each class of time series expression in
%   TRAINDATA, and performed gene selection of gene numbers in array
%   NFEATARR. TRAINDATA{1} are the negative training examples and
%   TRAINDATA{2} are the positive training examples. NSTATE is the number
%   of states. NFEATARR is an array of gene numbers, starting with total
%   number of genes. This function will go through NFEATARR to find an
%   optimal gene number. MODEL is the trained HMMs; MODEL{1} is the HMM for
%   negative examples and MODEL{2} is the HMM for positive examples.
%   MODEL{1}.selGenes are the selected genes; MODEL{1}.prior is the prior
%   probability of the negative class; MODEL{1}.tr are the transition
%   probabilities of the HMM; MODEL{1}.mu and MODEL{1}.sigma are the mean
%   and variance of the emission probability of the HMM. MODEL{2} has the
%   same structure.
%
%   TRAMGENTRAIN(...,'replicates',REPLICATES) performs the Baum-Welch
%   algorithm from different random initializations. REPLICATES is the
%   number of initializations.
%
%   TRAMGENTRAIN(...,'maxIterations',MAXITER) performs the Baum-Welch
%   algorithm for maximally MAXITER iterations.
%
%   TRAMGENTRAIN(...,'validateFold',MAXITER) specifies the fold number in
%   internal cross validation for optimal gene number. The default value is
%   10 folds.

[validateFold, algArg] = varArgRemove('validatefold', 10, varargin);
[mmieIterations, algArg] = varArgRemove('mmieIterations', 0, algArg);

isRFE = true;
nGene = size(trainData{1}{1}, 1);
nCate = length(trainData);
nFeatArr = sort(nFeatArr, 'descend');
featModels = cell(nCate, length(nFeatArr));
selGenes = 1 : nGene;

% Initial training with all features (genes)
mghmms = mmieGenTrainNew(trainData, nState, 'mmieIterations', 0, algArg{:});

% Feature ranking (gene ranking)
if ~isRFE && ~(length(nFeatArr) == 1 && nFeatArr(1) == nGene)
    geneOrder = rankFeature(mghmms, trainData, algArg{:});
end

% Evaluation
if length(nFeatArr) < 1, error('length(nFeatArr) < 1'); end
modSelAcc = zeros(length(nFeatArr), 1);
for iFeatArr = 1 : length(nFeatArr)
    nSelGene = nFeatArr(iFeatArr);
    if nSelGene == nGene
        mghmmFeatSel = mghmms;
        if isRFE, featSelTrainX = trainData; end % for next re-rank
        if validateFold > 0
            modSelAcc(iFeatArr) = calModelSelAcc(trainData, validateFold,...
                                  nState, 'mmieIterations', 0, algArg{:});
            %modSelAcc(iFeatArr) = parCalModelSelAcc(trainData, validateFold,...
            %                      nState, 'mmieIterations', 0, algArg{:});
        end
    else
        if isRFE, 
            geneOrder = rankFeature(mghmmFeatSel, featSelTrainX, algArg{:});
            selGenes = selGenes(geneOrder(1:nSelGene));
            featSelTrainX = selectFeature(featSelTrainX, geneOrder(1:nSelGene));
        else
            selGenes = geneOrder(1:nSelGene);
            featSelTrainX = selectFeature(trainData, selGenes);
        end
        if validateFold > 0
            modSelAcc(iFeatArr) = calModelSelAcc(featSelTrainX, validateFold,...
                                  nState, 'mmieIterations', 0, algArg{:});
            %modSelAcc(iFeatArr) = parCalModelSelAcc(featSelTrainX, validateFold,...
            %                      nState, 'mmieIterations', 0, algArg{:});
                              
        end
        mghmmFeatSel = mmieGenTrainNew(featSelTrainX, nState, ...
                                       'mmieIterations', 0, algArg{:});
    end
    mghmmFeatSel{1}.selGenes = selGenes; mghmmFeatSel{2}.selGenes = selGenes;
    featModels(:, iFeatArr) = mghmmFeatSel;
end
[topMsAcc, iFeatTop] = max(modSelAcc);
model = featModels(:, iFeatTop);
model.acc = topMsAcc;

%-----------------------------------------------------------Subroutines

function geneOrder = rankFeature(mghmms, trainData, varargin)
    trainLLR2 = patientGhmmGeneLLR(mghmms, trainData{2}, varargin{:});
    trainLLR1 = patientGhmmGeneLLR(mghmms, trainData{1}, varargin{:});%1 x G
    trainLLR = trainLLR2 - trainLLR1;%1 x G
    [sortTrainLLR, geneOrder] = sort(trainLLR, 'descend');

function msAcc1 = calModelSelAcc(fsTrnX, validateFold, nState, varargin)
    foldTrnX = cell(2, 1); foldValX = cell(2, 1);
    fsAcc = 0;
    logOdds = cell(1, 2);
    for fold = 1 : validateFold
        for m = 1 : 2
            q1 = floor(length(fsTrnX{m}) * (fold-1) / validateFold) + 1; 
            q2 = floor(length(fsTrnX{m}) * fold     / validateFold) + 1;
            foldValX{m} = fsTrnX{m}(q1 : q2-1);
            foldTrnX{m} = fsTrnX{m}([1:q1-1, q2:length(fsTrnX{m})]);
        end
        % discard logLikArr and logCondLiksArr
        cvMghmm = mmieGenTrainNew(foldTrnX, nState, varargin{:});
        logOdds{1} = tramPredict(foldValX{1},cvMghmm,varargin{:});
        logOdds{2} = tramPredict(foldValX{2},cvMghmm,varargin{:});
        fsAcc = fsAcc + (sum(logOdds{1} <= 0) + sum(logOdds{2} > 0)) ...
                        / (length(foldValX{1}) + length(foldValX{2}));
    end
    msAcc1 = fsAcc / validateFold;

% parallelized version of calModelSelAcc()
function msAcc1 = parCalModelSelAcc(fsTrnX, validateFold, nState, varargin)
    fsAcc = zeros(1,validateFold);
    parfor (fold = 1:validateFold, validateFold)
        logOdds = cell(1, 2);
        foldTrnX = cell(2, 1); 
        foldValX = cell(2, 1);
        for m = 1 : 2
            q1 = floor(length(fsTrnX{m}) * (fold-1) / validateFold) + 1; 
            q2 = floor(length(fsTrnX{m}) * fold     / validateFold) + 1;
            foldValX{m} = fsTrnX{m}(q1 : q2-1);
            foldTrnX{m} = fsTrnX{m}([1:q1-1, q2:length(fsTrnX{m})]);
        end
        % discard logLikArr and logCondLiksArr
        cvMghmm = mmieGenTrainNew(foldTrnX, nState, varargin{:});
        logOdds{1} = tramPredict(foldValX{1},cvMghmm,varargin{:});
        logOdds{2} = tramPredict(foldValX{2},cvMghmm,varargin{:});
        fsAcc(fold) = (sum(logOdds{1} <= 0) + sum(logOdds{2} > 0)) ...
                        / (length(foldValX{1}) + length(foldValX{2}));
    end
    msAcc1 = mean(fsAcc);

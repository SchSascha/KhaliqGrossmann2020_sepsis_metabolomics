function [mghmms,logLikArr,logCondLiksArr,trainAcc,testAcc]=mmieGenTrainNew(...
         trainX,nState,varargin)

replicateGiven = false;
for k = length(varargin)
    if ischar(varargin{k}) && strcmpi('replicates', varargin{k})
        replicateGiven = true;
    end
end
[nRepeatMmi, algArg] = varArgRemove('repeatmmi', 0, varargin);
[learnRtType, algArg] = varArgRemove('mmielearnrate', 0, algArg);
learnRtTypeTr = varArgRemove('mmietrlearnrate', 0, algArg);

% If there is no discriminative training
if isnumeric(learnRtType) && isnumeric(learnRtTypeTr) ...
                          && learnRtType==0 && learnRtTypeTr==0
    if nRepeatMmi~=0,error('RepeatMmi cannot be given without MMIE training');end
    [mghmms, logLikArr] = patientGhmmTrain(trainX, 'states', nState, algArg{:});
    logCondLiksArr = {}; trainAcc = []; testAcc = []; return;
elseif nRepeatMmi == 0 % random init chosen by generative training step
    [mghmms, logLikArr] = patientGhmmTrain(trainX, 'states', nState, algArg{:});
    if nargout >= 5
        [mghmms, logCondLiksArr, trainAcc, testAcc] = mmieGenTrain(trainX, ...
                                                      mghmms, algArg{:});
    else
        [mghmms, logCondLiksArr, trainAcc] = mmieGenTrain(trainX, ...
                                                          mghmms, algArg{:});
    end
    return;
end

if replicateGiven, error('RepeatMmi and replicates cannot both be given'); end
nCate = length(trainX);
logLikArr = cell(nCate, 1);
repMghmm = cell(nRepeatMmi, 1);
for iRepeat = 1 : nRepeatMmi, repMghmm{iRepeat} = cell(nCate, 1); end
nPatients = 0;
for c = 1 : nCate,
    [tr1 mu1 sigma1 logLikArr{c} mghmmArr] ...
        = mghmmTrainNew(trainX{c}, 'replicates', nRepeatMmi, varargin{:});
    for iRepeat = 1 : nRepeatMmi, repMghmm{iRepeat}{c} = mghmmArr{iRepeat}; end
    nPatients = nPatients + length(trainX{c});
end
for iRepeat = 1 : nRepeatMmi
    for c = 1 : nCate, 
        repMghmm{iRepeat}{c}.prior = length(trainX{c}) / nPatients; 
    end
end

topLogCondLik = -Inf;
logCondLiksArr = cell(nRepeatMmi, 1);
trainAcc = cell(nRepeatMmi, 1); testAcc = cell(nRepeatMmi, 1);
for iRepeat = 1 : nRepeatMmi
    [mghmms1, logCondLik1, trainAcc1, testAcc1] = mmieGenTrain(trainX, ...
        repMghmm{iRepeat}, varargin{:});
    logCondLiksArr{iRepeat} = logCondLik1;
    trainAcc{iRepeat} = trainAcc1; testAcc{iRepeat} = testAcc1;
    if max(logCondLik1) > topLogCondLik
        topLogCondLik = max(logCondLik1); mghmms = mghmms1;
    end
end

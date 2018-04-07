function [CrossValidationModel] = crossvalidation_sgtmtt(SDL, v, opt)
%CROSSVAL K-fold cross validation.
%   [vR, vC] = CROSSVAL(SDL, v, opt) returns the leave-one-out
%   cross-validation classification estimates of the observations
%   described by the feature matrix given in SDL.data using the labels in SDL.labels
%   SDL.data must be a N x T x D matrix, with N as the number of samples, 
%   T the number of time points and D the number of feature dimensions
%   SDL.labels must be given as N x T 
%   e.g. from a long list oSDL.labels = reshape(Labels,T,N)'
%   It is assumed that the time sequences have equal length and are aligned
%   v is the number of folds and opt contains (optional) arguments passed to sgtmtt.
%   with svm classifier
%   call: 
%   CrossValidationModel = crossvalidation_sgtmtt(oSDL, 5, {parameter_for_sgtmtt_comma_separated} )
%
X = SDL.data;
L = SDL.labels;

error(nargchk(3, 4, nargin))
[N,T,D] = size(X);
if nargin < 2
  v = [];
end

if nargin < 3
  opt = [];
end

if isempty(v)
  V = 1:N;
  v = N;
elseif length(v) == 1
  if ~isa(v, 'double') | ~isreal(v) | v <= 0 | round(v) ~= v
    error('V must be a positive, non-zero integer.')
  elseif v > N
    error('V must be less than the number of observations.')
  elseif v == N
    V = 1:N;
  else
    w = v*fix(N/v);
    r = randperm(v);
    V(randperm(N)) = [reshape(repmat(1:v, w/v, 1), w, 1); ...
		      r(1:N-w)'];
  end
else
  if ~isa(v, 'double') | ~isreal(v) | prod(size(v)) ~= length(v) | ...
	any(round(v) ~= v | v <= 0 | isinf(v))
    error(['Cross-validation set indeces V must be a vector of' ...
	   ' positive, finite, non-zero integers.'])
  elseif length(v) ~= N
    error(['Length of set indeces V must be same as number of' ...
	   ' observations in X and K.'])
  elseif ~all(sum(sparse(1:N, v, 1)))
    error('Set indeces V may not have empty groups.')
  end
  
  V = v;
  v = max(V);
end
  
c = zeros(N, 1);
vModelResults    = [];
vClassifyResults = [];
vClassifyResultsProb  = [];
vClassifyResultsProbWithReject = [];
mConfusionMatrix = [];
vUniqueLabels = unique(SDL.labels);
for i = 1:v
  s = warning;
  warning off
  SDLCurrent = SDL;
  SDLCurrent.data   = SDLCurrent.data(V ~= i, :, :);
  SDLCurrent.labels = SDLCurrent.labels(V ~= i);
  vCurrentLabelsReshaped = reshape(permute(SDLCurrent.labels,[2 1]),prod(size(SDLCurrent.labels)),1);
  if(length(opt)<3)
        [vModels,vUniqueLabels] = feval(@GTMTT_Supervised, SDLCurrent.data, vCurrentLabelsReshaped, opt{1}, opt{2});
  else
        [vModels,vUniqueLabels] = feval(@GTMTT_Supervised, SDLCurrent.data, vCurrentLabelsReshaped, opt{1}, opt{2},opt{3:end});
  end
  sModel = struct('Models',{vModels},'Labels',{vUniqueLabels},'SVMModel',[],'SVMModelLabels',[],'V',{V});  
  warning(s);
  
  % get the kernel for - all - data, but limit it later on again
  % pickup the log-likelihoods for the data
  [N,T,D] = size(SDL.data);
  mLiks   = zeros(N, length(vModels));
  Tab     = tabulate(SDL.labels);
  vPrior  = Tab(:,end)/100;
  vFrameIndices = [];
  mStateSequences = zeros(N,T*length(vModels));
  if(length(opt)>2 && ~isscalar(opt{3})) % time frame classification given ?
      vFrameIndices = opt{3};
  end
  for(p=1:length(vModels))
     for(k=1:N), 
        mCurrentSample = reshape(permute(SDL.data(k,:,:),[2 1 3]), T, D);
        if(~isempty(vFrameIndices))
            mCurrentSample = mCurrentSample(vFrameIndices,:);
        end
        [Lik,Resp] = test(vModels{p}.oModel1,mCurrentSample); 
        iStart = T*(p-1)+1;
        mStateSequences(k,[iStart:iStart+T-1]) = argmax(Resp{1});
        mLiks(k,p) = Lik;
     end
  end
  % mLiks are the log-probs of observing the given (sub-)pattern at t in state i running overall t
  % due to log its accumulated \sum_i p_1(i_1) * ... * \sum_i p_T(i_T) -->
  % log(prod(...) --> sum(log(sum_i p_1(t_i)), ..., log(sum_T p_T(i_T))
  % normalize model probabilities to sum to 1
  mProbs = bsxfun(@minus,mLiks,min(mLiks')');
  mProbs = bsxfun(@rdivide,mProbs,sum(mProbs,2));
  I=[1:N]';

  if(isnumeric(SDLCurrent.labels))
    vLabels = SDLCurrent.labels(:,1);
  else      
    vLabels = str2num(cell2mat(SDLCurrent.labels(:,1))); % only first entry
  end  
     if(isnumeric(SDL.labels))
        vLabels = SDL.labels(:,1);
      else      
        vLabels = str2num(cell2mat(SDL.labels(:,1))); % only first entry
      end  

      SDLTest.data   = I(V == i);
      SDLTest.labels = vLabels(V == i); 

      if(isnumeric(SDLTest.labels))
        vLabels = SDLTest.labels;
      else      
        vLabels = str2num(cell2mat(SDLTest.labels));
      end
      vClassifyResults = [vClassifyResults 0];
%  end
  vModelResults        = [vModelResults sModel];  
  'prob - acc'
  vCurrentTest            = find(V==i);
  PredictedPrototypeIndex = argmax(mLiks(V==i,:)');
  MeanLiks                = mean(mLiks);    
  IsSafePrediction        = arrayfun(@(x) (mLiks(x,PredictedPrototypeIndex(x))-(MeanLiks(PredictedPrototypeIndex(x))))>0, [1:numel(PredictedPrototypeIndex)])';  
  dAccProbWithReject      = sum(vUniqueLabels(PredictedPrototypeIndex) == vLabels & IsSafePrediction==1)./length(vLabels);
  dAccProb                = sum(vUniqueLabels(PredictedPrototypeIndex) == vLabels)./length(vLabels);
  IndexMatchingWithRejectClass = [ones(numel(PredictedPrototypeIndex),1)*numel(vUniqueLabels)+1 vUniqueLabels(PredictedPrototypeIndex)];
  ClassificationWithReject = IndexMatchingWithRejectClass(sub2ind(size(IndexMatchingWithRejectClass),[1:numel(PredictedPrototypeIndex)]', (vUniqueLabels(PredictedPrototypeIndex) == vLabels & IsSafePrediction)+1));
  vClassifyResultsProb    = [vClassifyResultsProb dAccProb]; 
  vClassifyResultsProbWithReject = [vClassifyResultsProbWithReject dAccProbWithReject]; 
  mConfusionMatrix        = [mConfusionMatrix; ConfusionMatrix(ClassificationWithReject,vLabels,vUniqueLabels)];   
  vWrongs{i}              = vCurrentTest(vUniqueLabels(argmax(mLiks(vCurrentTest,:)')) ~= vLabels);
end
CrossValidationModel      = struct('Models',vModelResults, 'Acc', vClassifyResultsProb, 'AccWithReject',vClassifyResultsProbWithReject, 'ConfusionMatrices', mConfusionMatrix, 'Wrongs',vWrongs);

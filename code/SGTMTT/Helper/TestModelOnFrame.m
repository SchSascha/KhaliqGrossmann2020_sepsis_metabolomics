% Helper function to test new data against the model 
% Assumes labeling of classes [1-C]
% Sample using a cross-validation run
% vAcc = []; for(l=1:10:800) k=1;vRange=[l:l+9];SDLTrain.data = SDL.data(vModelResults(k).V~=k,:,:); SDLTrain.labels = SDL.labels(vModelResults(k).V~=k);  SDLTest.data = SDL.data(vModelResults(k).V==k,:,:); SDLTest.labels = SDL.labels(vModelResults(k).V==k); [dAcc] = TestModelOnFrame(SDLTrain,SDLTest,vModelResults(k),vRange); vAcc = [vAcc;dAcc];end
% 
%
function [dAcc] = TestModelOnFrame(SDLTrain,SDLTest,oModelResults,vFrameIndices,vUniqueLabels)
    LTest   = SDLTest.labels;
    NTest   = length(LTest);
    LTrain  = SDLTrain.labels;
    NTrain  = length(LTrain);
    T       = numel(vFrameIndices);
    D       = size(SDLTest.data,3);
    DataAll = [SDLTrain.data; SDLTest.data];
    DataAll = DataAll(:,vFrameIndices,:);
    vModels = oModelResults.Models;
    NAll    = NTrain+NTest;
    
    % get the likelihoods from the given models for all data but only on
    % the limited time range
    mLiks   = zeros(NAll,length(vModels));    
    for(p=1:length(vModels))
        for(k=1:NAll),
            mCurrentSample = reshape(permute(DataAll(k,:,:),[2 1 3]), T, D);
            [Lik]          = test(vModels{p}.oModel1,mCurrentSample);
            mLiks(k,p)     = Lik;
        end
    end
    dAcc   = sum(vUniqueLabels(argmax(mLiks(NTrain+1:end,:)'))==LTest)/numel(LTest);

clear 

% for the simulated data
load('lin_sim.txt');

% reshape to 100 (persons) x 8 (time points) x 100 variables
iSamples    = 100;
iDims       = 100;
iTimePoints = 8;
K           = 9;             % number of states - in acc. to maximum setting in Lin but squared
M           = 4;             % number of basis functions 

for(k=1:iSamples),
    for(t=1:iTimePoints),
        X(k,t,:)=reshape(lin_sim(k,[3+(t-1)*iDims:3+(t-1)*iDims+iDims-1]),1,iDims);
    end
end
L = [zeros(50,1);ones(50,1);]+1; % labels
X = X(:,:,[1:30]);

% train the model
[vModels] = GTMTT_Supervised(X,L,K,M);
' model generation finished - press any key '

% plot the learned metric parameter
pause
E=get(vModels{1}.oModel1,'ometric')
bar(E.Par)

% plot the map and the path
for(c=1:2)
    Resp = get(vModels{c}.oModel1,'r');
    mod  = pmd(Resp{1});
    iTDim = size(X,2);
    %... now including paths
    membmap(K,mod,[1:iTDim]);
end

pause;
close all;

% plot the map and the data for the first class in 3 D (dim 3 matters most
% (30 dimensions)
iDims      = 30;
iTimepoint = 7; % most relevant time point
PlotGridAndData2D(vModels{1}.oModel1,6,7,iTimepoint,iTimePoints,iDims,0) % first class
PlotGridAndData2D(vModels{2}.oModel1,6,7,iTimepoint,iTimePoints,iDims,1) % second class data
pause;
figure;
iTimepoint = 3;
PlotGridAndData3D(vModels{1}.oModel1,6,7,9,iTimepoint,iTimePoints,iDims,0) % first class
PlotGridAndData3D(vModels{2}.oModel1,6,7,9,iTimepoint,iTimePoints,iDims,1) % second class data

pause;
% run a full 10-fold cross-validation
SDL.data   = X;
SDL.labels = L;
[oCVModel] = crossvalidation_sgtmtt(SDL,10,{K,M})
return;
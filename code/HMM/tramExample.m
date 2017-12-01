tic
rand('state', 10000);
% load simulated data (already log-transformed)
simData = load('sim.mat');
data = simData.data;
% only use 7 of the 8 time points
for m = 1 : 2
    for k = 1 : length(data{m}), data{m}{k} = data{m}{k}(:, 1:7); end
end
% the first fold of 4-fold cross validation
trainData = {data{1}(13:50), data{2}(13:50)};
testData = {data{1}(1:12), data{2}(1:12)};
% generative training
nState = 2;
nFeatArr = [100 50 40 30:-2:2];
model = tramGenTrain(trainData, nState, nFeatArr, 'replicates',10,'maxiterations',20);
% select genes
fsTrainData = selectFeature(trainData, model{1}.selGenes);
fsTestData  = selectFeature(testData , model{1}.selGenes);
% classify testing data
logOdds = cell(1, 2);
logOdds{1} = tramPredict(fsTestData{1}, model);
logOdds{2} = tramPredict(fsTestData{2}, model);
% calculate and print the accuracy 
accuracy = (sum(logOdds{1} <= 0) + sum(logOdds{2} > 0)) ...
           / (length(testData{1}) + length(testData{2}));
fprintf('Accuracy of generative HMM is %2.0f%%\n',100 * accuracy);

% discriminative training and evaluation
discModel = tramDiscTrain(fsTrainData, model, 'mmieIterations', 1000);
discLogOdds = cell(1, 2);
discLogOdds{1} = tramPredict(fsTestData{1}, discModel);
discLogOdds{2} = tramPredict(fsTestData{2}, discModel);
discAcc = (sum(discLogOdds{1} <= 0) + sum(discLogOdds{2} > 0)) ...
           / (length(testData{1}) + length(testData{2}));
fprintf('Accuracy of discriminative HMM is %2.0f%%\n',100 * discAcc);
toc

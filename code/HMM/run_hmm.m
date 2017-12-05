%% Load and clean the data, then run the loop-HMM from TRAM
% Includes a modification for loop-jump-HMM, where the loop HMM is extended
% by jumps over states. Results are OK so far.

tic;

%% Load human clinical data
filename = '../../data/measurements/Summary human sample data.csv';
[hd, header] = import_human_data(filename);
[hdm, hdp, hdmh, hdph] = split_human_metab_pheno(hd, header);

%% Filter columns with zero SD
remove_idx = zeros(length(hdmh),1);
for n = 6:length(hdm(1,:))
    col = hdm(:, n);
    col(cellfun(@isempty, col)) = [];
    numdat = cellfun(@(x) x(1), col);
    remove_idx(n) = std(numdat) == 0;
end
remove_idx = find(remove_idx);
hdm(:, remove_idx) = [];

%% Normalize columns
for n = 6:length(hdm(1,:))
    col = hdm(:, n);
    empty_bool = cellfun(@isempty, col);
    numdat = cellfun(@(x) x(1), col(~empty_bool));
    m = mean(numdat);
    s = std(numdat);
    numdat = (numdat - m) / s;
    hdm(~empty_bool, n) = num2cell(numdat);
end

%% Reorganize data to fit TRAM
% Prepare
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
head_pat_idx = find(cellfun(cellfind('Patient'), hdmh));
head_surv_idx = find(cellfun(cellfind('Survival'), hdmh));
pats = cell2mat(hdm(:, head_pat_idx));
surv = hdm(:, head_surv_idx);
[u_pats, u_idx] = unique(pats);
u_surv = surv(u_idx);
u_pat_excl_idx = diff(u_idx) <= 0;
u_pat_excl_idx(end+1) = (length(pats) - u_idx(end)) <= 0;
u_pat_excl_idx = find(u_pat_excl_idx);
% remove patients with only one time point
if ~isempty(u_pat_excl_idx)
    u_pats(u_pat_excl_idx) = [];
    u_surv(u_pat_excl_idx) = [];
end
% set NaNs
for n = 1:length(pats)
    hdm(n, find(cellfun(@isempty, hdm(n,:)))) = {nan};
end
% Build data structure
data = cell(1, 2);
data{1} = cell(1, sum(strcmp(u_surv, 'S')));
data{2} = cell(1, sum(strcmp(u_surv, 'NS')));
% Fill data structure
d1_idx = 1;
d2_idx = 1;
for n = 1:length(u_pats)
    if strcmp(u_surv(n), 'S')
        data{1}{d1_idx} = cell2mat(hdm(find(pats == u_pats(n)), 6:end))';
        d1_idx = d1_idx + 1;
    else
        data{2}{d2_idx} = cell2mat(hdm(find(pats == u_pats(n)), 6:end))';
        d2_idx = d2_idx + 1;
    end
end
% CV fold generation, randomly permute order of instances
nfolds = 5;
order1 = randperm(length(data{1}));
order2 = randperm(length(data{2}));
CV_set1 = cell(1,nfolds);
CV_set2 = cell(1,nfolds);
breaks1 = [0, floor((1:(nfolds-1)) * length(data{1}) / (nfolds)), length(data{1})];
breaks2 = [0, floor((1:(nfolds-1)) * length(data{2}) / (nfolds)), length(data{2})];
cvdata1 = data{1}(order1);
cvdata2 = data{2}(order2);
% for k = 1:nfolds
%     CV_set1{k} = cvdata1((breaks1(k)+1):breaks1(k + 1));
%     CV_set2{k} = cvdata2((breaks2(k)+1):breaks2(k + 1));
% end
%%% Cheap validation set selection
d1_start = length(data{1}) * 1 / 5;
d2_start = length(data{2}) * 1 / 5;
trainData = {data{1}(ceil(d1_start):length(data{1})), data{2}(ceil(d2_start):length(data{2}))};
testData = {data{1}(1:floor(d1_start)), data{2}(1:floor(d2_start))};

%% Set central parameters
% Number of system states
nState = 3;
% Number of features to try, has to start with all features
nFeatArr = [length(data{1}{1}(:,1)) 100 50 40 30:-4:2];
%nFeatArr = [length(data{1}{1}(:,1)) 100 50 20 10 6 2];

%% Train and test HMM
% generative training
parEval = license('test', 'distrib_computing_toolbox');
[model, int_top_acc] = tramGenTrain(trainData, nState, nFeatArr, 'replicates', 10,'maxiterations', 60, 'model', 'loophmm', 'parEval', parEval);
% select genes
fsTrainData = selectFeature(trainData, model{1}.selGenes);
fsTestData  = selectFeature(testData , model{1}.selGenes);
% classify testing data
logOdds1 = tramPredict(fsTestData{1}, model);
logOdds2 = tramPredict(fsTestData{2}, model);
% calculate and print the accuracy
acc = 100 * (sum(logOdds1 <= 0) + sum(logOdds2 > 0)) ...
    / (length(testData{1}) + length(testData{2}));
tpr = 100 * sum(logOdds1 <= 0)/length(testData{1});
tnr = 100 * sum(logOdds2 > 0)/length(testData{2});
% discriminative training and evaluation
discModel = tramDiscTrain(fsTrainData, model, 'mmieIterations', 5000);
% classify testing data
discLogOdds1 = tramPredict(fsTestData{1}, discModel);
discLogOdds2 = tramPredict(fsTestData{2}, discModel);
dacc = 100 * (sum(discLogOdds1 <= 0) + sum(discLogOdds2 > 0)) ...
    / (length(testData{1}) + length(testData{2}));
dtpr = 100 * sum(discLogOdds1 <= 0)/length(testData{1});
dtnr = 100 * sum(discLogOdds2 > 0)/length(testData{2});

% Report Accuracy
fprintf('Accuracy of generative HMM by internal CV is %2.0f%%\n', 100 * int_top_acc);

fprintf('Accuracy of generative HMM is %2.0f%%\n',acc);
fprintf('TPR of generative HMM is %2.0f%%\n',tpr);
fprintf('TNR of generative HMM is %2.0f%%\n',tnr);
   
fprintf('Accuracy of discriminative HMM is %2.0f%%\n',dacc);
fprintf('TPR of discriminative HMM is %2.0f%%\n',dtpr);
fprintf('TNR of discriminative HMM is %2.0f%%\n',dtnr);

toc;
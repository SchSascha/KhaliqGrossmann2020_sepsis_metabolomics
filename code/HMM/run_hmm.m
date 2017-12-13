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
hdmh(remove_idx) = [];

%% Filter rows with zero SD?
for n = 1:size(hdm,1)
    row = hdm(:, n);
    row(cellfun(@isempty, row)) = [];
    numdat = cellfun(@(x) x(1), row);
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
    %numdat = (numdat - m) / s;
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
pat = cell(1,2);
% Fill data structure
d1_idx = 1;
d2_idx = 1;
for n = 1:length(u_pats)
    if strcmp(u_surv(n), 'S')
        data{1}{d1_idx} = cell2mat(hdm(find(pats == u_pats(n)), 6:end))';
        pat{1}{d1_idx} = u_pats(n);
        d1_idx = d1_idx + 1;
    else
        data{2}{d2_idx} = cell2mat(hdm(find(pats == u_pats(n)), 6:end))';
        pat{2}{d2_idx} = u_pats(n);
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
pat{1} = pat{1}(order1);
pat{2} = pat{2}(order2);
for k = 1:nfolds
    CV_set1{k} = cvdata1((breaks1(k)+1):breaks1(k + 1));
    CV_set2{k} = cvdata2((breaks2(k)+1):breaks2(k + 1));
end

%% Set central parameters
% Number of system states
nState = 3;
% Number of features to try, has to start with all features
nFeatArr = [length(data{1}{1}(:,1)) 100 50 40 30:-2:2];
%nFeatArr = [length(data{1}{1}(:,1)) 100 50 20 10 6 2];
%nFeatArr = [length(data{1}{1}(:,1))];

%% Set performance measures
acc = zeros(1, nfolds);
tpr = zeros(1, nfolds);
tnr = zeros(1, nfolds);
dacc = zeros(1, nfolds);
dtpr = zeros(1, nfolds);
dtnr = zeros(1, nfolds);
sacc = zeros(1, nfolds);
stpr = zeros(1, nfolds);
stnr = zeros(1, nfolds);
int_top_acc = zeros(1, nfolds);
logOdds1 = zeros(1,length(cvdata1));
logOdds2 = zeros(1,length(cvdata2));
discLogOdds1 = zeros(1,length(cvdata1));
discLogOdds2 = zeros(1,length(cvdata2));
singleLogOdds1 = zeros(1,length(cvdata1));
singleLogOdds2 = zeros(1,length(cvdata2));


%% Train and test HMM
for n = 1:nfolds
%n = 2;
    train_cv_sel = 1:nfolds;
    train_cv_sel(n) = [];
    trainData = cell(1,2);
    for s = train_cv_sel
        for l = 1:length(CV_set1{s})
            trainData{1}{length(trainData{1})+1} = CV_set1{s}{l};
        end
        for l = 1:length(CV_set2{s})
            trainData{2}{length(trainData{2})+1} = CV_set2{s}{l};
        end
    end
    testData = {CV_set1{n}, CV_set2{n}};
    % generative training
    parEval = license('test', 'distrib_computing_toolbox');
    %parEval = false;
    [model, int_top_acc(n)] = tramGenTrain(trainData, nState, nFeatArr, 'replicates', 10,'maxiterations', 500, 'model', 'loopjumphmm', 'parEval', parEval);
    % select genes
    fsTrainData = selectFeature(trainData, model{1}.selGenes);
    fsTestData  = selectFeature(testData , model{1}.selGenes);
    % classify testing data
    ar_range1 = (breaks1(n)+1):breaks1(n + 1);
    ar_range2 = (breaks2(n)+1):breaks2(n + 1);
    logOdds1(ar_range1) = tramPredict(fsTestData{1}, model);
    logOdds2(ar_range2) = tramPredict(fsTestData{2}, model);
    % calculate and print the accuracy
    acc(n) = 100 * (sum(logOdds1(ar_range1) <= 0) + sum(logOdds2(ar_range2) > 0)) ...
        / (length(testData{1}) + length(testData{2}));
    tpr(n) = 100 * sum(logOdds1(ar_range1) <= 0)/length(testData{1});
    tnr(n) = 100 * sum(logOdds2(ar_range2) > 0)/length(testData{2});
    % discriminative training and evaluation
    %discModel = tramDiscTrain(fsTrainData, model, 'mmieIterations', 500);
    % classify testing data
    %discLogOdds1(ar_range1) = tramPredict(fsTestData{1}, discModel);
    %discLogOdds2(ar_range2) = tramPredict(fsTestData{2}, discModel);
    %dacc(n) = 100 * (sum(discLogOdds1(ar_range1) <= 0) + sum(discLogOdds2(ar_range2) > 0)) ...
    %    / (length(testData{1}) + length(testData{2}));
    %dtpr(n) = 100 * sum(discLogOdds1(ar_range1) <= 0)/length(testData{1});
    %dtnr(n) = 100 * sum(discLogOdds2(ar_range2) > 0)/length(testData{2});
    % reduce test set to first day
    fsTestDataSingle = cell(1, 2);
    for m = 1 : 2
        for k = 1 : length(fsTestData{m}), fsTestDataSingle{m}{k} = fsTestData{m}{k}(:, 1); end
    end
    % classify test data
    singleLogOdds1(ar_range1) = tramPredict(fsTestDataSingle{1}, model);
    singleLogOdds2(ar_range2) = tramPredict(fsTestDataSingle{2}, model);
    % calculate and print the accuracy
    sacc(n) = 100 * (sum(singleLogOdds1(ar_range1) <= 0) + sum(singleLogOdds2(ar_range2) > 0)) ...
        / (length(testData{1}) + length(testData{2}));
    stpr(n) = 100 * sum(singleLogOdds1(ar_range1) <= 0)/length(testData{1});
    stnr(n) = 100 * sum(singleLogOdds2(ar_range2) > 0)/length(testData{2});
end

% Report accuracy
fprintf('Accuracy of generative HMM by internal CV is %2.0f%%\n', mean(100 * int_top_acc));

fprintf('Accuracy of generative HMM is %2.0f%%\n',mean(acc));
fprintf('TPR of generative HMM is %2.0f%%\n',mean(tpr));
fprintf('TNR of generative HMM is %2.0f%%\n',mean(tnr));
   
fprintf('Accuracy of discriminative HMM is %2.0f%%\n',mean(dacc));
fprintf('TPR of discriminative HMM is %2.0f%%\n',mean(dtpr));
fprintf('TNR of discriminative HMM is %2.0f%%\n',mean(dtnr));

% print performance on single day samples
fprintf('Accuracy of generative HMM on first day only is %2.0f%%\n',mean(sacc));
fprintf('TPR of generative HMM is %2.0f%%\n',mean(stpr));
fprintf('TNR of generative HMM is %2.0f%%\n',mean(stnr));

% Retrain on all samples
trainData = {cvdata1, cvdata2};
parEval = license('test', 'distrib_computing_toolbox');
[model, int_top_acc_all] = tramGenTrain(trainData, nState, nFeatArr, 'replicates', 20,'maxiterations', 500, 'model', 'loopjumphmm', 'parEval', parEval);

% Print performance measures
fprintf('Accuracy of generative HMM by internal CV an all samples is %2.0f%%\n', 100 * int_top_acc_all);

toc;

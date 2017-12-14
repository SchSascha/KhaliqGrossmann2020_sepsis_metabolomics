%% Decode TRAM HMM and calculate Viterbi path

% Decode survivor time courses
subset = 1;
decoded_S_seqs = cell(1, length(data{subset}));
for set = 1:length(data{subset})
    seq = cvdata1{set}(model{subset}.selGenes,:);
    tr = model{subset}.tr;
    mu = model{subset}.mu;
    sigma = model{subset}.sigma;
    [pStates, loglik, fs, bs, s] = mghmmDecode(seq, tr, mu, sigma);
    decoded_S_seqs{set} = {pStates, loglik, fs, bs, s};
end

% Decode nonsurvivor time courses
subset = 2;
decoded_NS_seqs = cell(1, length(data{subset}));
for set = 1:length(data{subset})
    seq = cvdata2{set}(model{subset}.selGenes,:);
    tr = model{subset}.tr;
    mu = model{subset}.mu;
    sigma = model{subset}.sigma;
    [pStates, loglik, fs, bs, s] = mghmmDecode(seq, tr, mu, sigma);
    decoded_NS_seqs{set} = {pStates, loglik, fs, bs, s};
end

% Calculate survivor Viterbi paths
subset = 1;
viterbi_S_seqs = cell(1, length(data{subset}));
for set = 1:length(data{subset})
    seq = cvdata1{set}(model{subset}.selGenes,:);
    tr = model{subset}.tr;
    mu = model{subset}.mu;
    sigma = model{subset}.sigma;
    [estStates, loglikSeq] = mghmmViterbi(seq, tr, mu, sigma);
    viterbi_S_seqs{set} = {estStates, loglikSeq};
end

% Calculate nonsurvivor Viterbi paths
subset = 2;
viterbi_NS_seqs = cell(1, length(data{subset}));
for set = 1:length(data{subset})
    seq = cvdata2{set}(model{subset}.selGenes,:);
    tr = model{subset}.tr;
    mu = model{subset}.mu;
    sigma = model{subset}.sigma;
    [estStates, loglikSeq] = mghmmViterbi(seq, tr, mu, sigma);
    viterbi_NS_seqs{set} = {estStates, loglikSeq};
end

%% Calculate state as if all samples are from single patients with unknown label
% Calculate survivor Viterbi paths
viterbi_blind_seqs = cell(2, size(hdm, 1));
row = 1;
for part = 1:2
for patient = 1:length(data{part})
for sample = 1:size(data{part}{patient}, 2)
    seq = data{part}{patient}(model{1}.selGenes,sample);
    tr = model{1}.tr;
    mu = model{1}.mu;
    sigma = model{1}.sigma;
    [estStates, loglikSeq] = mghmmViterbi(seq, tr, mu, sigma);
    viterbi_blind_seqs{1, row} = estStates;
    tr = model{2}.tr;
    mu = model{2}.mu;
    sigma = model{2}.sigma;
    [estStates, loglikSeq] = mghmmViterbi(seq, tr, mu, sigma);
    viterbi_blind_seqs{2, row} = estStates;
    row = row + 1;
end
end
end

%% Prepare long format table for reading and plotting in R
% Prepare file and format
fileID = fopen('tram_viterbi_paths.csv', 'w');
% Write table header
%viterbi_header = ['Status', 'State', 'Case', hdmh{model{1}.selGenes + 5}]';
viterbi_header = {'Status', 'State', 'Case'};
formatspec = [repmat('%s\t', 1, length(viterbi_header) - 1), '%s', '\n'];
fprintf(fileID, formatspec, viterbi_header{1}, viterbi_header{2}, viterbi_header{3});
% Write table
for n = 1:length(viterbi_S_seqs)
    for l = 1:length(viterbi_S_seqs{n}{1})
        fprintf(fileID, formatspec, 'S', num2str(viterbi_S_seqs{n}{1}(l)), num2str(pat{1}{n}));%, cvdata1{n}(model{1}.selGenes, l)');
    end
end
for n = 1:length(viterbi_NS_seqs)
    for l = 1:length(viterbi_NS_seqs{n}{1})
        fprintf(fileID, formatspec, 'NS', num2str(viterbi_NS_seqs{n}{1}(l)), num2str(pat{2}{n}));%, cvdata2{n}(model{2}.selGenes, l)');
    end
end
fclose(fileID);

% Prepare file and format
fileID = fopen('tram_viterbi_blind_paths.csv', 'w');
% Write table header
viterbi_header = {'Status', 'State_S', 'State_NS', 'SampleID'};
formatspec = [repmat('%s\t', 1, length(viterbi_header) - 1), '%s', '\n'];
fprintf(fileID, formatspec, viterbi_header{1}, viterbi_header{2}, viterbi_header{3}, viterbi_header{4});
% Write table
for n = 1:size(viterbi_blind_seqs, 2)
        fprintf(fileID, formatspec, hdm{n, 4}, num2str(viterbi_blind_seqs{1, n}), num2str(viterbi_blind_seqs{2, n}), num2str(hdm{n, 1}));
end
fclose(fileID);

%% Write metabolite names selected by HMM training
fileID = fopen('tram_selected_features.csv', 'w');
% Write table header
formatspec = '%s\n';
% Write table
for l = 1:length(model{1}.selGenes)
    fprintf(fileID, formatspec, hdmh{model{1}.selGenes(l) + 5}{1});
end
fclose(fileID);

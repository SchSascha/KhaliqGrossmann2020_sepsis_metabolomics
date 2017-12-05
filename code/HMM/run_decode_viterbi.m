%% Decode TRAM HMM and calculate Viterbi path

% Decode survivor time courses
subset = 1;
decoded_S_seqs = cell(1, length(data{subset}));
for set = 1:length(data{subset})
    seq = data{subset}{set}(model{subset}.selGenes,:);
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
    seq = data{subset}{set}(model{subset}.selGenes,:);
    tr = model{subset}.tr;
    mu = model{subset}.mu;
    sigma = model{subset}.sigma;
    [pStates, loglik, fs, bs, s] = mghmmDecode(seq, tr, mu, sigma);
    decoded_NS_seqs{set} = {pStates, loglik, fs, bs, s};
end

% Calculate survivor Viterbi paths
subset = 1;
viterbi_S_seqs = cell(1, length(data{subset}));
for set = 1:length(data{subset})Combine Multiple Plots

You can combine plots in several ways. Combine plots in the same axes, or create multiple axes in a figure.
Combine Plots in Same Axes

Use the hold on command to combine multiple plots in the same axes. For example, plot two lines and a scatter plot. Then reset the hold state to off.
    seq = data{subset}{set}(model{subset}.selGenes,:);
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
    seq = data{subset}{set}(model{subset}.selGenes,:);
    tr = model{subset}.tr;
    mu = model{subset}.mu;
    sigma = model{subset}.sigma;
    [estStates, loglikSeq] = mghmmViterbi(seq, tr, mu, sigma);
    viterbi_NS_seqs{set} = {estStates, loglikSeq};
end

%% Prepare long format table for reading and plotting in R
% Prepare file and format
fileID = fopen('tram_viterbi_paths.csv', 'w');
% Write table header
viterbi_header = ['Status', 'State', 'Case', hdmh{model{1}.selGenes + 5}]';
formatspec = [repmat('%s\t', 1, length(viterbi_header) - 1), '%s', '\n'];
fprintf(fileID, formatspec, viterbi_header{:});
% Write table
for n = 1:length(viterbi_S_seqs)
    for l = 1:length(viterbi_S_seqs{n}{1})
        fprintf(fileID, formatspec, 'S', num2str(viterbi_S_seqs{n}{1}(l)), num2str(n), data{1}{n}(model{1}.selGenes, l)');
    end
end
for n = 1:length(viterbi_NS_seqs)
    for l = 1:length(viterbi_NS_seqs{n}{1})
        fprintf(fileID, formatspec, 'NS', num2str(viterbi_NS_seqs{n}{1}(l)), num2str(n + length(viterbi_S_seqs)), data{2}{n}(model{2}.selGenes, l)');
    end
end
fclose(fileID);
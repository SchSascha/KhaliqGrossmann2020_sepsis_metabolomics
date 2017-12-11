%% Write loglikelyhoods to file for ROC plot

% Write generative HMM logliks
formatspec = '%s\t%s\n';
fileID = fopen('gen_logliks.csv', 'w');
fprintf(fileID, formatspec, repmat('S', length(logOdds1), 1), logOdds1);
fprintf(fileID, formatspec, repmat('NS', length(logOdds2), 1), logOdds2);
fclose(fileID);

%Write discriminative HMM logliks
formatspec = '%s\t%s\n';
fileID = fopen('disc_logliks.csv', 'w');
fprintf(fileID, formatspec, repmat('S', length(discLogOdds1), 1), discLogOdds1);
fprintf(fileID, formatspec, repmat('NS', length(discLogOdds2), 1), discLogOdds2);
fclose(fileID);

%Write generative HMM logliks of first day samples
formatspec = '%s\t%s\n';
fileID = fopen('gen_day0_logliks.csv', 'w');
fprintf(fileID, formatspec, repmat('S', length(singleLogOdds1), 1), singleLogOdds1);
fprintf(fileID, formatspec, repmat('NS', length(singleLogOdds2), 1), singleLogOdds2);
fclose(fileID);
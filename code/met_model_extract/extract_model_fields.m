% Saves parts of the metabolic models model to flat files.
%
% author: Peter Grossmann
% date: 06.09.2017
% update: 20.12.2017 (processing of iHsa)
% update: 16.03.2018 (processing of iRno and Recon3D)

%% Initialize COBRA toolbox
code_dir = cd('../cobratoolbox');
initCobraToolbox;
cd(code_dir);

%% Set file path
out_dir = '../../data/template_models/';

%% Load base human metabolism model
% iHSA
%iHsa = xls2model('../../data/iHsa/ncomms14250-s10.xlsx'); % ditched in
%favor of convert-and-forget
%iHsa.rev = iHsa.lb < 0 & iHsa.ub > 0;
% Write iHsa as mat file
%save('../../data/iHsa/iHsa.mat', 'iHsa');
load('../../data/template_models/iHsa.mat');
% iRno
load('../../data/template_models/iRno.mat');
% Recon3D 3.01
load('../../data/template_models/Recon3D_301.mat');
Recon3D.metChEBIID = Recon3D.metCHEBIID;

models = {iHsa, iRno, Recon3D};
mnames = {'iHsa', 'iRno', 'Recon3D'};


%% Write gene list to disk

for m = 1:length(models)
    fileID = fopen(strcat(out_dir, mnames{m}, '_gene_ids.txt'), 'w');
    for n = 1:length(models{m}.genes)
        fprintf(fileID, '%s\n', models{m}.genes{n});
    end
    fclose(fileID);
end

%% Write metabolite names to disk

for m = 1:length(models)
    fileID = fopen(strcat(out_dir, mnames{m}, '_mets.txt'), 'w');
    formatstr = '%s\t%s\t%s\t%s\n';
    fprintf(fileID, formatstr, 'Name', 'ChEBI ID', 'KEGG ID', 'PUBCHEM ID');
    for n = 1:length(models{m}.metNames)
        fprintf(fileID, formatstr, models{m}.metNames{n}, models{m}.metChEBIID{n}, models{m}.metKEGGID{n}, models{m}.metPubChemID{n});
    end
    fclose(fileID);
end

%% Write HMDB IDs to disk
% Only Recon3D has HMDB IDs

fileID = fopen(strcat(out_dir, 'Recon3D_HMDB_IDs.txt'), 'w');
fprintf(fileID, '%s\n', 'HMDB ID');
for n = 1:length(Recon3D.metHMDBID)
    fprintf(fileID, '%s\n', Recon3D.metHMDBID{n});
end
fclose(fileID);

%% Write stoichiometric matrix to disk

for m = 1:length(models)
    S = full(models{m}.S);
    fileID = fopen(strcat(out_dir, mnames{m}, '_stoich_mat.txt'), 'w');
    formatspec = [repmat('%.0f\t', 1, size(S, 2) - 1), '%.0f\n'];
    for n = 1:size(S,1)
        fprintf(fileID, formatspec, S(n ,:));
    end
    fclose(fileID);
end

%% Write reaction rules and properties to disk

for m = 1:length(models)
    fileID = fopen(strcat(out_dir, mnames{m}, '_rxn_properties.txt'), 'w');
    formatstr = '%s\t%s\t%s\n';
    fprintf(fileID, formatstr, 'Gene reaction rule', 'EC', 'Subsystem');
    for n = 1:length(models{m}.grRules)
        fprintf(fileID, formatstr, models{m}.grRules{n}, models{m}.rxnECNumbers{n}, models{m}.subSystems{n}{1});
    end
    fclose(fileID);
end
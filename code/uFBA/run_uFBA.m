%% Read metabolite data and run uFBA
%
% author: Peter Grossmann
% date: 14.06.2018 (initial version)

out_dir = '../../results/uFBA/';

% Prepare out dir
if (exist(out_dir, 'dir') == 0)
    mkdir(out_dir);
end

% Init COBRA toolbox
cobra_dir = '../cobratoolbox/';
code_dir = cd(cobra_dir);
initCobraToolbox;
changeCobraSolver('gurobi', 'ALL');
cd(code_dir);

% Read metabolic template model
load('../../data/template_models/iRno.mat');
load('../../data/template_models/iHsa.mat');
load('../../data/template_models/Recon3D_301.mat');
iRno.mets = strrep(iRno.mets, '[s]', '[e]'); % in iRno [s] is for extracellular compartment
iHsa.mets = strrep(iRno.mets, '[s]', '[e]'); % in iHsa [s] is for extracellular compartment
iRno.comps{strcmp(iRno.comps, 's')} = 'e';
iHsa.comps{strcmp(iHsa.comps, 's')} = 'e';
iRno.compNames{strcmp(iRno.compNames, 's')} = 'e';
iHsa.compNames{strcmp(iHsa.compNames, 's')} = 'e';
models = {iRno, iHsa, Recon3D};

% Collect input files to run uFBA on
dir_cont = dir('../../data/uFBA/*.csv');
in_files = cell(1, length(dir_cont));
in_files = arrayfun(@(x) strcat(x.folder, '/', x.name), dir_cont, 'UniformOutput', false);

mo = arrayfun(@(x) strsplit(x.name, '_'), dir_cont, 'UniformOutput', false);
mo = cellfun(@(x) x{1}, mo, 'UniformOutput', false);
mo = cellfun(@(x) find(ismember({'iRno', 'iHsa', 'Recon3D'}, x)), mo, 'UniformOutput', false);

uFBAmodels = cell(length(in_files), 1);

%neededSinks = {'ascb_L[c]', 'urate[e]', 'gthrd[e]'}; from uFBA example
neededSinks = {'10fthf5glu[c]'};

% Run uFBA on input files
for n = 1%:length(in_files)
    data = read_input_file_for_uFBA(in_files{n});
    mo_struct = models{mo{n}};
    vars.metNames = mo_struct.mets; % order of metabolite critical, has to match order in model
    vars.changeSlopes = data.slope;
    vars.changeIntervals = data.confint;
    vars.neededSinks = neededSinks;
    vars.ignoreSlopes = data.ignore;
    vars.ignoreSlopes(isnan(vars.ignoreSlopes)) = true;
    vars.solvingStrategy = 'case2';
    um = buildUFBAmodel(models{mo{n}}, vars);
    uFBAmodels{n} = um; % The target fields are f (obj) and v (rxn fluxes)
end
uFBAresults = cellfun(@(x) optimizeCbModel(x), uFBAmodels);

% Save built models
save(strcat(out_dir, 'uFBA_models.mat'), 'uFBAmodels', 'uFBAresults');

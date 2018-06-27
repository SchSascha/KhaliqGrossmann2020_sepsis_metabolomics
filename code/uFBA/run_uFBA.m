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
models = {iRno, iHsa, Recon3D};

% Collect input files to run uFBA on
in_files = ls('../../data/uFBA/*.csv');
in_files = strsplit(in_files, '\n|\t');

mo = cellfun(@(x) strsplit(x, '_'), in_files, 'UniformOutput', false);
mo = cellfun(@(x) find(x == ['iRno', 'iHsa', 'Recon3D']), mo);

uFBAmodels = cell(length(in_files));

% Run uFBA on input files
for n = 1:length(in_files)
    data = read_input_file_for_uFBA(in_files{n});
    vars.metNames = data.met;
    vars.changeSlopes = data.slope;
    vars.changeIntervals = data.confint;
    vars.ignoreSlope = data.ignore;
    um = buildUFBAmodel(models{mo{n}}, vars);
    uFBAmodels{n} = um; % The target fields are f (obj) and v (rxn fluxes)
end

% Save built models
save(strcat(out_dir, 'uFBA_models.mat'), 'uFBAmodels');

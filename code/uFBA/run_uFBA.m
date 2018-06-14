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
model = iRno;

% Collect input files to run uFBA on
in_files = ls('../../data/uFBA/');
in_files = strsplit(in_files, '\n|\w');

uFBAmodels = cell(length(in_files));

% Run uFBA on input files
for n = 1:length(in_files) 
    data = read_input_file_for_uFBA(in_files{n});
    vars.metNames = data.met;
    vars.changeSlopes = data.slope;
    vars.changeIntervals = data.confint;
    um = buildUFBAmodel(model, vars);
    uFBAmodels{n} = um; % The target fields are f (obj) and v (rxn fluxes)
end

% Save built models
save(strcat(out_dir, 'models.mat'), 'uFBAmodels');

%% As example:

cobra_dir = '../cobratoolbox/';
code_dir = cd(cobra_dir);
initCobraToolbox;
changeCobraSolver('gurobi', 'ALL');
cd(code_dir);
tutorialPath = fileparts(which('tutorial_uFBA.mlx'));
load([tutorialPath filesep 'sample_data.mat']);
model = readCbModel([tutorialPath filesep 'sample_data.mat'],'modelName','model');
changeSlopes = zeros(length(met_IDs), 1);
changeIntervals = zeros(length(met_IDs), 1);
for i = 1:length(met_IDs)
    [tmp1, tmp2] = regress(met_data(:, i), [time ones(length(time), 1)], 0.05);
    changeSlopes(i, 1) = tmp1(1);
    changeIntervals(i, 1) = abs(changeSlopes(i, 1) - tmp2(1));
end
tmp1 = changeSlopes - changeIntervals;
tmp2 = changeSlopes + changeIntervals;
ignoreSlopes = double(tmp1 < 0 & tmp2 > 0);
uFBAvariables.metNames = met_IDs;
uFBAvariables.changeSlopes = changeSlopes;
uFBAvariables.changeIntervals = changeIntervals;
uFBAvariables.ignoreSlopes = ignoreSlopes;
uFBAoutput = buildUFBAmodel(model, uFBAvariables);
model_ufba = optimizeCbModel(uFBAoutput.model);
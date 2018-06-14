%% Read metabolite data and run uFBA
%
% author: Peter Grossmann
% date: 14.06.2018 (initial version)

% Init COBRA toolbox
cobra_dir = '../cobratoolbox/';
code_dir = cd(cobra_dir);
initCobraToolbox;
changeCobraSolver('gurobi', 'ALL');
cd(code_dir);

% Read metabolic template model
iRno = load('../../data/template_models/iRno.mat');
model = iRno.model;

% Collect input files to run uFBA on
in_files = ls('../../data/uFBA/');
in_files = strsplit(in_files, '\n|\w');

uFBAmodels = cell();

% Run uFBA on input files
for n = 1:length(in_files) 
    data = read_input_file_for_uFBA(in_files{n});
    vars.metNames = data.met;
    vars.changeSlopes = data.slope;
    vars.changeIntervals = data.confint;
    um = buildUFBAmodel(model, vars);
    uFBAmodels{n} = um;
end

% Save built models
save('../../results/uFBA/models.mat', 'uFBAmodels');

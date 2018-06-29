out_dir = '../../results/uFBA/';

% Init COBRA toolbox
cobra_dir = '../cobratoolbox/';
code_dir = cd(cobra_dir);
initCobraToolbox;
changeCobraSolver('gurobi', 'ALL');
cd(code_dir);


tutorialPath = fileparts(which('tutorial_uFBA.mlx'));

load([tutorialPath filesep 'sample_data.mat']);

% We load the model by readCbModel to make sure it fits to the specifications.

model = readCbModel([tutorialPath filesep 'sample_data.mat'],'modelName','model') 

changeSlopes = zeros(length(met_IDs), 1);

changeIntervals = zeros(length(met_IDs), 1);

for i = 1:length(met_IDs)

[tmp1, tmp2] = regress(met_data(:, i), [time ones(length(time), 1)], 0.05);

changeSlopes(i, 1) = tmp1(1);

changeIntervals(i, 1) = abs(changeSlopes(i, 1) - tmp2(1));

end

uFBAvariables.metNames = met_IDs;

uFBAvariables.changeSlopes = changeSlopes;

uFBAvariables.changeIntervals = changeIntervals;

uFBAvariables.ignoreSlopes = ignoreSlopes;

uFBAoutput = buildUFBAmodel(model, uFBAvariables);

model_ufba = optimizeCbModel(uFBAoutput.model)

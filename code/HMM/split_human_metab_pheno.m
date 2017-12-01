function [data_metab, data_pheno, header_metab, header_pheno] = split_human_metab_pheno(data, header)
%split_human_metab_pheno splits the data set into metabolic and
%phenomenological variables
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
pheno_start = find(cellfun(cellfind('Urea'), header));
header_metab = header(1:(pheno_start-1));
data_metab = data(:, 1:(pheno_start-1));
header_pheno = header([1:5, pheno_start:end]);
data_pheno = data(:, [1:5, pheno_start:end]);
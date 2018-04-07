function g = init(g)
% Initialize the parameters Y, W, alpha and beta of the vgtm model
%
% If g.initMeth == 'rand' then W is initialized randomly.  beta is
% initialized through its updating equation using an equal responsability
% matrix.
%
% Else if g.initMeth == 'pca' then W is initialized based on PCA method.
%
% Ivan Olier, 2007      last modification: apr 3/2007

% Compute U
bManualGrid = true;
if isempty(g.U)
    SqrK = sqrt(g.K);
    if fix(SqrK) ~= SqrK
        error('The number of latent points must be squared or, alternatively, the latent space must be manually set.');
    end
    g.U = rctg(g.K);
    bManualGrid = false;
end

if isempty(g.s)
    g.s = 1;
    if g.Verb
        disp('MSG: WidthBF has been automaticaly set to 1.');
    end
end

if isempty(g.M)
    g.M = (round(sqrt(g.K)/2))^2;
    if(g.M<4)
        g.M=4;
    end
    if g.Verb
        disp(['MSG: The number of RBFs has been automaticaly set to ' int2str(g.M)]);
    end
end

% Compute FI
if isempty(g.FI)
    SqrM = sqrt(g.M);
    if fix(SqrM) ~= SqrM
        error('The number of RBFs must be squared or, alternatively, the RBF matrix must be manually set.');
    end
    if(~bManualGrid)
        MU = rctg(g.M); % define a rectangular grid
        MU = MU*(g.M/(g.M-1));
        sigma = g.s*(MU(1,2)-MU(2,2));
    else
        if(isempty(g.MU))
            error('You must specify a grid (MU) for the basis functions in accordance to the latent points grid');
        end
        MU=g.MU;
        MU = MU*(g.M/(g.M-1));
        sigma = g.s * mean(MU(:));
    end    
    g.FI = gbf(MU, sigma, g.U);
end

% Compute XMatFormat
if isempty(g.Data.X)
    error('There is no defined a training dataset. Please, use the ADDDATA function to create it');
end

XMatFormat    = cell2mat(g.Data.X');
if (~isempty(g.Par.W) || ~isempty(g.Par.A) || ~isempty(g.Par.PI)) && g.Verb
    disp('WARNING: Some SGTMTT parameters are not empty');
end

% initialize metric object - just calculate some dummy distances
g.Par.oMetric.distance(XMatFormat(1,:),XMatFormat);

if isempty(g.InitMeth)
    g.InitMeth = 'pca';
end
if strcmp(g.InitMeth,'rnd') 
    [g.Par.W g.Par.Beta] = rndi(g.K,XMatFormat,g.FI,g.Par.oMetric);
elseif strcmp(g.InitMeth,'pca')
    [g.Par.W g.Par.Beta] = pci(XMatFormat, g.U, g.FI,g.Par.oMetric);
else
    error([g.InitMeth 'is not a valid initialization method']);
end

% fix beta g.Par.Beta = 2;

if g.Verb; disp(['MSG: W and BETA initialized using "' g.InitMeth '" method']); end

g.Par.A = trprobi(g.U);
if g.Verb; disp('MSG: Transition Matrix initialized using Gaussians'); end

% make it a loop forward HMM
if(0)
    Z = ones(size(g.Par.A));
    S = spdiags(Z,[0 1],size(Z,1),size(Z,2));
    g.Par.A = g.Par.A .* (full(S));
    g.Par.A = bsxfun(@times,g.Par.A,1.0./sum(g.Par.A,2));
end

g.Par.PI = 1/g.K*ones(g.K,1);
if g.Verb; disp('MSG: Initial state probabilities set to 1/K');end

% Squared distances between X and Y
g.Par.Y=g.FI*g.Par.W;
g.DistXY = g.Par.oMetric.distance(g.Par.Y,XMatFormat);

if isempty(g.Tolerance)
    g.Tolerance = 1e-5;
    if g.Verb
        disp('MSG: The termination tolerance has been set to 1e-5');
    end
end
g.MaxIter = 200;
if isempty(g.MaxIter) && g.Verb
    disp('WARNING: MaxIter is not defined, the training will be stopped only if the convergence is reached.');
end

g.Status='initialized';

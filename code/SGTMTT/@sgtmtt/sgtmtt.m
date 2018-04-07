%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructor of the class sgtm-tt
%   MODEL = SGTMTT(K) build a class for the Standard GTM-TT model with K 
%   (number of latent points)
%   labels as integer vector 1 - first class, 2 - second class a.s.o. hidden states.
%
% (C) Ivan Olier, 2006-2009, GTM-TT code
%     Frank-M. Schleif 2011-2012, SGTM-TT supervised + relevance
% see Demo.m and GTMTT_Supervised.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model1=sgtmtt(K, data, labels,oNorm)
error(nargchk(0, 4, nargin, 'struct'));
if(~exist('oNorm','var'))
    oNorm = [];
end
% Dataset structure
if nargin==0
    Data.X = [];
    Data.D = [];
    Data.N = [];
    labels = [];
    model1.K = [];
elseif nargin>=1
    sqrtK = sqrt(K);
    if sqrtK~=floor(sqrtK)
        error('The number of latent points must be a squared number.');
    else
        model1.K = K;
    end
    if nargin==1
        Data.X = [];
        Data.D = [];
        Data.N = [];
    elseif nargin==3 || nargin==4
        if isstruct(data)
            if ~isfield(data,'X') || ~isfield(data,'D') || ~isfield(data,'N')
                error('Bad format of the input dataset structure');
            end
            Data = data;
        else
            % transpose to full data matrix
            XMat = reshape(permute(data,[2 1 3]),size(data,1) * size(data,2),size(data,3)); 
            % normalize full matrix and store normalization parameters
            if(~isempty(oNorm))
                par1 = oNorm{1};
                par2 = oNorm{2};
                par3 = oNorm{3};
                iMethod = oNorm{4};
            else
                [~, par1, par2, par3, iMethod] = normdata(XMat);            
            end
            %par1 = 0;
            %par2 = 1;
            Data.X = [];
            Data.N = [];
            Data.D = size(XMat,2);
            Data.Norm = true;
            Data.NormPar1 = par1;
            Data.NormPar2 = par2;
            Data.NormPar3 = par3;
            Data.iMethod = iMethod;            
            N = size(data,1);
            % add all data sequence by sequence (+ normalization)
            for(k=1:N)
                XMat = reshape(permute(data(k,:,:),[2 1 3]),size(data,2),size(data,3)); 
                Data = adddata(Data,XMat);
            end
        end
    end
end
model1.Data    = Data;
model1.Labels  = labels;
model1.ProtoLabels  = [];
model1.oMetric = []; % final metric model - dont forget to update

% Model parameters
Par.PI = [];        % priors
Par.A = [];         % transition probabilities
Par.W = [];         % weigts matrix
Par.Y = [];         % projections of latent points to data space
Par.Beta = [];      % inverse variance of the gaussian model in the data space
Par.oMetric       = []; % metric 
Par.fMetricUpdate = @mupdate; % update algorithm for metric parameter
model1.Par = Par;

% Expectation step
model1.RespMat = [];
model1.LLhood = [];
model1.LLhoodxObs = [];

% Latent space
model1.U = []; % sample points in the latent space - controlled by K
model1.L = []; % --- not used ---

% Max time series duration
model1.Tmax = [];

% RBF grid
model1.M = [];  % number of basis functions
model1.MU= [];  % basis function grid - i.g. automatically set
model1.s = [];  % width / sigma of the basis functions - default 1
model1.FI= [];  % output of the linear basis functions Phi(x)

model1.DistXY    = []; % squared distances between X and Y
model1.MaxIter   = []; % Max iterations
model1.Tolerance = []; % to decide when the convergency is reached
model1.InitMeth  = []; % Initialization method: 'pca','rnd'

model1.Verb = true; %Verbose (true by default)
model1.Status = 'created'; %to control whether the model is already trained
model1.Date = date; % date when the model is created.

model1 = class(model1, 'sgtmtt');




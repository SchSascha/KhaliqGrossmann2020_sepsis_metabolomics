function DataNew = adddata(Data, XMat)
% DATANEW = ADDDATA(DATA, XMAT)
% Adds a new data matrix XMAT to DATA structure.  The resulting structure is
% stored in DATANEW.
%
% DATANEW = ADDDATA(XMAT) creates the structure DATANEW with the data matrix
% XMAT.
%
% The fields of DATA and DATANEW structures are:
%   X   : The dataset.  It is a cell array of data matrices: {[N1xD],
%   [N2xD], ... [NNxD]}.
%   D   : Number of features (columns)
%   N   : [N1, N2, ... NN]
%   Norm: {true or false} - Normalized dataset?
%

error(nargchk(1, 2, nargin, 'struct'));

if nargin==1 % normalize one sequence
    [DataNew.X{1}, par1, par2, par3, iMethod] = normdata(Data);
    [DataNew.N(1) DataNew.D] = size(Data);
    DataNew.Norm = true;
    DataNew.NormPar1 = par1;
    DataNew.NormPar2 = par2;
    DataNew.NormPar3 = par3;
    DataNew.iMethod = iMethod;
elseif nargin==2
    if ~isfield(Data,'X') || ~isfield(Data,'D') || ~isfield(Data,'N')
        error('Bad format of the input dataset structure');
    elseif size(XMat,2) ~= Data.D
        error('The number of columns of the dataset structure and the data matrix are mismatched.');
    else
        DataNew = Data;
        if ~DataNew.Norm
            DataNew.X{end+1} = XMat;
        else % take the original normalization parameters 
            if(DataNew.iMethod == 1) 
                DataNew.X{end+1} = bsxfun(@times,bsxfun(@minus,XMat,DataNew.NormPar1),1.0./DataNew.NormPar2);
            else
                'no such normalisation defined - fix code first'
                return;
            end
        end
        DataNew.N(end+1) = size(XMat,1);
    end
end


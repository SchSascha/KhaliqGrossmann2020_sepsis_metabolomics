function [XNorm, par1, par2, par3, iMethod]=normdata(X, vmax, vmin,iMethod)
% Data normalization:
%   (option 1) by removing of the mean and then setting the max-abs
% value to vmax:
%       [XNorm, maxX, meanX] = normdata(X, vmax) 
% -or-
%       XNorm = normdata(X,vmax)
%
% if vmax is omitted then vmax is assumed to 1
%
%   (option 2) by scaling the dataset to vmax - vmin
%       [XNorm, maxX, minX] = normdata(X, vmax, vmin)
%
% -or-
%       XNorm = normdata(X,vmax, vmin)
%
% Author: Ivan Olier, 2007
% Last Modification: Nov 2/2007

nreg=size(X,1);
XNorm = 0;
par1  = 0;
par2  = 0;
par3  = 0;
iMethod = 1;
if(~exist('iMethod','var'))
    iMethod = nargin;
end

switch(iMethod)
    case 1        
        par1 = mean(X);
        par2 = std(X);
        XNorm = bsxfun(@rdivide,bsxfun(@minus,X,par1),par2);        
    case 2
        meanX = mean(X);
       % size(meanX)
        XNorm = X - repmat(meanX,nreg,1);
        maxX = max(abs(XNorm));
        XNorm = XNorm./repmat(maxX,nreg,1)*vmax;
        par1 = meanX;
        par2 = maxX;
        par3 = vmax;
    case 3
        minX = min(X);
        maxX = max(X);
        ratio = (vmax-vmin)./(maxX-minX);
        XNorm = (X-repmat(minX,nreg,1)).*repmat(ratio,nreg,1) + vmin;
        par1  = minX;
        par2  = ratio;
        par3  = vmin;
    otherwise 
        error('The number of input arguments is wrong');        
end

if nargout >5, error('The number of output arguments is wrong'); end;

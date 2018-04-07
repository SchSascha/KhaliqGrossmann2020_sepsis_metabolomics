% Distance method for generic metric functions
% first  return value is the squared, scaled euclidean distance of X to Y
% second return value is the dim wise distance without metric scaling and unsquared
function [dValue,mDistDim] = distance(cClass, X, Y)      
    dValue = 0;
    if(isempty(cClass.Par))
        warning('Lambda will be initialized to default');       
        [~,D]       = size(Y);   %take dimension from second argument
        S.type      = '.';
        S.subs      = {'Par'};
        L=ones(D,1)./D;
        L=L./sqrt(sum(L.^2));
        cClass      = builtin('subsasgn',cClass,S,L);
    end
        
    %dValue  = sdist(Y,X)';
    %mDistDim= 0;
    %return;
      
    vLam    = cClass.Par;
    %vLam = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0;]';
    %cClass.Par = vLam;
    %vLam(:) = rand(length(vLam),1)*1e-6;
    %vLam([6,7,28]) = 1;
    %vLam    = vLam./sqrt(sum(vLam.^2));
    if(nargout==2)
        [dValue,mDistDim] = sdistNeu(X,Y,vLam);
    else
        [dValue]          = sdistNeu(X,Y,vLam);
        mDistDim          = [];
    end
return

function [DIST,mDistDim] = sdistNeu(W,X,vLam)
% Calculate the squared distances between two sets of data points. 
%
%		This function calculates distances between all data 
%		points in the two data sets T and Y and returns
%		them in a matrix.
%
% Synopsis:	[DIST, minDist, maxDist] = gtm_dist(T, Y, m)
% 		[DIST] = gtm_dist(T, Y)
%
% Arguments:	T, Y -	data set matrices in which each row is a
%			data point; dimensions N-by-D and K-by-D
%			respectively
%
%		m - 	mode of calculation; iff m > 0, min- and  
%			maxDist (below) are calculated; the 
%			default mode is 0 
%
% Return:	DIST -	matrix containing the calculated distances;
%			dimension K-by-N; DIST(k,n) contains the
%			squared distance between T(n,:) and Y(k,:).
%
%		minDist,
%		  maxDist -	vectors containing the minimum and 
%				maximum of each column in DIST, 
%				respectively; 1-by-N; required 
%				iff m > 0.
%			
% Notes:	This m-file provides this help comment and a MATLAB
%		implementation of the distance calculation. If, 
%		however, a mex-file with the same name is present
%		in the MATLABPATH, this will be called for doing
%		the calculation. As this is a computationally 
%		demanding step of the algorithm, an efficient
%		mex-file implementation will improve the performance 
%		of the GTM training algorithm.
%
% See also:	gtm_dstg
%

% Version:	The GTM Toolbox v1.0 beta
%
% Copyright:	The GTM Toolbox is distributed under the GNU General Public 
%		Licence (version 2 or later); please refer to the file 
%		licence.txt, included with the GTM Toolbox, for details.
%
%		(C) Copyright Markus Svensen, 1996

% if (nargin == 3)
%   if (m < 0)
%     error(['Invalid value for mode: ', num2str(m)]);
%   else
%     mode = m;
%   end
% elseif (nargin == 2)
%   mode = 0;
% else
%   error('Wrong number of input arguments');
% end
% 
% if (mode > 0)
%   if (nargout < 3)
%     error('Calculation mode > 0 requires 3 output arguments!');
%   end
% end

%[N, tD] = size(T);
%[K, yD] = size(Y);

% if (yD ~= tD)
%   error('Mismatch in number of columns between T and Y.');
% end

% Summing over components of matrices, we can make use of a 'matrix 
% version' of the identity: (a - b)^2 == a^2 + b^2 - 2*a*b, 
% which yields faster execution. When T and Y consist of single columns 
% (which may be the case with calls from gtm_gbf), this has to be handled 
% as a special case. 
% if yD > 1
  %% Modification to accept missing data coded as NaN %%
%   Known = ~isnan(T);
%   T(find(~Known)) = 0; 
  %% end modif.

% if nargin==1
%     X=g.Data;
%     W=g.W;
% end

[N D] = size(X);
[K D] = size(W);
mDistDim =  [];
% lambda is squared as well
vLam = vLam.^2;
if D>1
    %DIST = W * diag(vLam) * X';
    DIST = W*bsxfun(@times,vLam',X)';
    %%xx = sum((X.^2 * diag(vLam))');
    xx = sum(bsxfun(@times,X.^2,vLam'),2)';
    %%ww = sum((W.^2 * diag(vLam))');
    ww = sum(bsxfun(@times,W.^2,vLam'),2)';
    %%DIST = ww'*ones(1,N) + ones(K,1)*xx - 2*DIST;
    DIST = bsxfun(@minus,bsxfun(@plus,ww,xx'),2*DIST')';
    if(nargout==2)
        mDistDim = zeros(K, N, D);
        for(k=1:K)
            mDistDim(k,:,:) = bsxfun(@minus,X,W(k,:));
        end
    end
else
    DIST     = zeros(K, N);
    mDistDim = zeros(K, N, D);
    for n=1:N
        for k=1:K
            mDistDim(k,n,:) = X(n,:) - W(k,:);              % plain distance unscaled and unsquared
            DIST(k,n) = sum((vLam .* mDistDim(k,n,:).^2));  % lambda is already squared - see l:115
        end
    end
end

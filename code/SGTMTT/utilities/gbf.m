function FI = gbf(MU, sigma, X, varargin)

% Calculates the output of Gaussian basis functions for a given set of input
%
% Synopsis:	FI = gtm_gbf(MU, sigma, X)
%
% Arguments:	MU - 	a M-by-L matrix containing the centers of the 
%			basis functions
%
%		sigma - a scalar giving the standard deviation of the
%			radii-symmetric Gaussian basis functions,
%
%		X -	the latent variable sample forming the set of
%			inputs; K-by-L
%	
% Return:	FI -	the matrix of basis functions output values;
%			K-by-(M+1), "+1" for a bias basis function
%			with a fixed value of 1.0
%
% See also:	gtm_lbf	
%

% Version:	The GTM Toolbox v1.0 beta
%
% Copyright:	The GTM Toolbox is distributed under the GNU General Public 
%		Licence (version 2 or later); please refer to the file 
%		licence.txt, included with the GTM Toolbox, for details.
%
%		(C) Copyright Markus Svensen, 1996

[K L] = size(X);
[M L2] = size(MU);

if (L~=L2 | size(sigma)~=[1 1])
  error('Mismatch in dimensions of input argument matrices - likely latent space and rbf space are incompatible.'); % You cannot have, say, a 2D latent space AND a 3D BFs grid 
end

% Calculate outputs of the Gaussian basis functions
if(isempty(varargin))
    DIST = sdist(MU, X); % Calulate distances between all points of MU and X
else
    oMetric = varargin{1};
    DIST    = oMetric.distance(MU, X);
end
FI = exp((-1/(2*sigma^2))*DIST);

% Add bias basis function
FI = [FI' ones(K,1)];

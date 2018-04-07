function grid = rctg(xDim, yDim)
% Produces a 2D grid with points arranged in a rectangular lattice.
%
%		The grid is centered on the origin and scaled so the 
%		dimension (X or Y) with largest number of points 
%		ranges from -1 to 1.
%
% Usage (1st form):
%               g = rctg(g)
%
% where: g - gtm class
%
% Usage (2nd form):
%               grid = rctg(xDim, yDim)
%
% Input arguments:	
%               xDim, yDim -	number of points along the X and Y dimensions; must be >=2.
%
% Output argument:	
%               grid -	a (xDim*yDim)-by-2 matrix of grid points with the first point being 
%                       the top-left corner and subsequent points following column-wise.
%
% This code is based on 'gtm_rctg' function written by M. Svensen 1996
% Author: Ivan Olier, 2007
% Last modification: mar 31/2007
% 

if nargin>2
    error('The number of arguments is wrong');
elseif nargin==1
    xDim=sqrt(xDim);
    yDim=xDim;
end

if (xDim<2 | yDim<2 | (yDim ~= fix(yDim)) | (xDim ~= fix(xDim)))
  error(['Invalid grid dimensions ', ...
		40, num2str(xDim), 44, num2str(yDim), 41, 46]); % 40,44,41,46 are ascii for (,). respectively
end

% Produce a grid with the right number of rows and columns
[X, Y] = meshgrid([0:1:(xDim-1)], [(yDim-1):-1:0]);

% Change grid representation from mesh to vector - just to get position of
% the rbf functions
grid = m2r(X, Y);

% Scale grid to correct size (max.abs. value 2, so that when centered around 0, it spans from -1 to 1)
maxVal = max(max(grid));
grid = grid*(2/maxVal);

% Shift grid to correct position ( centering around 0 )
maxXY= max(grid);
grid(:,1) = grid(:,1) - maxXY(1)/2;
grid(:,2) = grid(:,2) - maxXY(2)/2;






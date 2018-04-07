function  I = argmax(X)

%ARGMAX -- Index of the maximal element.
%
%  I = ARGMAX(X)   is equivalent to  [_,I] = max(X)
%
%  Example:  
%    argmax([4 5 8 1 6 3  9 2 7 3 6 1])  ==>  7
%    argmax([4 5 8 1 6 3; 9 2 7 3 6 1])  ==>  [2 1 1 2 1 1]
%
%  See also MAX, ARGMIN.

%  Original coding by Alexander Petrov, Ohio State University.
%  $Revision: 1.0 $  $Date: 2004/03/03 09:20 $
%
% Part of the utils toolbox version 1.1 for MATLAB version 5 and up.
% http://alexpetrov.com/softw/utils/
% Copyright (c) Alexander Petrov 1999-2006, http://alexpetrov.com
% Please read the LICENSE and NO WARRANTY statement in ../utils_license.m

[foo,I] = max(X) ;

%%  Return I
%%%%%% End of file ARGMAX.M

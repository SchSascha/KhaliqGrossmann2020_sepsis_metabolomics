function [mod indx]=pmd(RespMat,model)
% Calculate the mode projections
% gives the mapped grid positions in the latent space
% for each point in the original data space based on
% the responsibilities
% - useful for plots
%   mod=pmd(RespMat)
if(nargin==1)
    [mod indx]=pmd_std(RespMat)
else
    [mm, indx]=max(RespMat);
    U   = get(model,'u')
    mod = U(indx',:);
end


function [mod indx]=pmd_std(RespMat)
%Calculate the mode projections
%   mod=pmd(RespMat)

[mm, indx]=max(RespMat);
U = rctg(size(RespMat,1));
mod = U(indx',:);
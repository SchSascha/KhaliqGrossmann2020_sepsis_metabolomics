% Metric update according to the Soft Robust LVQ cost function
% DistPerTime Protos x Time x Points x Dims
% RespPerTime Protos x Time x Points
function mupdate(Parameters,X,X_Labels,Y_Labels,Resp,distXYDim,vProb,eta,bFlag)
oMetric = Parameters.oMetric;
beta    = Parameters.Beta;
[N, D] = size(X);
K      = numel(Y_Labels);
L     = oMetric.Par;

% robust soft metric learning update
lab_eq = bsxfun(@eq, X_Labels, Y_Labels');
lab_neq=~lab_eq;
sum_R  = sum(lab_eq.*Resp);

% sum over r_kn on nth column
% over k for those latent points with equal label to the data only
%    
R_y=bsxfun(@rdivide, Resp, sum_R);
R_y(isnan(R_y))=0;
R_y(R_y==inf)=0;

R_eq_tmp=zeros(K, N);
R_eq_tmp(lab_eq)=(R_y(lab_eq)-Resp(lab_eq));

R_neq_tmp=zeros(K, N);
R_neq_tmp(lab_neq)=Resp(lab_neq);   
%Tmp = DistPerTime;
%distXYDim = reshape(Tmp,size(Tmp,1), size(Tmp,2) * size(Tmp,3),size(Tmp,4) ) ;
distXYDim = permute(distXYDim.^2,[2,3,1]);
if(bFlag)
    R_eq_tmp = reshape(R_eq_tmp,K,size(Resp,1) / size(distXYDim,1),size(distXYDim,1));
    R_eq_tmp = squeeze(sum(R_eq_tmp,2));
    R_neq_tmp = reshape(R_neq_tmp,K,size(Resp,1) / size(distXYDim,1),size(distXYDim,1));
    R_neq_tmp = squeeze(sum(R_neq_tmp,2));    
end
mMat      = bsxfun(@times, permute(R_eq_tmp - R_neq_tmp,[2,3,1]),bsxfun(@minus, (1./L)',beta * bsxfun(@times, L', distXYDim)));   
for(k=1:length(vProb))
    mMat(:,:,k) = vProb(k) * mMat(:,:,k);
end
L_tmp = sum(sum( mMat,3),1);

% normalization
L=max(1e-6,L+eta*L_tmp'); 
L=L./sqrt(sum(L.^2));
oMetric.Par = L;

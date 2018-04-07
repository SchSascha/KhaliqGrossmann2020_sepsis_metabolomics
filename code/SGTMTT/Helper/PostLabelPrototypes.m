function [mGamma,PrototypeLabels,dAcc,P,vLik,vTimeRelevance] = PostLabelPrototypes(oModel1,oModel2,X1,XR,L1,LR);
    N1     = size(X1,1);
    NR     = size(XR,1);
    T      = size(X1,2);
    D      = size(X1,3);
    Lik1_against_1 = zeros(1,N1);
    Lik1_against_2 = zeros(1,N1);
    Lik2_against_1 = zeros(1,NR);
    Lik2_against_2 = zeros(1,NR); 
    mLhoodV1a1= zeros(N1,T);
    mLhoodV1a2= zeros(N1,T);
    mLhoodV2a1= zeros(NR,T);
    mLhoodV2a2= zeros(NR,T);

    K       = get(oModel1,'K');
    mResp   = zeros(K*2,(N1+NR)*T); % store the responsibilities for each point for all protototypes (first K labeled 1, second K labeled 2 or R)  
    mGamma  = mResp;
    for(k=1:N1),         
        [Lik,Res,~,LhoodV1]         = test(oModel1,reshape(permute(X1(k,:,:),[2,1,3]),T,D));         
        mLhoodV1a1(k,:)   = log(LhoodV1); % log-likelihoods per time point for this sample
        iStart            = (k-1)*T+1;
        mResp([1:K],[iStart:iStart+T-1]) = Res{1};
        mGamma([1:K],[iStart:iStart+T-1]) = Res{1};
        Lik1_against_1(k) = Lik;
        [Lik,Res,~,LhoodV2]         = test(oModel2,reshape(permute(X1(k,:,:),[2,1,3]),T,D)); 
        mLhoodV1a2(k,:)   = log(LhoodV2);
        mResp([K+1:end],[iStart:iStart+T-1]) = -Res{1};        
        mGamma([K+1:end],[iStart:iStart+T-1]) = Res{1};
        Lik1_against_2(k) = Lik;
    end    
    for(k=1:NR),         
        [Lik,Res,~,LhoodV1]         = test(oModel1,reshape(permute(XR(k,:,:),[2,1,3]),T,D));         
        mLhoodV2a1(k,:)   = log(LhoodV1);
        iStart            = N1*T + (k-1)*T+1;
        mResp([1:K],[iStart:iStart+T-1]) = mResp([1:K],[iStart:iStart+T-1])-Res{1};
        mGamma([1:K],[iStart:iStart+T-1]) = Res{1};
        Lik2_against_1(k) = Lik;
        [Lik,Res,~,LhoodV2]         = test(oModel2,reshape(permute(XR(k,:,:),[2,1,3]),T,D)); 
        mLhoodV2a2(k,:)   = log(LhoodV2);
        mResp([K+1:end],[iStart:iStart+T-1]) = mResp([K+1:end],[iStart:iStart+T-1])+Res{1}; 
        mGamma([K+1:end],[iStart:iStart+T-1]) = Res{1}; 
        Lik2_against_2(k) = Lik;
    end; 
    vTimeRelevance = 0.5 * (mean((mLhoodV1a2 - mLhoodV1a1 <= 0)) + mean(mLhoodV2a1 - mLhoodV2a2 <= 0) );
    vLik = [Lik1_against_1 Lik2_against_2];
    dAcc = (sum(Lik1_against_2 - Lik1_against_1 <= 0) + sum(Lik2_against_1 - Lik2_against_2 <= 0)) / (N1+NR);
    dAcc
    X_resp                  = [sum(mResp(:,[1:N1*T]),2),sum(mResp(:,[N1*T+1:end]),2)];
    % normalize
    X_resp                  = bsxfun(@minus,X_resp,min(X_resp(:)));
    [~, PrototypeLabels]    = max(X_resp, [], 2);
    vNormalizer             = sum(X_resp,2);
    vAbsoluteResponsibility = vNormalizer;
    vAbsoluteResponsibility = vAbsoluteResponsibility./max(vAbsoluteResponsibility(:)); % which prototype gots high responsibility
    X_resp                  = X_resp./vNormalizer(:,ones(1,size(X_resp,2)));
    % get weightings of the responsibilities
    R=sort( X_resp,2,'descend');    % sort responsibilities for each prototype
    P=abs(R(:,1)-R(:,2));           % get absolute distance of the first two largest responsibilites - useful if #classes >2
    P=P./sum(P);                    % normalize to 1
    P=vAbsoluteResponsibility .* P; % prob of being responsible & prob of having good separation

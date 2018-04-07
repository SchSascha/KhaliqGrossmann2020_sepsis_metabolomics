function [dAcc] = PredictSGTMTT(sModel,X,L)
    N     = size(X,1);
    T     = size(X,2);
    D     = size(X,3);
    C     = length(sModel); % number of classes
    vLikelihoods = zeros(C,1);
    dAcc  = 0;
    if(C > 1)
        for(k=1:N),         
            vLabeling = zeros(C,1);
            for(c=1:C)
                sCurrentModel = sModel(c).Models.oModel1;
                vLabeling(c)  = sModel(c).Models.vLabels(1);
                [Lik,Res]  = test(sCurrentModel,reshape(permute(X(k,:,:),[2,1,3]),T,D));         
                vLikelihoods(c) = Lik;
            end
            dAcc=dAcc+(L(k,1)==vLabeling(argmax(vLikelihoods)));
        end
    else
        sCurrentLabels = sModel.Labels;
        for(k=1:N),         
            sCurrentModel = sModel(1).Models.oModel1;
            [Lik1,Res1] = test(sCurrentModel,reshape(permute(X(k,:,:),[2,1,3]),T,D));         
            sCurrentModel = sModel(1).Models.oModelR;
            [Lik2,Res2] = test(sCurrentModel,reshape(permute(X(k,:,:),[2,1,3]),T,D));                         
            vLikelihoods(1) = Lik1;
            vLikelihoods(2) = Lik2;                       
            dAcc=dAcc+(L(k,1)==sCurrentLabels(argmax(vLikelihoods)));
        end
    end
    dAcc = dAcc./N;

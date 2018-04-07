%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supervised Generative Topographic Mapping - Through Time (SGTM-TT)
%
% (C) Frank-Michael Schleif (2010-2012) - SGTM-TT
% fschleif@techfak.uni-bielefeld.de
%
% and Ivan Olier (2006-2009) - GTM-TT implementation
% ivan.olier@manchester.ac.uk
%
% GTM-TT is original proposed by
% I. G. D. Strachan
% IGD.Strachan@gmail.com
% 
% for GTM see:
% @article{DBLP:journals/ijon/BishopSW98,
%   author    = {Christopher M. Bishop and
%                Markus Svens{\'e}n and
%                Christopher K. I. Williams},
%   title     = {Developments of the generative topographic mapping},
%   journal   = {Neurocomputing},
%   volume    = {21},
%   number    = {1-3},
%   year      = {1998},
%   pages     = {203-224},
%   ee        = {http://dx.doi.org/10.1016/S0925-2312(98)00043-5},
%   bibsource = {DBLP, http://dblp.uni-trier.de}
% }
%
% for GTM-TT see:
% GTM Through Time, Proceedings IEE Fifth International Conference on
%   Artificial Neural Networks, Cambridge, U.K., (1997) 111-116
% 
% @article{DBLP:journals/nn/OlierV08,
%   author    = {Iv{\'a}n Olier and
%                Alfredo Vellido},
%   title     = {Advances in clustering and visualization of time series
%                using GTM through time},
%   journal   = {Neural Networks},
%   volume    = {21},
%   number    = {7},
%   year      = {2008},
%   pages     = {904-913},
%   ee        = {http://dx.doi.org/10.1016/j.neunet.2008.05.013},
%   bibsource = {DBLP, http://dblp.uni-trier.de}
% }
%
% for SGTM-TT see:
% @inproceedings{DBLP:conf/ijcnn/SchleifGH12,
%   author    = {Frank-Michael Schleif and
%                Andrej Gisbrecht and
%                Barbara Hammer},
%   title     = {Relevance learning for short high-dimensional time series
%                in the life sciences},
%   booktitle = {IJCNN},
%   year      = {2012},
%   pages     = {1-8},
%   ee        = {http://dx.doi.org/10.1109/IJCNN.2012.6252653},
%   crossref  = {DBLP:conf/ijcnn/2012},
%   bibsource = {DBLP, http://dblp.uni-trier.de}
% }
% @proceedings{DBLP:conf/ijcnn/2012,
%   title     = {The 2012 International Joint Conference on Neural Networks
%                (IJCNN), Brisbane, Australia, June 10-15, 2012},
%   booktitle = {IJCNN},
%   publisher = {IEEE},
%   year      = {2012},
%   isbn      = {978-1-4673-1488-6},
%   ee        = {http://ieeexplore.ieee.org/xpl/mostRecentIssue.jsp?punumber=6241467},
%   bibsource = {DBLP, http://dblp.uni-trier.de}
% }
%
% mData   - Data : N x T x D 
% vLabels - Labels (1,2,...) N x 1
% K - number of states (for each GTM)
% M - number of basis functions to calculate Y = Phi * W 
% varargin - {metric(true),functional(true),iter(100)} - provide all if used
%
% Sample call:
% [vModels] = GTMTT_Supervised(X,L,K,4);
%
% see also Demo.m and Demo2.m 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vModels,vUnique] =  GTMTT_Supervised(mData,vLabels,K,M,varargin)
    % for experiments
    %s = RandStream('mt19937ar','Seed',1);
    %RandStream.setGlobalStream(s);
    % inits 
    vUnique             = unique(vLabels);
    D                   = size(mData,3);    
    T                   = size(mData,2);    % assumes equal number of time points for each point
    iFolds              = numel(vUnique);   % number of 1 vs rest folds
    eta0                = 0.1;              % learn rate for metric update
	iSamplingSize       = 5  ;              % for the distance calculation use each (1), or each third (3) point ...
    Deff                = D+1;    
    iterMin             = 10;               % minimal number of cycles
    iterMax             = 100;
    bWithMetricLearning = true;             % set this to false to skip metric adaptation
    bFunctional         = false;             % used metric (euclidean (false) or lee-functional (true)  
    bThresholdFeatures  = false;             % automatic thresholding of the features in late learning
    if(nargin>4)
        iterMax             = varargin{1};
        bWithMetricLearning = varargin{2};
        bFunctional         = varargin{3};
    end
    if(iFolds == 2)
        iFolds = 1;
    end
    
    % normalize the data and store the normalization info
    [~, par1, par2, par3, iMethod] = normdata(reshape(permute(mData,[2 1 3]),size(mData,1) * size(mData,2),size(mData,3))); % get normalization constants
    oNorm = {par1, par2, par3, iMethod};
    
    % for each 1 vs rest scheme
    for(i=1:iFolds)
        vLabel          = vLabels; % backup will be changed below        
        iOriginalLabelClass1      = vUnique(i);        
        vCurrentEntries = (vLabel == vUnique(i));
        N               = length(vCurrentEntries);
        vCurrentEntriesR= (vLabel ~= vUnique(i));        
        L1              = ones(sum(vCurrentEntries),1);
        X1              = mData(vCurrentEntries,:,:);
        LR              = 2*ones(sum(~vCurrentEntries),1);
        XR              = mData(vCurrentEntriesR,:,:);
        vLabel          = [L1;LR];
        oModel1         = sgtmtt(K,X1,L1,oNorm);         % define the model
        oModel1         = set(oModel1,'m',M);      % number of basis functions
        oModelR         = sgtmtt(K,XR,LR,oNorm);         % define the rest model
        oModelR         = set(oModelR,'m',M);      % number of basis functions
        % set common metric parameters (global)
        oMetric         = smetric;
        oModelR         = set(oModelR,'oMetric',oMetric); % thats a handle - so both have the same metric object
        oModel1         = set(oModel1,'oMetric',oMetric); % thats a handle - so both have the same metric object
        oModel1         = set(oModel1,'deff',D); % effective dimension (=D without relevance)
        oModelR         = set(oModelR,'deff',D); %
        XMatFormat1 = cell2mat(get(oModel1,'x')');
        XMatFormatR = cell2mat(get(oModelR,'x')');
        % (!) block order
        XMat        = [XMatFormat1;XMatFormatR];        
        XPerTime    = reshape(XMat,T,N,D);                
        bMetricLearning = false;    % will be activate later (see below)  - deactivated for the first 10 cycles
        for(iter=1:iterMax)            
            % train both models
            oModel1     = train_one_step(oModel1,iter); %training for one step
            oModelR     = train_one_step(oModelR,iter); %training for one step
            if(iter==1) % take some inits for fix
                oMetric = get(oModel1,'oMetric'); % take the initialized metric back
                Deff    = sum(oMetric.Par); % effective dimension incorporating scaling                
                oModel1 = set(oModel1,'deff',Deff); % update effective dimension
                oModelR = set(oModelR,'deff',Deff); %
            end
            % evaluate model convergence
            bConverged1 = CheckConvergence(oModel1);
            bConverged2 = CheckConvergence(oModelR);
            if(bConverged1 && bConverged2 && iter>iterMin || iter == iterMax)
                break;
            end
            
            if(bWithMetricLearning && iter > 10)
                bMetricLearning = true; % not yet
            end  
            
            % simple probabilistic max-likelihood classification
            [mProtos,~,dAcc,~,vLik] = PostLabelPrototypes(oModel1,oModelR,X1,XR,L1,LR);                
            dAcc
            Y1          = get(oModel1,'y');
            YR          = get(oModelR,'y');                
            % Postlabeling - will also merge the two models   
            
            if(bMetricLearning)
                vRandPerm       = randperm(N);
                oLeeMetric      = metric_lee;
                iNrOfUpdates    = ceil(N/iSamplingSize);
                RespPerTime     = reshape(mProtos,K*2,T,N);    % Protos x Time x Points
                mDimDistLee1    = zeros(iNrOfUpdates,D);
                mDimDistLee2    = zeros(iNrOfUpdates,D);
                for(k=1:iNrOfUpdates) % for each x point for each dimension
                    iCurrent = vRandPerm(k);
                    vRespK   = argmax(RespPerTime([1:K],:,iCurrent));     % responsibilities with respect to the first (true) model
                    vRespKp  = argmax(RespPerTime([K+1:end],:,iCurrent)); % responsibilities with respect to the alternative (rest) model
                    
                    % fst model - correct class
                    if(vLabel(iCurrent) ~= 1)
                        mRecSignals2 = Y1(vRespK,:)';
                        mRecSignals1 = YR(vRespKp,:)';                    
                    else
                        mRecSignals1 = Y1(vRespK,:)';
                        mRecSignals2 = YR(vRespKp,:)';                    
                    end
                    
                    if(bFunctional)
                        for(d=1:D),                        
                            mDimDistLee1(k,d) = sum(distance(oLeeMetric,mRecSignals1(d,:),XPerTime(:,iCurrent,d)'));
                            mDimDistLee2(k,d) = sum(distance(oLeeMetric,mRecSignals2(d,:),XPerTime(:,iCurrent,d)'));
                        end
                    else
                        for(d=1:D),                        
                            mDimDistLee1(k,d) = sum((mRecSignals1(d,:)-XPerTime(:,iCurrent,d)').^2);
                            mDimDistLee2(k,d) = sum((mRecSignals2(d,:)-XPerTime(:,iCurrent,d)').^2);
                        end                        
                    end
                end
                clear mRecSignals; % free space
                
                % try to discriminate
                ranking                             = zeros(size(X1,1),1);
                [~,ILiks]                           = sort(vLik(1:size(X1,1)),'descend');
                ranking(ILiks)                      = [0:length(ILiks)-1];
                vWeights                            = exp(-ranking/(length(ranking)*0.5));

                ranking                             = zeros(length(vLik) - size(X1,1),1);
                [~,ILiks]                           = sort(vLik(size(X1,1)+1:end),'descend');
                ranking(ILiks)                      = [0:length(ILiks)-1];
                vWeights(size(X1,1)+1:length(vLik)) = exp(-ranking/(length(ranking)*0.5));
                vUpdateIndices                      = vRandPerm(1:iNrOfUpdates);
     
                % calculate the metric update in acc. to the margin crit of the GLVQ cost function
                DJ=sum(bsxfun(@times,mDimDistLee1,oMetric.Par'),2);
                DI=sum(bsxfun(@times,mDimDistLee2,oMetric.Par'),2);
                PSIP = bsxfun(@times,2*DI,1.0./(DI+DJ).^2);
                PSIM = bsxfun(@times,2*DJ,1.0./(DI+DJ).^2);                
                V= (sum(bsxfun(@times,vWeights(vUpdateIndices), -(bsxfun(@times,PSIP,mDimDistLee1) - bsxfun(@times,PSIM, mDimDistLee2))))./sum(vWeights(vUpdateIndices)))';
                Vo = oMetric.Par;

                Vo=max(1e-6,Vo+eta0*V);                 
                if(bThresholdFeatures && iter/iterMin > 0.8) % at 50%
                    vIRelevant = find(Vo>(mean(Vo) + std(Vo)));
                    if(isempty(vIRelevant))
                        iRequestedFeatures = 7;
                        [~,vI] = sort(Vo,'descend');
                        Vo(vI(iRequestedFeatures+1:end)) = 0;           % remove irrelevant entries
                    else
                        vTemp = Vo(vIRelevant);
                        Vo(:) = 0;
                        Vo(vIRelevant) = vTemp;
                    end
                    bWithMetricLearning = false; % 
                    bMetricLearning     = false; % metric learning switched of for next step
                end
                Vo=Vo./sqrt(sum(Vo.^2));
                oMetric.Par = Vo;
                oModelR = set(oModelR,'oMetric',oMetric); % thats a handle - so both have the same metric object
                oModel1 = set(oModel1,'oMetric',oMetric); % thats a handle - so both have the same metric object
                Deff    = sum(oMetric.Par); % effective dimension incorporating scaling                
                if(floor(Deff)<=1)
                    Deff = D;
                end    
                Deff   = 1+Deff;
            end
            
            % [ BETA update - and synchronization]
            % we also have to care our self for a modified beta update
            % (the internal one is switched off now)
            [DistXY] = oMetric.distance([Y1;YR], XMat); 
            iPoints1 = size(X1,1) * size(X1,2);
            dBeta1 = size(XMatFormat1,1)*Deff/sum(sum(DistXY(:,[1:iPoints1]).*mProtos(:,[1:iPoints1])));          % plausible 
            dBetaR = size(XMatFormatR,1)*Deff/sum(sum(DistXY(:,[iPoints1+1:end]).*mProtos(:,[iPoints1+1:end])));  % plausible 
            dBeta  = 0.5*(dBeta1 + dBetaR);
            oModel1 = set(oModel1,'beta',dBeta);
            oModelR = set(oModelR,'beta',dBeta);                        

            oModel1         = set(oModel1,'deff',Deff); % update effective dimension
            oModelR         = set(oModelR,'deff',Deff); %
            
            if(bMetricLearning && mod(iter,10)==0)
                oMetric.Par'
            end
        end
        [~,~,dAcc,~,~,vTimeRelevance] = PostLabelPrototypes(oModel1,oModelR,X1,XR,L1,LR);     
        dAcc

        sModel = struct('oModel1',oModel1,'oModelR',oModelR,'vLabels',[iOriginalLabelClass1 -1],'vTimeRelevance',vTimeRelevance);
        vModels{i} = sModel;  
        if(iFolds==1)
            vModels{i+1} = struct('oModel1',oModelR,'vLabels',[vUnique(i+1) -1],'vTimeRelevance',vTimeRelevance);
        end
    end
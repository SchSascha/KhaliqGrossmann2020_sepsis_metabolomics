function g=train(g)
% Train the GTM-TT model
% modifications: Frank-Michael Schleif, March-June, 2011

if ~strcmp(g.Status,'initialized')
    if g.Verb
        disp('MSG: Initializing the model...');
    end
    g = init(g);
end

FI_T       = g.FI';
N          = g.Data.N;
NFullData  = sum(N);
XMatFormat = cell2mat(g.Data.X');

if isempty(g.MaxIter) 
    MaxIter = intmax;
else
    MaxIter = g.MaxIter;
end
MaxIter = 500;
DataLabels      = g.Labels;
eta0            = 1e-5; % learn rate for metric update
bLearnMetric    = false; 
[M,D] = size(g.Par.W);
Deff = D;

bFigure = true;
for iter=1:MaxIter
    % Forward/backward procedure (E-step)        
    % calculate emission probabilities
    dBeta = g.Par.Beta;
    BjnMat = ((dBeta/(2*pi))^(Deff/2)*exp(-0.5*dBeta*g.DistXY));
    Bjn = mat2cell(BjnMat,g.K,N);
    
    % HMM forward-backward algorithm
    % Gamma - resp. matrix, gamma-init - for PI reestimate , Xi - expected number of transitions
    [Gamma, GammaInit, Xi, LLhood, LLhoodV]=forwback(g.Par.PI, g.Par.A, Bjn, N);
    g.RespMat = Gamma;    
    if g.Verb
        fprintf('It: %g\tlogLH: %g\tBeta: %g\n', iter, LLhood, g.Par.Beta);
    end
    g.LLhood(iter) = LLhood;
    g.LLhoodxObs(iter,:) = LLhoodV;
    
    if(0)
        if iter == 1 LLini = LLhood;
        elseif (LLhood-LLini)<(1 + g.Tolerance)*(g.LLhood(iter-1)-LLini) && iter>100
            if g.Verb fprintf('MSG: The algorithm has reached the convergence.\n'); end;
            g.Status = 'trained';
            break;
        elseif ~isfinite(LLhood)
            if g.Verb fprintf('MSG: The training has been early stopped becouse the log-likelihood tends to infinite\n'); end;
            g.Status = 'failed';
            break;
        end
    end
    
    if iter==MaxIter
        if g.Verb fprintf('MSG: The algorithm has not reached the convergence.\n'); end;
        g.Status = 'undertrained';
        break;
    end
    
    %M-step    
    g.Par.A = Xi./repmat(sum(Xi,2),[1 g.K]);    % re-estimate of the transition prob matrix (just normalize Xi)
    g.Par.A(isnan(g.Par.A))=0;                  % cleanup
    g.Par.PI = GammaInit/sum(GammaInit);        % update priors (just normalize gammainit)   
    
    %Resolve equation Phi'*Gs*Phi*W' = Phi' * Rs *X
    % Rs - the responsibility matrix r_kn = p_nk' / sum p_nk'
    % Gs - diag(sum(Rs))
    % phis (M+1 x K) * Gs (K x T x N) * phis (K x M+1) * weights (D x M+1) == phis (M+1 x K) * Rs * Data (T x N x D)
    %Update W
    GammaMat = cell2mat(Gamma);
    
    %  Phi'*G*Phi 
    A1 = full(FI_T*spdiags(sum(GammaMat')', 0, g.K, g.K)*g.FI);
    % get upper triangle matrix R such that R * R' = A1
    [cholDcmp singular] = chol(A1);

    if singular
        if g.Verb
            fprintf('WARNING: M-Step matrix singular, using pinv.\n');
        end
        g.Par.W = pinv(A1)*(FI_T*(GammaMat*XMatFormat));
    else
        % R' * R * W' = Phi' * Rs *X  
        g.Par.W = cholDcmp \ (cholDcmp' \ (FI_T*(GammaMat*XMatFormat)));
    end

    if(mod(iter,30)==0)
        g.Par.oMetric.Par'
    end
    
    % Calculate Y new
    g.Par.Y    = g.FI*g.Par.W;  % y(x_k, W^), latend points mapped to the orig D space using basis functions FI and updated weights    
    
    % calculate labels    
    X_resp=[];
    for i=1:max(DataLabels)    
      X_resp(:,i)=sum(GammaMat(:,DataLabels==i),2);
    end
    vNormalizer             = sum(X_resp,2);
    vAbsoluteResponsibility = vNormalizer;
    vAbsoluteResponsibility = vAbsoluteResponsibility./max(max(vAbsoluteResponsibility)); % which prototype gots high responsibility

    X_resp      = X_resp./vNormalizer(:,ones(1,size(X_resp,2)));
    [~, PrototypeLabels]=max(X_resp, [], 2);
    if(mod(iter,10)==0)
        dAcc = sum(PrototypeLabels(argmax(GammaMat)) == DataLabels)./length(DataLabels)    
    end
        
    if(iter>10 && i<MaxIter && bLearnMetric) % skip last update
        % get weightings of the responsibilities
        R=sort( X_resp,2,'descend');    % sort responsibilities for each prototype
        P=abs(R(:,1)-R(:,2));           % get absolute distance of the first two largest responsibilites - useful if #classes >2
        P=P./sum(P);                    % normalize to 1
        P=vAbsoluteResponsibility .* P; % prob of being responsible & prob of having good separation
        
        % update metric parameter (distance object is a handle)
        g.Par.fMetricUpdate(g.Par,XMatFormat,DataLabels,PrototypeLabels,GammaMat',mDimDist,P,eta0);     
        Deff   = sum(g.Par.oMetric.Par); % effective dimension incorporating scaling
        if(floor(Deff)<=1)
            Deff = g.Data.D;
        end    
        Deff   = 1+Deff;        
    end
    g.Par.Deff = Deff;
    
    [g.DistXY,mDimDist] = g.Par.oMetric.distance(g.Par.Y, XMatFormat); % includes relevances here
    
    if(iter>10 && bLearnMetric)
        if(1) % re-estimate the beta parameter but scale the influence of the single dimensions, but wait a bit
            g.Par.Beta = betai(bsxfun(@times,(g.FI' * g.Par.Y),g.Par.oMetric.Par')); 
        else % reestimate based on prototype distances in the projected back to data space
            g.Par.Beta = betai(g.FI' * g.Par.Y);
        end
    else
        g.Par.Beta = NFullData*Deff/sum(sum(g.DistXY.*GammaMat));          % plausible 
    end
    % g.Par.Beta = betai(g.FI' * g.Par.Y); % alternatively reestimate
    
    if(0 & mod(iter,10)==0)        %3d
        [Y,I] = sort(g.Par.oMetric.Par,'descend');
        plot_data_grid_2d(XMatFormat, DataLabels, g.Par.Y, PrototypeLabels, M,bFigure,I(1),I(2));
        sprintf('Current best metric 2d [%i = %f %i = %f]',I(1),Y(1),I(2),Y(2))
        if(bFigure)
            bFigure = ~bFigure;
        end
        pause(5)
    end
end
% calculate labels    
X_resp=[];
for i=1:max(DataLabels)    
  X_resp(:,i)=sum(GammaMat(:,DataLabels==i),2);
end
[~, PrototypeLabels]=max(X_resp, [], 2);
sum(PrototypeLabels(argmax(GammaMat)) == DataLabels)./length(DataLabels)
g.ProtoLabels = PrototypeLabels;
g.oMetric     = g.Par.oMetric; % finish metric object
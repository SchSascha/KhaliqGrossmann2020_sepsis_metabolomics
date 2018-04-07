% Distance method derived to lam for generic metric functions
% for the lambda scaling we consider only a global lambda
function vValue = distance_deriv_lam(cClass, vX, vY,varargin) 
    if(isvector(vX))
        iNrDims = length(vX);
        iNrVectors = 1;
    else
        iNrDims    = size(vX,2);
        iNrVectors = size(vX,1);
    end
	if(length(varargin) == 1 && isvector(varargin{1}) && length(varargin{1}) == iNrDims)
        vLam = varargin{1};
    else
        vLam = cClass.vLam;
    end
	if ~(length(vLam) == iNrDims || size(vLam,2) == iNrDims)
        vLam = ones(iNrVectors,iNrDims);        
    end
    if(iNrVectors == 1)
        vLam = vLam(1,:); %  matrix use only first dimension (only global lambda)
    end % since subsequently prototype updates are calculated individual no repmat(lam,i,1) needed

    X_ALL = vX - vY;
    vValue = zeros(iNrVectors,iNrDims);
    C = cClass.tau.^2/2; % remain square due to tau in c_j
    vIndicesWithOutBorder = [2:iNrDims-1];   
    vIndices_1  = [2:iNrDims];
    vIndicesM1  = vIndices_1-1;
    vIndices_2  = [1:iNrDims-1];
    vIndicesP1  = vIndices_2+1;
     %% !!! we are taking here p == 2
    vCp         = 1/2 * ones(iNrVectors,1); %1.0/cClass.p * distance(cClass, vX, vY,0); %% outer derive to lam
    A           = zeros(1,iNrDims-1);
    B           = zeros(1,iNrDims);          
    C_j         = zeros(1,iNrDims);       
    T1          = zeros(1,iNrDims);    
    T2          = zeros(1,iNrDims-1);    
    %% including 
    vLam = vLam(1,:);
    for i=1:iNrVectors                      
        X       = X_ALL(i,:);
        absX    = abs(X);
        XSquare = X.^2;
        I1     = find(0 <= X(vIndices_1) .* X(vIndicesM1)); %
        I1_NEG = find(0 > X(vIndices_1) .* X(vIndicesM1)); % xor is slow setxor(I1,[1:length(vIndicesM1)]);

        I2     = find(0 <= X(vIndices_2) .* X(vIndicesP1)); %
        I2_NEG = find(0 > X(vIndices_2) .* X(vIndicesP1)); % setxor(I2,[1:length(vIndicesP1)]);  
        
        %% calculate const parts c_p^j f. j = k-1,k,k+1  - as matrix C_0,
        %% C_n+1 has to be handled separately

        vI1Neg  = vIndices_1(I1_NEG);
        vI2Neg  = vIndices_2(I2_NEG);
        vIM1Neg = vIndicesM1(I1_NEG);
        vIP1Neg = vIndicesP1(I2_NEG);
        vI1     = vIndices_1(I1);
        vI2     = vIndices_2(I2);
        vLamI1Neg       = vLam(vI1Neg);
        vLamI1NegSquare = vLamI1Neg.^2;
        vLamIP1Neg      = vLam(vIP1Neg);
        vLamIP1NegSquare = vLamIP1Neg.^2;
        vLamIM1Neg      = vLam(vIM1Neg);
        vLamIM1NegSquare = vLamIM1Neg.^2;
        vLamI2Neg       = vLam(vI2Neg);
        vLamI2NegSquare = vLamI2Neg.^2;
        
        A(vIndicesM1(I1))        = vLam(vI1) .* absX(vI1);
        A(vIM1Neg)               = vLamI1NegSquare .* XSquare(vI1Neg) ...
                                   ./ ...
                                   (vLamI1Neg .* absX(vI1Neg) + vLamIM1Neg .* absX(vIM1Neg));
  
        B(vIndicesP1(I2))        = vLam(vI2) .* absX(vI2);
        B(vIP1Neg)    = vLamI2NegSquare .* XSquare(vI2Neg) ...
                                   ./ ...
                                   (vLamI2Neg .* absX(vI2Neg) + vLamIP1Neg .* absX(vIP1Neg));
        % for p == 2; cClass.p *
        C_j(1)                     = 2.0 * (vLam(1) .* absX(1)     + B(2));                               %.^(cClass.p-1);        
        C_j(vIndicesWithOutBorder) = 2.0 * (A(vIndicesWithOutBorder-1) + B(vIndicesWithOutBorder+1));     %.^(cClass.p-1);        
        C_j(iNrDims)               = 2.0 * (A(iNrDims-1)               + vLam(iNrDims) .* absX(iNrDims)); %.^(cClass.p-1);        
        
        %% calculate collected formula terms for vk-1*vk and vk+1 * vk -
        %% valid f indices 2-N (index k = 1 has to be handled separately)
        T1(vI1)     = zeros(1,length(I1)) + C_j(vI1) .* absX(vI1);
        T1(vI1Neg) = ( vLamI1NegSquare .* C_j(vI1Neg) .* XSquare(vI1Neg) .* absX(vI1Neg) ...
                                 - ...
                                 vLamIM1NegSquare .* C_j(vIM1Neg) .* absX(vI1Neg) .* XSquare(vIM1Neg) ...
                                 + ...
                                 2 * vLamI1Neg .* C_j(vI1Neg) .* XSquare(vI1Neg) .* absX(vIM1Neg) .* vLamIM1Neg) ...
                                 / ...
                                 (vLamI1Neg .* absX(vI1Neg) + absX(vIM1Neg) .* vLamIM1Neg).^2;

        %% valid f indices 1-N-1 (index k = N has to be handled separately)
        T2(vI2)     = zeros(1,length(I2)) + C_j(vI2) .* absX(vI2);
        T2(vI2Neg) = ( vLamI2NegSquare .* C_j(vI2Neg) .* XSquare(vI2Neg) .* absX(vI2Neg) ...
                                 - ...
                                 vLamIP1NegSquare .* C_j(vIP1Neg) .* absX(vI2Neg) .* XSquare(vIP1Neg) ...
                                 + ...
                                 2 * vLamI2Neg .* C_j(vI2Neg) .* XSquare(vI2Neg) .* absX(vIP1Neg) .* vLamIP1Neg) ...
                                 / ...
                                 (vLamI2Neg .* absX(vI2Neg) + absX(vIP1Neg) .* vLamIP1Neg).^2;                               

        %% all together
        vValue(i,1)                     = C * vCp(i) * (C_j(1) .* absX(1) + T2(1)); % first dimension is border
        vValue(i,vIndicesWithOutBorder) = C * vCp(i) * (T1(vIndicesWithOutBorder) + T2(vIndicesWithOutBorder)); % index-1 due to internal shift
        vValue(i,iNrDims)               = C * vCp(i) * (T1(iNrDims) + C_j(iNrDims) .* absX(iNrDims) ); % first dimension is border index-1 due to internal shift
    end
return
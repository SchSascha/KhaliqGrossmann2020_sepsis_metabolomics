% Distance method derived to x for generic metric functions
% for the lambda scaling we consider only a global lambda
function vValue = distance_deriv_x(cClass, vX, vY, varargin) 
    if(isvector(vX)) % X are prototypes, Y are datapoints
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

    X_ALL = vX-vY;
    vValue = zeros(iNrVectors,iNrDims);
    C = cClass.tau^2/2;
    vIndicesWithOutBorder = [2:iNrDims-1];   
    vIndices_1  = [2:iNrDims];
    vIndicesM1  = vIndices_1-1;
    vIndices_2  = [1:iNrDims-1];
    vIndicesP1  = vIndices_2+1;
    U1          = zeros(1,iNrDims-1);    
    V1          = zeros(1,iNrDims-1);    
    U2          = zeros(1,iNrDims);        
    V2          = zeros(1,iNrDims);     
    vLam = vLam(1,:);    
    for i=1:iNrVectors
        X    = X_ALL(i,:);
        absX = abs(X);
        I1     = find(0 <= X(vIndices_1) .* X(vIndicesM1)); %
        I1_NEG = find(0 > X(vIndices_1) .* X(vIndicesM1)); % xor is slow setxor(I1,[1:length(vIndicesM1)]);

        I2     = find(0 <= X(vIndices_2) .* X(vIndicesP1)); %
        I2_NEG = find(0 > X(vIndices_2) .* X(vIndicesP1)); % setxor(I2,[1:length(vIndicesP1)]);

        U1(vIndicesM1(I1))     = zeros(1,length(I1));
        U1(vIndicesM1(I1_NEG)) = (vLam(vIndicesM1(I1_NEG)) .* X(vIndicesM1(I1_NEG))    ./  (vLam(vIndices_1(I1_NEG)) .* absX(vIndices_1(I1_NEG)) + vLam(vIndicesM1(I1_NEG)) .* absX(vIndicesM1(I1_NEG)))).^2;

        V1(vIndicesM1(I1))     = vLam(I1);
        V1(vIndicesM1(I1_NEG)) = vLam(vIndices_1(I1_NEG)) .* abs(X(vIndices_1(I1_NEG))) ./ (vLam(vIndices_1(I1_NEG)) .* absX(vIndices_1(I1_NEG)) + vLam(vIndicesM1(I1_NEG)) .*  absX(vIndicesM1(I1_NEG)));

        U2(vIndicesP1(I2))      = zeros(1,length(I2));
        U2(vIndicesP1(I2_NEG))  = (vLam(vIndicesP1(I2_NEG)) .* X(vIndicesP1(I2_NEG))   ./  (vLam(vIndices_2(I2_NEG)) .* absX(vIndices_2(I2_NEG)) + vLam(vIndicesP1(I2_NEG)) .* absX(vIndicesP1(I2_NEG)))).^2;
    
        V2(vIndicesP1(I2))      = vLam(I2);
        V2(vIndicesP1(I2_NEG))  = vLam(vIndices_2(I2_NEG)) .* abs(X(vIndices_2(I2_NEG))) ./(vLam(vIndices_2(I2_NEG)) .* absX(vIndices_2(I2_NEG)) + vLam(vIndicesP1(I2_NEG)) .* absX(vIndicesP1(I2_NEG)));

        vValue(i,1)                     = C * ( 2- 0  - U2(2)) .* (1 + V2(2)) .* X(1) ; % first dimension is border
        vValue(i,vIndicesWithOutBorder) = C * ( 2- U1(vIndicesWithOutBorder-1) - U2(vIndicesWithOutBorder+1)) .* (V1(vIndicesWithOutBorder-1) + V2(vIndicesWithOutBorder+1)) .* X(vIndicesWithOutBorder); % index-1 due to internal shift
        vValue(i,iNrDims)               = C * ( 2- U1(iNrDims-1) - 0) .* (V1(iNrDims-1) + 1) .* X(iNrDims); % first dimension is border index-1 due to internal shift
    end    
return
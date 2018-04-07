% Distance method for generic metric functions
function vValue = distance(cClass, vX, vY, varargin)
	iNrDims = length(vX);
    %vValue  = zeros(1,iNrDims);
    if(nargin<4)        
        vLam = cClass.vLam;
    else
        vLam = varargin{1};
    end
    if ~(length(vLam) == iNrDims || size(vLam,2) == iNrDims)
        vLam = ones(1,iNrDims);        
    end
    
    X = vX - vY;
    C = cClass.tau^2/2;
    vIndicesWithOutBorder = [2:iNrDims-1];
    vIndices_1  = [2:iNrDims];
    vIndicesM1  = vIndices_1-1;
    vIndices_2  = [1:iNrDims-1];
    vIndicesP1  = vIndices_2+1;
    A           = zeros(1,iNrDims-1);    
    B           = zeros(1,iNrDims);       
    
    if(sum(vLam) ~= iNrDims) % no true lambda update
        vLamSquare  = vLam.^2;

        absX    = abs(X);
        XSquare = X.^2;
        I1     = (0 <= X(vIndices_1) .* X(vIndicesM1)); %
        I1_NEG = ~I1;

        I2     = (0 <= X(vIndices_2) .* X(vIndicesP1)); %
        I2_NEG = ~I2; % -- > slow setxor(I2,vVecP1);

        vIndI1NEG = vIndices_1(I1_NEG);
        vIndI2NEG = vIndices_2(I2_NEG);

        A(vIndicesM1(I1))        = vLam(vIndices_1(I1)) .* absX(vIndices_1(I1));
        A(vIndicesM1(I1_NEG))    = vLamSquare(vIndI1NEG) .* XSquare(vIndI1NEG) ...
                                   ./ ... 
                                   (vLam(vIndI1NEG) .* absX(vIndI1NEG) + vLam(vIndicesM1(I1_NEG)) .* absX(vIndicesM1(I1_NEG)));


        B(vIndicesP1(I2))        = vLam(vIndices_2(I2)) .* absX(vIndices_2(I2));
        B(vIndicesP1(I2_NEG))    = vLamSquare(vIndI2NEG) .* XSquare(vIndI2NEG) ...
                                    ./ ...
                                   (vLam(vIndI2NEG) .* absX(vIndI2NEG) + vLam(vIndicesP1(I2_NEG)) .* absX(vIndicesP1(I2_NEG)));

        vValue = (C .* [vLam(1) .* absX(1) + B(2), A(vIndicesWithOutBorder-1) + B(vIndicesWithOutBorder+1),A(iNrDims-1) + vLam(iNrDims) .* absX(iNrDims)]).^2 ;
        %vValue(1)                      = (C * (vLam(1) .* absX(1) + B(2))).^2;                              
        %vValue(vIndicesWithOutBorder)  = (C * (A(vIndicesWithOutBorder-1) + B(vIndicesWithOutBorder+1))).^2;
        %vValue(iNrDims)                = (C * (A(iNrDims-1) + vLam(iNrDims) .* absX(iNrDims))).^2;          
    else
        absX    = abs(X);
        XSquare = X.^2;
        I1     = (0 <= X(vIndices_1) .* X(vIndicesM1)); %
        I1_NEG = ~I1;

        I2     = (0 <= X(vIndices_2) .* X(vIndicesP1)); %
        I2_NEG = ~I2; % -- > slow setxor(I2,vVecP1);

        vIndI1NEG = vIndices_1(I1_NEG);
        vIndI2NEG = vIndices_2(I2_NEG);
        vIM1Neg   = vIndicesM1(I1_NEG);
        vIP1Neg   = vIndicesP1(I2_NEG);

        A(vIndicesM1(I1))        = absX(vIndices_1(I1));
        A(vIM1Neg)               = XSquare(vIndI1NEG) ./ (absX(vIndI1NEG) + absX(vIM1Neg));
        B(vIndicesP1(I2))        = absX(vIndices_2(I2));
        B(vIP1Neg)               = XSquare(vIndI2NEG) ./ (absX(vIndI2NEG) + absX(vIP1Neg));

        vValue = (C .* [absX(1) + B(2), A(vIndicesWithOutBorder-1) + B(vIndicesWithOutBorder+1),A(iNrDims-1) + absX(iNrDims)]).^2 ;        
    end

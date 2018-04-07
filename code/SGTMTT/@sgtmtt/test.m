function [TestLLH, Resp, LLhood,LhoodV]=test(g, X)

error(nargchk(1, 2, nargin, 'struct'));

if nargin==1
    X = g.Data.X{1};
else
    if(1)
        if(g.Data.iMethod == 1) 
            X = bsxfun(@times,bsxfun(@minus,X,g.Data.NormPar1),1.0./g.Data.NormPar2);
        else
            'no such normalisation defined - fix code first'
            return;
        end
    end
end

[N D] = size(X);
if D ~= g.Data.D
    error('Dimensions of the training and test dataset are mismatched.');
end
Y = g.Par.Y;
[DistXY] = g.Par.oMetric.distance(Y, X); %include relevances here
BjnMat = (g.Par.Beta/(2*pi))^(g.Par.Deff/2)*exp(-0.5*g.Par.Beta*DistXY);
Bjn = mat2cell(BjnMat,g.K,N);
[Resp1, ~, ~, LLhood, TestLLH, LhoodV]=forwback(g.Par.PI, g.Par.A, Bjn, N);
if nargout>=2, Resp = Resp1; end
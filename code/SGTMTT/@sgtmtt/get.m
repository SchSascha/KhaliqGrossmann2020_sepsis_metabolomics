function val = get(g, propName)
% Get sgtmtt properties from the specified object
% and return the value

switch lower(propName)
    case 'k'
        val = g.K;
    case 'n'
        val = g.Data.N;
    case 'd'
        val = g.Data.D;
    case 'l'
        val = g.L;
    case 'par'
        val = g.Par;
    case 'pi'
        val = g.Par.PI;
    case 'a'
        val = g.Par.A;
    case 'y'
        val = g.Par.Y;
    case 'w'
        val = g.Par.W;
    case 'beta'
        val = g.Par.Beta;
    case 'x'
        val = g.Data.X;
    case 'llhood'
        val = g.LLhood;
    case 'llhoodxobs'
        val = g.LLhoodxObs;
    case 'r'
        val = g.RespMat;
    case 'u'
        val = g.U;
    case 'm'
        val = g.M;
    case 's'
        val = g.s;
    case 'fi'    
        val = g.FI;
    case 'maxiter'
        val = g.MaxIter;
    case 'distxy'
        val = g.DistXY;
    case 'status'
        val = g.Status;
    case 'verb'
        val = g.Verb;
    case 'ometric'
        val = g.Par.oMetric;        
    case 'tolerance'
        val = g.Tolerance;        
    case 'labels'
        val = g.Labels;     
    otherwise
        error('Unknown property');
end
function g=set(g,varargin)
%Set SGTMTT properties

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch lower(prop)
%         case 'k'
%             g.K = val;
        case 'data'
            g.Data = val;
        case 'l'
            g.L = val;
        case 'm'
            g.M = val;
        case 'u'
            g.U = val;            
            g.K = size(g.U,1);
        case 'mu'
            g.MU = val;            
        case 'ometric'
            g.Par.oMetric = val;
        case 's'
            g.s = val;
        case 'maxiter'
            g.MaxIter = val;
            %         case 'ModelPar'
            %             g.Par = val;
            %         case 'W'
            %             g.Par.W = val;
        case 'pi'
            g.Par.PI = val;
        case 'a'
            g.Par.A = val;
        case 'beta'
            g.Par.Beta = val;
        case 'deff'
            g.Par.Deff = val;
        case 'tolerance'
            g.Tolerance = val;
        case 'initmeth'
            g.InitMeth = val;
        case 'status'
            g.Status = val;
        case 'verbose'
            if islogical(val)
                g.Verb = val;
            else
                error('Only TRUE or FALSE');
            end
        otherwise
            error('Unknown property');
    end
end


% Lee distance class - based on metric class
function cLeeMetric = metric_lee(varargin)
    % Lee metric constructor
    switch nargin
        case 0
            %if no inputargs are given create default
            cMetric      = metric('Lee metric', 'metric_lee');
            cLeeMetric.p    = 2;
            cLeeMetric.tau  = 1;
            cLeeMetric.vLam = [];
            cLeeMetric.bWithAdaptation = 1;            
            cLeeMetric   = class(cLeeMetric,'metric_lee',cMetric);
        case 1
           % if single argument of class euclidean, return it
           if (isa(varargin{1},'metric_euclidean'))
              s = varargin{1}; 
           else
              error('Input argument is not an euclidean object')
           end
        case 3
           % create object using specified values
           cMetric = metric(varargin);
           cLeeMetric.p = 2;           
           cLeeMetric.tau = 1;           
           cLeeMetric.vLam = [];          
           cLeeMetric.bWithAdaptation = 1;           
           cLeeMetric   = class(cLeeMetric,'metric_lee',cMetric);
        otherwise
           error('Wrong number of input arguments')
        end            
return
% Basis class for generic metric functions
function cMetric = metric(varargin) 
    % Metric constructor function for metric object
    % cMetric = metric(parameters) 
switch nargin
    case 0
        % if no input arguments, create a default object 
        %- with default members
        cMetric.descriptor = 'none'
        cMetric.type       = 'none'
        cMetric = class(cMetric, 'metric');
    case 1
        %if single argument of class metric, return it
        if (isa(varargin{1},'metric'))
            cMetric = varargin{1};
        else
            error('Wrong argument type')
        end 
    case 2
        % create object with specified values
        cMetric.descriptor = varargin{1};
        cMetric.type       = varargin{2};        
        cMetric = class(cMetric, 'metric');
otherwise 
        error('Wrong number of input arguments')
end
return
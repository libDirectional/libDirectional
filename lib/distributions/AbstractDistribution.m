classdef (Abstract) AbstractDistribution
    % Abstract base class for all distributions
    
    properties
        dim {mustBeGreaterThanOrEqual(dim,1)}
    end
    
    methods (Abstract)
        % Evaluate pdf at positions stored in xa
        pdf(this, xa);
    end
    
end
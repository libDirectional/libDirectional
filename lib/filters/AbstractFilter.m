classdef (Abstract) AbstractFilter < handle & matlab.mixin.Copyable
    % Abstract base class for all filters
 
    methods (Abstract)
        setState(this, state)
        est = getEstimate(this)
        est = getPointEstimate(this)
    end
    
    methods
        function d = dim(this)
            % Covenience function to get dimension of the filter. Overwrite
            % if filter is not directly based on a distribution
            d = this.getEstimate().dim;
        end
    end
end


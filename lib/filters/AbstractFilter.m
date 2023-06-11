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
            arguments
                this (1,1) AbstractFilter
            end
            d = this.getEstimate().dim;
        end
        function plotFilterState(this)
            arguments
                this (1,1) AbstractFilter
            end
            this.getEstimate().plot();
        end
    end
end


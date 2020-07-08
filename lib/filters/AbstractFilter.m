classdef (Abstract) AbstractFilter < handle & matlab.mixin.Copyable
    % Abstract base class for all filters
    properties
        % See conventions of underlying manifold in
        % AbstractHypertoroidalDistribution or AbstractHypertoroidalDistribution.
        % For Axial filter see AbstractAxialFilter.
        dim
    end
    
    methods (Abstract)
        setState(this, state)
        est = getEstimate(this)
    end
    methods
        function mean=getEstimateMean(this)
            dist = this.getEstimate;
            mean = dist.meanDirection;
        end
    end
end


classdef (Abstract) AbstractHypertoroidalFilter < AbstractFilter
    % Abstract base class for filters on the hypertorus
    
    properties
    end
    
	methods (Abstract)
        setState(this, state)
        est = getEstimate(this)
    end
    methods
        function mean=getEstimateMean(this)
            dist=this.getEstimate;
            mean=dist.circularMean;
        end
    end
end


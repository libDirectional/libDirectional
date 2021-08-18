classdef (Abstract) AbstractHypertoroidalFilter < AbstractFilter
    % Abstract base class for filters on the hypertorus
    methods
        function est = getPointEstimate(this)
            arguments
                this AbstractHypertoroidalFilter
            end
            est = this.getEstimate().meanDirection();
        end
        function mean = getEstimateMean(this)
            arguments
                this AbstractHypertoroidalFilter
            end
            mean = this.getEstimate().meanDirection();
        end
    end
end
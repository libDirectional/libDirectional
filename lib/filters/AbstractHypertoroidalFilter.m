classdef (Abstract) AbstractHypertoroidalFilter < AbstractFilter
    % Abstract base class for filters on the hypertorus
    methods
        function est = getPointEstimate(this)
            arguments
                this (1,1) AbstractHypertoroidalFilter
            end
            est = this.getEstimate().meanDirection();
        end
        function mean = getEstimateMean(this)
            % Just an alternative interface. Fall back to getPointEstimate,
            % which should be overwritten by inherting functions.
            arguments
                this (1,1) AbstractHypertoroidalFilter
            end
            mean = this.getPointEstimate();
        end
    end
end
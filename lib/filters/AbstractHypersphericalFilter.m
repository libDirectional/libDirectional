classdef (Abstract) AbstractHypersphericalFilter < AbstractFilter
    % Abstract base class for filters on the hypersphere (without antipodal
    % symmetry). Use AbstractAxialFilter or AbstractHyperhemisphericalDistribution
    % for antipodally symmetric problems!
    methods
        function est = getPointEstimate(this)
            arguments
                this (1,1) AbstractHypersphericalFilter
            end
            est = this.getEstimate().meanDirection;
        end
        function mean = getEstimateMean(this)
            % Just an alternative interface. Fall back to getPointEstimate,
            % which should be overwritten by inherting functions.
            arguments
                this (1,1) AbstractHypersphericalFilter
            end
            mean = this.getPointEstimate();
        end
    end
end
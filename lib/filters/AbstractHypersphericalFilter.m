classdef (Abstract) AbstractHypersphericalFilter < AbstractFilter
    % Abstract base class for filters on the hypersphere (without antipodal
    % symmetry). Use AbstractAxialFilter or AbstractHyperhemisphericalDistribution
    % for antipodally symmetric problems!
    methods
        function est = getPointEstimate(this)
            arguments
                this AbstractHypersphericalFilter
            end
            est = this.getEstimate().meanDirection;
        end
        function mean = getEstimateMean(this)
            arguments
                this AbstractHypersphericalFilter
            end
            mean = this.getEstimate().meanDirection;
        end
    end
end
classdef (Abstract) AbstractHyperhemisphericalFilter < AbstractFilter
    methods
        function est = getPointEstimate(this)
            arguments
                this AbstractHyperhemisphericalFilter
            end
            est = this.getEstimate().mode;
        end
    end
end
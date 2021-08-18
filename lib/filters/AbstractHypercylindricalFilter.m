classdef (Abstract) AbstractHypercylindricalFilter < AbstractLinPeriodicFilter
    methods
        function mu = getPointEstimate(this)
            mu = this.getEstimate().hybridMean();
        end
    end
end


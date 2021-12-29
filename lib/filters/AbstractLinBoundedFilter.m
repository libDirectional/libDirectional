classdef (Abstract) AbstractLinBoundedFilter < AbstractFilter
    % Class for Cartesian products of bounded (can be periodic or non-periodic) and linear domains
    methods
        function est = getPointEstimate(this)
            arguments
                this (1,1) AbstractLinBoundedFilter
            end
            est = this.getEstimate().hybridMean;
        end
    end
end


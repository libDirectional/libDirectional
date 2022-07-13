classdef (Abstract) AbstractLinBoundedDistribution < AbstractCartProdDistribution
    % For Cartesian products of linear and bounded (periodic or parts of
    % Eculidean spaces) domains. Assumption is that it is bounded x R^n (in
    % this order)
    properties (SetAccess = protected)
        linD {mustBeInteger,mustBeNonnegative} % number of linear dimensions
        boundD {mustBeInteger,mustBeNonnegative} % number of bounded (e.g., periodic or hyperrectangular) dimensions
    end
    methods
        function m = mean(this)
            m = this.hybridMean();
        end
        function m = hybridMean(this)
            m = [this.marginalizeLinear().mean(); this.marginalizePeriodic().mean()];
        end
    end
end
classdef (Abstract) AbstractLinBoundedDistribution < AbstractCartProdDistribution
    % For Cartesian products of linear and bounded (periodic or parts of
    % Eculidean spaces) domains. Assumption is that it is bounded x R^n (in
    % this order)
    properties (SetAccess = protected)
        linD {mustBeInteger,mustBeNonnegative} % number of linear dimensions
        boundD {mustBeInteger,mustBeNonnegative} % number of bounded (e.g., periodic or hyperrectangular) dimensions
        periodicManifoldType char ...
            {mustBeMember(periodicManifoldType,{'unspecified','hypertorus','hypersphere','hyperhemisphere'})} = 'unspecified'
    end
end
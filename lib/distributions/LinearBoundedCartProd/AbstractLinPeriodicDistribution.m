classdef (Abstract) AbstractLinPeriodicDistribution < AbstractLinBoundedDistribution
    % For Cartesian products of linear and periodic domains. Assumption is
    % that it is bounded x R^n (in this order)
    properties (SetAccess = protected)
        periodicManifoldType char ...
            {mustBeMember(periodicManifoldType,{'unspecified','hypertorus','hypersphere','hyperhemisphere'})} = 'unspecified'
    end
end
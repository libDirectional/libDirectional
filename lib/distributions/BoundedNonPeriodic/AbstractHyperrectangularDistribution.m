classdef (Abstract) AbstractHyperrectangularDistribution < AbstractDistribution
    properties
        bounds (:,2)
    end
    methods
        function s = getManifoldSize(this)
            s = prod(diff(this.bounds,1,2));
        end
    end
end
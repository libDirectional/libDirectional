classdef (Abstract) AbstractNonConditionalDistribution < AbstractDistribution
    methods (Abstract)
        getManifoldSize(this);
    end
end


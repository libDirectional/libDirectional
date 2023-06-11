classdef (Abstract) AbstractPeriodicDistribution < AbstractNonConditionalDistribution
    methods
        function m = mean(this)
            m = this.meanDirection();
        end
    end
end
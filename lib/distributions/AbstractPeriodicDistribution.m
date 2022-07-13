classdef (Abstract) AbstractPeriodicDistribution < AbstractDistribution
    methods
        function m = mean(this)
            m = this.meanDirection();
        end
    end
end
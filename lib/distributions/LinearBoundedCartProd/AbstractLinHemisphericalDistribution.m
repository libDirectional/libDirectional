classdef (Abstract) AbstractLinHemisphericalDistribution < AbstractLinPeriodicDistribution
 
    methods
        function this = AbstractLinHemisphericalDistribution()
            this.periodicManifoldType = 'hyperhemisphere';
        end
    end
end


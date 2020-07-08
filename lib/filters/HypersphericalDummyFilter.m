classdef HypersphericalDummyFilter < AbstractDummyFilter & AbstractHypersphericalFilter
    methods
        function this = HypersphericalDummyFilter(dim)
            assert(dim>=2,'Dim must be the dimension of the Euclidean space that the densities are embedded in. Use at least 2.');
            this.dim=dim;
            % Do nothing
        end
        
        function hfd = getEstimate(this)
            hfd=HypersphericalUniformDistribution(this.dim);
        end
        
        function mean=getEstimateMean(this)
            mean=HypersphericalUniformDistribution(this.dim).sample(1);
        end
        
    end
end
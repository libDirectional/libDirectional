classdef HypersphericalDummyFilter < AbstractDummyFilter & AbstractHypersphericalFilter & AbstractDummyFilter
    methods
        function this = HypersphericalDummyFilter(dim)
            arguments
                dim (1,1) double {mustBeGreaterThanOrEqual(dim,2)} % 'Dim must be the dimension of the Euclidean space that the densities are embedded in. Use at least 2.'
            end
            this=this@AbstractDummyFilter(dim);
        end
        
        function hfd = getEstimate(this)
            hfd=HypersphericalUniformDistribution(this.dim);
        end
        
        function mean=getEstimateMean(this)
            mean=HypersphericalUniformDistribution(this.dim).sample(1);
        end
        
    end
end
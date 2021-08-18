classdef HypersphericalDummyFilter < AbstractDummyFilter & AbstractHypersphericalFilter
    methods
        function this = HypersphericalDummyFilter(dim)
            arguments
                dim (1,1) double {mustBeGreaterThanOrEqual(dim,2)} % 'Dim must be the dimension of the Euclidean space that the densities are embedded in. Use at least 2.'
            end
            this.dist = HypersphericalUniformDistribution(dim);
        end
        
        function est = getPointEstimate(this)
            est = getPointEstimate@AbstractDummyFilter(this);
        end 
    end
end
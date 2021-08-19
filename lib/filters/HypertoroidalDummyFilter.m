classdef HypertoroidalDummyFilter < AbstractDummyFilter & AbstractHypertoroidalFilter
    
    methods
        function this = HypertoroidalDummyFilter(dim)
            arguments
                dim (1,1) double {mustBeGreaterThanOrEqual(dim,1)}
            end
            this.dist = HypertoroidalUniformDistribution(dim);
        end
        
        function est = getPointEstimate(this)
            est = getPointEstimate@AbstractDummyFilter(this);
        end   
    end
end
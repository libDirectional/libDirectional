classdef HypertoroidalDummyFilter < AbstractDummyFilter & AbstractHypertoroidalFilter
    
    methods
        function this = HypertoroidalDummyFilter(dim)
            this.dim=dim;
        end
        
        function hfd = getEstimate(this)
            hfd=HypertoroidalUniformDistribution(this.dim);
        end
        
        function mean=getEstimateMean(this)
            mean=HypertoroidalUniformDistribution(this.dim).sample(1);
        end
        
    end
end
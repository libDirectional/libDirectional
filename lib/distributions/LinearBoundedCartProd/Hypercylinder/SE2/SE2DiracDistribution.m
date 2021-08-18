classdef SE2DiracDistribution < HypercylindricalDiracDistribution & AbstractSE2Distribution

    methods
        function this = SE2DiracDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (3 x L)
            %       Dirac locations - first bounded then linear dimensions
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                d_ (3,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2); % All Diracs have equal weights by default
            end
            this@HypercylindricalDiracDistribution(1, d_, w_);
        end
    end
    
     methods (Static)
        function ddist = fromDistribution(dist, nParticles)
            arguments
                dist (1,1) AbstractHypercylindricalDistribution
                nParticles (1,1) {mustBeInteger,mustBePositive}
            end
            assert(dist.boundD==1,dist.linD==2);
            ddist = SE2DiracDistribution(dist.sample(nParticles), 1/nParticles * ones(1,nParticles));
        end
     end
    
end


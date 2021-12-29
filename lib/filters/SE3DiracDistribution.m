classdef SE3DiracDistribution < LinBoundedDiracDistribution & AbstractSE3Distribution
    methods
        function this = SE3DiracDistribution(d_, w_)
            arguments
                d_ (7,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (7 x L)
            %       Dirac locations (S^3 x R^2)
            %   w_ (1 x L)
            %       weights for each Dirac
            assert ( max(abs(vecnorm(d_(1:4,:))-1)) < 1E-5);
            this@LinBoundedDiracDistribution(4, d_, w_);
        end
        function dist = marginalizeLinear(this)
            arguments
                this (1,1) SE3DiracDistribution
            end
            dist = HyperhemisphericalDiracDistribution(this.d(1:4,:), this.w);
        end
        function m = mean(this)
            arguments
                this (1,1) SE3DiracDistribution
            end
            m = this.hybridMean();
        end
    end    
    
    methods (Static)
        function ddist = fromDistribution(dist, nParticles)
            arguments
                dist (1,1) AbstractSE3Distribution
                nParticles (1,1) {mustBeInteger,mustBePositive}
            end
            ddist = SE3DiracDistribution(dist.sample(nParticles), 1/nParticles * ones(1,nParticles));
        end
    end
end
classdef SE2ParticleFilter < HypercylindricalParticleFilter & AbstractSE2Filter
    % SIR particle filter for the SE2
    methods
        function this = SE2ParticleFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            arguments
                nParticles (1,1) double {mustBeInteger, mustBePositive}
            end
            this@HypercylindricalParticleFilter(nParticles, 1, 2);
        end
        
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractHypertoroidalDistribution)
            %       new state
            arguments
                this (1,1) HypercylindricalParticleFilter
                dist_ (1,1) AbstractHypercylindricalDistribution
            end
            if isa(dist_,'HypercylindricalDiracDistribution')
                dist_ = SE2PWDDistribution(dist_.d, dist_.w);
            elseif ~isa(dist_, 'SE2PWDDistribution')
                dist_ = SE2PWDDistribution.fromDistribution(dist_, numel(this.dist.w));
            end
            this.dist = dist_;
        end        
    end
    
end


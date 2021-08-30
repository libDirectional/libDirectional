classdef ToroidalParticleFilter < AbstractToroidalFilter & HypertoroidalParticleFilter
    % SIR Particle filter on the torus    
    
    methods
        function this = ToroidalParticleFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use
            arguments
                nParticles (1,1) double {mustBeInteger, mustBePositive}
            end
            this@HypertoroidalParticleFilter(nParticles,2);
            this.dist = ToroidalWDDistribution(this.dist.d,this.dist.w);
        end
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractToroidalDistribution)
            %       new state
            arguments
                this (1,1) ToroidalParticleFilter
                dist_ (1,1) AbstractToroidalDistribution
            end
            if ~isa(dist_, 'ToroidalWDDistribution')
                dist_ = ToroidalWDDistribution(dist_.sample(length(this.dist.d)));
            end
            this.dist = dist_;
        end
    end
    
end


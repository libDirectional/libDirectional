classdef ToroidalParticleFilter < AbstractToroidalFilter & HypertoroidalParticleFilter
    % SIR Particle filter on the torus    
    
    methods
        function this = ToroidalParticleFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use  
            this@HypertoroidalParticleFilter(nParticles,2);
            this.wd = ToroidalWDDistribution(this.wd.d,this.wd.w);
        end
        function setState(this, wd_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractToroidalDistribution)
            %       new state            
            assert (isa (wd_, 'AbstractToroidalDistribution'));
            if ~isa(wd_, 'ToroidalWDDistribution')
                wd_ = ToroidalWDDistribution(wd_.sample(length(this.wd.d)));
            end
            this.wd = wd_;
        end
    end
    
end


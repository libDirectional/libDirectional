classdef CircularParticleFilter < AbstractCircularFilter & HypertoroidalParticleFilter
    % A sequential importance resampling (SIR) particle filter on the circle 
    % based on the wrapped Dirac distribution.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Bayesian Filtering in Circular State Spaces
    % arXiv preprint: Systems and Control (cs.SY), January 2015.
    
    methods
        function this = CircularParticleFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            arguments
                nParticles (1,1) double {mustBeInteger,mustBePositive}
            end
            this@HypertoroidalParticleFilter(nParticles, 1);
            dist_ = WDDistribution((0:nParticles-1)/nParticles*2*pi);
            this.setState(dist_);
        end
        
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractCircularDistribution)
            %       new state         
            arguments
                this (1,1) CircularParticleFilter
                dist_ (1,1) AbstractCircularDistribution
            end
            if ~isa(dist_, 'WDDistribution')
                dist_ = WDDistribution(dist_.sample(length(this.dist.d)));
            end
            this.dist = dist_;
        end
        
        function likelihoodVal = associationLikelihood(this, likelihood)
            likelihoodVal = associationLikelihood@AbstractParticleFilter(this, likelihood);
        end
    end
    
end


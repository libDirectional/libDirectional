classdef CircularParticleFilter < AbstractCircularFilter & HypertoroidalParticleFilter
    % A sequential importance resampling (SIR) particle filter on the circle 
    % based on the wrapped Dirac distribution.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Bayesian Filtering in Circular State Spaces
    % arXiv preprint: Systems and Control (cs.SY), January 2015.
        
    properties
    end
    
    methods
        function this = CircularParticleFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use         
            assert(isscalar(nParticles));
            assert(nParticles >= 1);
            this@HypertoroidalParticleFilter(nParticles, 1);
            wd_ = WDDistribution((0:nParticles-1)/nParticles*2*pi);
            this.setState(wd_)
        end
        
        function setState(this, wd_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractCircularDistribution)
            %       new state            
            assert (isa (wd_, 'AbstractCircularDistribution'));
            if ~isa(wd_, 'WDDistribution')
                wd_ = WDDistribution(wd_.sample(length(this.wd.d)));
            end
            this.wd = wd_;
        end
        
        function predictIdentity(this, noiseDistribution)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            assert (isa (noiseDistribution, 'AbstractCircularDistribution'));
            this.predictNonlinear(@(x) x, noiseDistribution);
        end
        
        function predictNonlinear(this, f, noiseDistribution)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi)
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            
            assert(isa (noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(f,'function_handle'));
            %apply f
            wdF = this.wd.applyFunction(f);           
            %calculate effect of (additive) noise 
            n = length(this.wd.d);
            noise = noiseDistribution.sample(n);
            for i=1:n
                wdF.d(i) = mod(wdF.d(i) + noise(i),2*pi);
            end
            this.wd = wdF;
        end
        
        function predictNonlinearNonAdditive(this, f, samples, weights)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k), w(k))    mod 2pi,
            % where w(k) is non-additive noise given by samples and weights.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) x W to [0,2pi) (W is the space
            %       containing the noise samples)
            %   noiseSamples (d x n matrix)
            %       n samples of the noise as d-dimensional vectors
            %   noiseWeights (1 x n vector)
            %       weight of each sample
            
            %(samples, weights) are discrete approximation of noise
            assert(size(weights,1) == 1, 'weights most be row vector')
            assert(size(samples,2) == size(weights,2), 'samples and weights must match in size');
            assert(isa(f,'function_handle'));
            
            weights = weights/sum(weights); %ensure normalization
            n = length(this.wd.d);
            noiseIds = discretesample(weights, n);
            d = zeros(1, n);
            
            for i=1:n
                d(i) = f(this.wd.d(i),samples(noiseIds(i)));
            end
            this.wd = WDDistribution(d,this.wd.w);
        end
        
        function updateIdentity(this, noiseDistribution, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            assert(isa(noiseDistribution, 'AbstractCircularDistribution'));
            assert(isscalar(z));
            this.updateNonlinear(LikelihoodFactory.additiveNoiseLikelihood(@(x) x, noiseDistribution), z);
        end
        
        function updateNonlinear(this, likelihood, z)
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi) to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            
            % You can either use a likelihood depending on z and x
            % and specify the measurement as z or use a likelihood that
            % depends only on x and omit z.
            if nargin==2
                this.wd = this.wd.reweigh(likelihood);
            else
                this.wd = this.wd.reweigh(@(x) likelihood(z,x));
            end
            this.wd = WDDistribution(this.wd.sample(length(this.wd.d))); %use SIR.
        end
        
        function likelihoodVal=associationLikelihood(this,likelihood)
            likelihoodVal=sum(likelihood.pdf(this.getEstimate.d).*this.getEstimate.w);
        end
    end
    
end


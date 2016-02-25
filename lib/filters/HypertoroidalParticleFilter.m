classdef HypertoroidalParticleFilter < AbstractToroidalFilter
    % SIR particle filter on the hypertorus
    
    properties
        wd
    end
    
    methods
        function this = HypertoroidalParticleFilter(nParticles,dim)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            this.wd = HypertoroidalWDDistribution(repmat((0:nParticles-1)/nParticles*2*pi,[dim,1]));
        end
        
        function setState(this, wd_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractHypertoroidalDistribution)
            %       new state            
            assert (isa (wd_, 'AbstractHypertoroidalDistribution'));
            if ~isa(wd_, 'HypertoroidalWDDistribution')
                wd_ = HypertoroidalWDDistribution(wd_.sample(length(this.wd.d)));
            end
            this.wd = wd_;
        end
        
        function predictIdentity(this, noiseDistribution)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            this.predictNonlinear(@(x) x, noiseDistribution);
        end
        
        function predictNonlinear(this, f, noiseDistribution)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   f (function handle)
            %       system function
            %   noiseDistribution (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            
            assert (isa (noiseDistribution, 'AbstractHypertoroidalDistribution'));
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
            %       function from [0,2pi)^dim x W to [0,2pi)^dim (W is the space
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
            this.wd.d=d;
        end
        
        function updateIdentity(this, noiseDistribution, z)
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
            %       function from Z x [0,2pi)^dim to [0, infinity), where Z is
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
            this.wd.d = this.wd.sample(length(this.wd.d));
            this.wd.w = 1/size(this.wd.d,2)*ones(1,size(this.wd.d,2));
        end
        
        function wd = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   wd (HypertoroidalWDDistribution)
            %       current estimate            
            wd = this.wd;
        end
    end
    
end


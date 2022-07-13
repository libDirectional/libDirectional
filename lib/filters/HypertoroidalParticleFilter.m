classdef HypertoroidalParticleFilter < AbstractHypertoroidalFilter & AbstractParticleFilter
    % SIR particle filter on the hypertorus
    
    methods
        function this = HypertoroidalParticleFilter(nParticles, dim)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            arguments
                nParticles (1,1) double {mustBeInteger,mustBePositive}
                dim (1,1) double {mustBeInteger,mustBePositive}
            end
            this.dist = HypertoroidalWDDistribution(repmat((0:nParticles-1)/nParticles*2*pi,[dim,1]));
        end
        
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractHypertoroidalDistribution)
            %       new state            
            assert (isa (dist_, 'AbstractHypertoroidalDistribution'));
            if ~isa(dist_, 'HypertoroidalWDDistribution')
                dist_ = HypertoroidalWDDistribution(dist_.sample(length(this.dist.d)));
            end
            this.dist = dist_;
        end
        
        function predictNonlinear(this, f, noiseDistribution, functionIsVectorized)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   f (function handle)
            %       system function
            %   noiseDistribution (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            arguments
                this (1,1) HypertoroidalParticleFilter
                f (1,1) function_handle
                % Can be empty for no noise, therefore do not enforce (1,1)
                noiseDistribution AbstractHypertoroidalDistribution = HypercylindricalDiracDistribution.empty
                functionIsVectorized (1,1) logical = true
            end
            predictNonlinear@AbstractParticleFilter(this, f, noiseDistribution, functionIsVectorized, false);
            this.dist.d = mod(this.dist.d,2*pi);
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
            n = length(this.dist.d);
            noiseIds = discretesample(weights, n);
            d = zeros(this.dist.dim, n);
            for i=1:n
                % This loop is necessary for non-additve noise
                d(:,i) = f(this.dist.d(:,i),samples(noiseIds(i)));
            end
            this.dist.d=d;
        end
    end
    
end


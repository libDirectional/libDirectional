classdef (Abstract) AbstractParticleFilter < AbstractFilter
    % SIR particle filter on the hypertorus
    
    properties
        dist
    end
    
    methods
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractHypertoroidalDistribution)
            %       new state
            arguments
                this (1,1) AbstractParticleFilter
                dist_ (1,1) AbstractDistribution
            end
            % All automatic conversion can be done in the inherting classes
            assert(isa(dist_, class(this.dist)),'New distribution has to be of the same class as (or inhert from) the previous density.');
            this.dist = dist_;
        end
        
        function predictIdentity(this, noiseDistribution)
            arguments
                this (1,1) AbstractParticleFilter
                noiseDistribution (1,1) AbstractDistribution
            end
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            this.predictNonlinear(@(x) x, noiseDistribution);
        end
        
        function predictNonlinear(this, f, noiseDistribution, functionIsVectorized, shiftInsteadOfAdd)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k),
            % where w(k) is additive noise given by noiseDistribution.
            % Overwrite this function when it is no valid to use + in the
            % respective manifold!
            %
            % Parameters:
            %   f (function handle)
            %       system function
            %   noiseDistribution (AbstractHypertoroidalDistribution)
            %       distribution of additive noise
            arguments
                this (1,1) AbstractParticleFilter
                f (1,1) function_handle
                % Provide distribution or keep empty
                noiseDistribution = []
                functionIsVectorized (1,1) logical = true
                shiftInsteadOfAdd (1,1) logical = ~(isa(noiseDistribution,'AbstractHypertoroidalDistribution')||isa(noiseDistribution,'AbstractHypercylinddricalDistribution'))
            end
            assert(isempty(noiseDistribution) || this.dist.dim == noiseDistribution.dim);
            % Apply f
            if functionIsVectorized
                this.dist.d = f(this.dist.d);
            else
                this.dist = this.dist.applyFunction(f);           
            end
            if ~isempty(noiseDistribution)
                if ~shiftInsteadOfAdd
                    % Calculate effect of (additive) noise. Bounded dimensions have
                    % to be wrapped in subclasses
                    noise = noiseDistribution.sample(numel(this.dist.w));
                    this.dist.d = this.dist.d + noise;
                else
                    % We leave it up to the subclasses to check the noise
                    % is "neural".
                    % This means: All 0s for hypertoroidal, linear domains 
                    % and Cartesian products thereof; [0;...;0;1] for 
                    % hyperspherical/hyperhemispherical domains 
                    % (can be Cartesian product with linear domains)
                    for i=1:length(this.dist.d)
                        noiseCurr = noiseDistribution.setMode(this.dist.d(:,i));
                        this.dist.d(:,i) = noiseCurr.sample(1);
                    end
                end
            end 
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
            arguments
                this (1,1) HypercylindricalParticleFilter
                f (1,1) function_handle
                samples (:,:) double
                weights (1,:) double
            end
            %(samples, weights) are discrete approximation of noise
            assert(size(samples,2) == size(weights,2), 'samples and weights must match in size');
            
            weights = weights/sum(weights); %ensure normalization
            n = numel(this.dist.w);
            noiseIds = discretesample(weights, n);
            d = zeros(this.dist.dim, n);
            for i=1:n
                % This loop is necessary for non-additve noise
                d(:,i) = f(this.dist.d(:,i),samples(noiseIds(i)));
            end
            this.dist.d=d;
        end
        
        function updateIdentity(this, noiseDistribution, z, shiftInsteadOfAdd)
            arguments
                this (1,1) AbstractParticleFilter
                noiseDistribution (1,1) AbstractDistribution
                z (:,1) double
                shiftInsteadOfAdd (1,1) logical = ~(isa(noiseDistribution,'AbstractHypertoroidalDistribution')||isa(noiseDistribution,'AbstractHypercylinddricalDistribution'))
            end
            assert(isempty(z) || size(z,1)==noiseDistribution.dim);
            if ~shiftInsteadOfAdd
                likelihood = LikelihoodFactory.additiveNoiseLikelihood(@(x) x, noiseDistribution);
                this.updateNonlinear(likelihood, z);
            else
                noiseForLikelihood = noiseDistribution.setMode(z);
                likelihood = @(x)noiseForLikelihood.pdf(x);
                this.updateNonlinear(likelihood);
            end
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
            arguments
                this (1,1) AbstractParticleFilter
                likelihood {mustBeA(likelihood,{'function_handle','AbstractDistribution'})}
                z (:,1) double = []
            end
            % You can either use a likelihood depending on z and x
            % and specify the measurement as z or use a likelihood that
            % depends only on x and omit z.
            if isa(likelihood,'AbstractDistribution')
                assert(isempty(z), 'Cannot pass a density and a measurement. To assume additive noise, use updateIdentity.');
                likelihood = @(x)likelihood.pdf(x);
            end
                
            if isempty(z)
                this.dist = this.dist.reweigh(likelihood);
            else
                this.dist = this.dist.reweigh(@(x) likelihood(z,x));
            end
            % Resample
            this.dist.d = this.dist.sample(length(this.dist.d));
            % Reset weights to equal weights
            this.dist.w = 1/size(this.dist.d,2)*ones(1,size(this.dist.d,2));
        end
        
        function dist = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   dist (HypertoroidalWDDistribution)
            %       current estimate  
            arguments
                this (1,1) AbstractParticleFilter
            end
            dist = this.dist;
        end
        
        function likelihoodVal=associationLikelihood(this, likelihood)
            arguments
                this (1,1) AbstractParticleFilter
                likelihood (1,1) AbstractDistribution % Must be given as Distribution only depending on x (could be changed because may not be a density)
            end
            % see Florian Pfaff, Kailai Li, and Uwe D. Hanebeck,
            % Association Likelihoods for Directional Estimation
            % Proceedings of the 2019 IEEE International Conference on
            % Multisensor Fusion and Integration for Intelligent Systems (MFI 2019),
            % Taipei, Republic of China, May, 2019.
            likelihoodVal = sum(likelihood.pdf(this.getEstimate.d).*this.getEstimate.w);
        end
    end
    
end


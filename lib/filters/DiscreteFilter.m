classdef DiscreteFilter < AbstractCircularFilter
    % A discrete filter on the circle based on the wrapped
    % Dirac distribution with a grid of equidistant samples.
    %
    % Gerhard Kurz, Florian Pfaff, Uwe D. Hanebeck,
    % Discrete Recursive Bayesian Filtering on Intervals and the Unit Circle
    % Proceedings of the 2016 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI 2016), 
    % Baden-Baden, Germany, September 2016.    
    %
    % also mentioned in 
    % Florian Pfaff, Gerhard Kurz, and Uwe D. Hanebeck,
    % Multimodal Circular Filtering Using Fourier Series
    % Proceedings of the 18th International Conference on Information Fusion (Fusion 2015),
    % Washington, D.C., USA, July 2015.
    
    properties
        wd WDDistribution
    end
    
    methods
        function this = DiscreteFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use
            wd_ = WDDistribution((0:nParticles-1)/nParticles*2*pi);
            this.setState(wd_)
        end
        
        function setState(this, distribution)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (AbstractCircularDistribution)
            %       new state
            assert (isa (distribution, 'AbstractCircularDistribution'));
            if ~isa(distribution, 'WDDistribution')
                this.wd = WDDistribution(this.wd.d, distribution.pdf(this.wd.d));
            else
               nParticles = size(distribution.d,2);
               assert(all(distribution.d == (0:nParticles-1)/nParticles*2*pi));
               this.wd = distribution;               
            end
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
        
        function predictNonlinear(this, f, noiseDistribution, useProportional)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi)
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            %   useProportional (boolean)
            %       use proportional weight assignment method
            assert (isa (noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(f,'function_handle'));
            
            if nargin<=3
                useProportional = true;
            end
            
            if isa(noiseDistribution, 'WDDistribution')
                % fall back to nonadditive case, because noise is given by
                % samples anyway
               this.predictNonlinearNonAdditive(@(x,y) f(x) + y, noiseDistribution.d, noiseDistribution.w)      
               return
            end
            
            %apply f
            wdF = this.wd.applyFunction(f); 
            nParticles = length(wdF.d);
            particleDistance = 2*pi/nParticles;
            w_ = zeros(size(this.wd.w)); % set all weights to zero
            for i=1:nParticles
                id1 = 1+floor(wdF.d(i)/particleDistance); % id of particle to the left
                id2 = id1+1; %id of particle to the right
                if id2>nParticles %check for wrapping
                    id2 = 1;
                end
                distance1 = angularError(wdF.d(i),this.wd.d(id1));
                distance2 = angularError(wdF.d(i),this.wd.d(id2));
                if useProportional
                    w_(id1) = w_(id1) + distance2 / (distance1+distance2) * wdF.w(i);
                    w_(id2) = w_(id2) + distance1 / (distance1+distance2) * wdF.w(i);
                else
                    if distance1<distance2
                        w_(id1) = w_(id1) + wdF.w(i);
                    else
                        w_(id2) = w_(id2) + wdF.w(i);
                    end
                end
            end
            
            %add (additive) noise 
            wdNoise = WDDistribution(this.wd.d, arrayfun(@(x) noiseDistribution.pdf(x), this.wd.d));
            % this.wd = this.wd.convolution(wdNoise);
            % take advantage of the fact that number of components does not
            % increase
            
            % w_ = cconv (w_, wdNoise.w, nParticles);
            % workaround for bugs in cconv?!
            w_= conv (w_, wdNoise.w);
            w_ = w_(1:nParticles) + [w_(nParticles+1:end) 0];
            
            this.wd = WDDistribution(this.wd.d, w_);
        end
        
        function predictNonlinearNonAdditive(this, f, noiseSamples, noiseWeights)          
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
            
            % (samples, weights) are discrete approximation of noise
            assert(size(noiseWeights,1) == 1, 'weights most be row vector')
            assert(size(noiseSamples,2) == size(noiseWeights,2), 'samples and weights must match in size');
            assert(isa(f,'function_handle'));
            
            noiseWeights = noiseWeights/sum(noiseWeights); % ensure normalization

            % apply f
            nParticles = length(this.wd.d);
            particleDistance = 2*pi/nParticles;
            oldWeights = this.wd.w;
            this.wd.w = zeros(size(this.wd.w)); % set all weights to zero
            for i=1:nParticles
                for j=1:length(noiseSamples)
                    d = mod(f(this.wd.d(i),noiseSamples(j)), 2*pi);
                    id1 = 1+floor(d/particleDistance); % id of particle to the left
                    id2 = id1+1; %id of particle to the right
                    if id2>nParticles %check for wrapping
                        id2 = 1;
                    end
                    distance1 = angularError(d,this.wd.d(id1));
                    distance2 = angularError(d,this.wd.d(id2));
                    this.wd.w(id1) = this.wd.w(id1) + distance2 / (distance1+distance2) * oldWeights(i)*noiseWeights(j) ;
                    this.wd.w(id2) = this.wd.w(id2) + distance1 / (distance1+distance2) * oldWeights(i)*noiseWeights(j) ;
                end
            end            
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
            
            assert(isa(likelihood,'function_handle'));
            
            % You can either use a likelihood depending on z and x
            % and specify the measurement as z or use a likelihood that
            % depends only on x and omit z.
            if nargin==2
                this.wd = this.wd.reweigh(likelihood);
            else
                this.wd = this.wd.reweigh(@(x) likelihood(z,x));
            end
        end
        
        function wd = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   wd (WDDistribution)
            %       current estimate
            wd = this.wd;
        end
    end
    
end


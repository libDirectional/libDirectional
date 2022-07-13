classdef HypersphericalParticleFilter < AbstractHypersphericalFilter & AbstractParticleFilter
    % SIR particle filter on the hypersphere
    
    methods
        function this = HypersphericalParticleFilter(nParticles,dim)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            arguments
                nParticles (1,1) {mustBeInteger,mustBePositive}
                dim (1,1) {mustBeInteger,mustBePositive}
            end
            this.dist = HypersphericalDiracDistribution(repmat(eye(dim,1),1,nParticles));
        end
        
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   dist_ (HypersphericalDiracDistribution)
            %       new state            
            assert (isa (dist_, 'AbstractHypersphericalDistribution'));
            if ~isa(dist_, 'HypersphericalDiracDistribution')
                dist_ = HypersphericalDiracDistribution(dist_.sample(length(this.dist.d)));
            end
            this.dist = dist_;
        end
        
        function predictNonlinear(this, f, noiseDistribution, functionIsVectorized)
            % Predicts assuming a nonlinear system model.
            % Currently only supports VMFDistributions. mu is not used in
            % the calculation. To avoid a warning, use mu=[0;0;1] to indicate that
            % your noise is intended to be "neutral" concerning its mean.
            %
            % Parameters:
            %   f (function handle)
            %       system function. The function may include the adding of
            %       noise, in which case noiseDistribution should be set to
            %       'none'
            %   noiseDistribution (VMFDistribution)
            %       distribution of additive noise
            arguments
                this (1,1) AbstractParticleFilter
                f (1,1) function_handle
                % Provide distribution or keep empty
                noiseDistribution {mustBeA(noiseDistribution,{'VMFDistribution','WatsonDistribution','double'})} = VMFDistribution.empty
                functionIsVectorized (1,1) logical = true
            end
            assert(isempty(noiseDistribution) || this.dist.dim == noiseDistribution.dim);
            if ~isempty(noiseDistribution) && ~isequal(noiseDistribution.mu,[zeros(noiseDistribution.dim-1,1);1])
                warning('mu of noiseDistribution is not "neutral", which we define as [0;0;...;0;1].');
            end
            predictNonlinear@AbstractParticleFilter(this, f, noiseDistribution, functionIsVectorized, true);
        end
        
        function predictCustom(this, genNewSamples)
            % Takes custom function that generates new samples out of the old
            % ones. If there is system noise, respecting it is part of the
            % responsibility of genNewSamples.
            this.dist.d = genNewSamples(this.dist.d);
        end
        
        function updateIdentity(this, measNoise, z)
            % Updates assuming identity system model. Currently only 
            % supports VMF and WatsonDistribution and will be evaluated
            % according to
            % f(z_k|x_k) = VMF(x_k, kappa),
            % where kappa is the concentration of the VMF/Watson distribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractHypersphericalDistribution)
            %       distribution of additive noise
            %   z (column vector)
            %       measurement
            assert(isprop(measNoise,'mu'), 'Currently, only VMF- and Watson-distributed noise terms are supported.');
            
            if nargin==3
                if norm(measNoise.mu-[zeros(this.dim-1,1); 1]) > 1E-6
                    error('UpdateIdentity:UnexpectedMeas', 'z needs to be [0;...; 0; 1] to use updateIdentity.');
                end
                measNoise.mu=z;
            end
            disttmp=this.dist.reweigh(@(x)measNoise.pdf(x));
            
            this.dist.d = disttmp.sample(length(this.dist.d));
            this.dist.w = 1/size(this.dist.d,2)*ones(1,size(this.dist.d,2));
        end
        
    end
    
end

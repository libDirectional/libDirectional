classdef HypercylindricalParticleFilter < LinBoundedParticleFilter & AbstractHypercylindricalFilter
    % SIR particle filter for the hypercylinder
    methods
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
            if ~isa(dist_, 'HypercylindricalDiracDistribution')
                dist_ = HypercylindricalDiracDistribution.fromDistribution(dist_, numel(this.dist.w));
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
            %   noiseDistribution (AbstractHypercylindricalDistribution)
            %       distribution of additive noise
            arguments
                this (1,1) HypercylindricalParticleFilter
                f (1,1) function_handle
                % Can be empty for no noise, therefore do not enforce (1,1)
                noiseDistribution AbstractHypercylindricalDistribution {mustBeScalarOrEmpty} = HypercylindricalDiracDistribution.empty
                functionIsVectorized (1,1) logical = true
            end
            predictNonlinear@AbstractParticleFilter(this, f, noiseDistribution, functionIsVectorized, false);
            % Wrap bounded dimensions
            this.dist.d(1:this.dist.boundD,:) = mod(this.dist.d(1:this.dist.boundD,:),2*pi);
        end 
    end
end


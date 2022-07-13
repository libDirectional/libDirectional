classdef SE3ParticleFilter < LinBoundedParticleFilter & AbstractSE3Filter
    methods
        function this = SE3ParticleFilter(nParticles)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            arguments
                nParticles (1,1) double {mustBeInteger, mustBePositive}
            end
            this@LinBoundedParticleFilter(nParticles, 4, 3);
        end
        
        function setState(this, dist_)
            % Sets the current system state
            %
            % Parameters:
            %   distribution
            %       new state
            arguments
                this (1,1) SE3ParticleFilter
                dist_ (1,1) AbstractLinBoundedDistribution
            end
            if ~isa(dist_, 'SE3DiracDistribution')
                assert(dist_.boundD==4&&dist_.linD==3);
                dist_ = SE3DiracDistribution.fromDistribution(dist_, numel(this.dist.w));
            end
            this.dist = dist_;
        end

        function predictNonlinear(this, f, noiseDistribution, functionIsVectorized)
            % Predicts assuming a nonlinear system model.
            % Currently only supports VMFDistributions. mu is not used in
            % the calculation. To avoid a warning, use mu=[0;0;1] for the hemispherical dimesnion
            % to indicate that your noise is intended to be "neutral" concerning its mean.
            %
            % Parameters:
            %   f (function handle)
            %       system function. The function may include the adding of
            %       noise, in which case noiseDistribution should be set to
            %       'none'
            %   noiseDistribution
            %       distribution of additive noise
            arguments
                this (1,1) AbstractParticleFilter
                f (1,1) function_handle
                % Provide distribution or keep empty
                noiseDistribution = []
                functionIsVectorized (1,1) logical = true
            end
            assert(isempty(noiseDistribution) || ismethod(noiseDistribution,'setMode'));
            if ~isempty(noiseDistribution) && ~isequal(noiseDistribution.mode(),[zeros(this.dist.boundD-1,1);1;zeros(this.dist.linD,1)])
                warning('mu of noiseDistribution is ignored...');
            end
            predictNonlinear@AbstractParticleFilter(this, f, noiseDistribution, functionIsVectorized, true);
        end

        function est = getPointEstimate(this)
            arguments
                this (1,1) SE3ParticleFilter
            end
            estLin = this.dist.marginalizePeriodic().mean();

            distHemisphere = this.dist.marginalizeLinear();
            hpfTmp = HyperhemisphericalParticleFilter(size(this.dist.d,2), 4);
            hpfTmp.setState(distHemisphere);
            
            est = [hpfTmp.getPointEstimate() ;estLin];
        end
    end
    
end


classdef HypersphericalParticleFilter < AbstractHypersphericalFilter
    % SIR particle filter on the hypersphere
    
    properties
        wd
    end
    
    methods
        function this = HypersphericalParticleFilter(nParticles,dim)
            % Constructor
            %
            % Parameters:
            %   nParticles (integer > 0)
            %       number of particles to use 
            %   dim (integer >0)
            %       dimension
            this.wd = HypersphericalDiracDistribution(repmat(eye(dim,1),1,nParticles));
            this.dim=dim;
        end
        
        function setState(this, wd_)
            % Sets the current system state
            %
            % Parameters:
            %   wd_ (HypersphericalDiracDistribution)
            %       new state            
            assert (isa (wd_, 'AbstractHypersphericalDistribution'));
            if ~isa(wd_, 'HypersphericalDiracDistribution')
                wd_ = HypersphericalDiracDistribution(wd_.sample(length(this.wd.d)));
            end
            this.wd = wd_;
        end
        
        function predictIdentity(this, noiseDistribution)
            % Predicts assuming identity system model. Currently only 
            % supports VMFDistributions. i.e.,
            % f(x(k+1)|x_k) = VMF(x_k, kappa),
            % where kappa is the concentration of the VMF distribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractHypersphericalDistribution)
            %       distribution of additive noise
            this.predictNonlinear(@(x) x, noiseDistribution);
        end
        
        function predictNonlinear(this, f, noiseDistribution)
            % Predicts assuming a nonlinear system model.
            % Currently only supports VMFDistributions. mu is not used in
            % the calculation. To avoid a warning, use mu=[0;0;1] to indicate that
            % your noise is intended to be "zero mean".
            %
            % Parameters:
            %   f (function handle)
            %       system function
            %   noiseDistribution (VMFDistribution)
            %       distribution of additive noise
            
            assert (isa (noiseDistribution, 'VMFDistribution'),'Currently, only VMF distributed noise terms are allowed.');
            if ~isequal(noiseDistribution.mu,[0;0;1])
                warning('mu of noiseDistribution is being ignored...');
            end
            assert(isa(f,'function_handle'));
            %apply f
            wdF = this.wd.applyFunction(f);
            %calculate effect of (additive) noise
            noiseModified=noiseDistribution;
            for i=1:length(this.wd.d)
                noiseModified.mu=this.wd.d(:,i);
                wdF.d(:,i)=noiseModified.sample(1);
            end
            this.wd = wdF;
        end
        
        function updateIdentity(this, noiseDistribution, z)
            % Updates assuming identity system model. Currently only 
            % supports VMFDistributions. i.e.,
            % f(z_k|x_k) = VMF(x_k, kappa),
            % where kappa is the concentration of the VMF distribution.
            %
            % Parameters:
            %   noiseDistribution (AbstractHypersphericalDistribution)
            %       distribution of additive noise
            %   z (column vector)
            %       measurement
            assert (isa (noiseDistribution, 'VMFDistribution'),'Currently, only VMF distributed noise terms are supported.');
            
            if nargin==3
                noiseDistribution.mu=z;
                warning('updateIdentity:muReplaced','mu of noiseDistribution is replaced by measurement...');
            end
            wdtmp=this.wd.reweigh(@(x)noiseDistribution.pdf(x));
            
            this.wd.d = wdtmp.sample(length(this.wd.d));
            this.wd.w = 1/size(this.wd.d,2)*ones(1,size(this.wd.d,2));
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
            % 
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
            %   wd (HypersphericalDiracDistribution)
            %       current estimate            
            wd = this.wd;
        end
        
        function mean=getEstimateMean(this)
            % Return mean of current estimate 
            %
            % Returns:
            %   mean (column vector)
            %       current estimate mean (normalized)
            vecSum = sum(this.wd.d.*repmat(this.wd.w,this.dim,1),2);
            mean = vecSum/norm(vecSum);
        end
    end
    
end

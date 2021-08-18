classdef VMFFilter < AbstractHypersphericalFilter
    % Filter based von the VMF distribution.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Unscented von Mises-Fisher Filtering
    % IEEE Signal Processing Letters, 2016.    
    
    properties
        state (1,1) VMFDistribution = VMFDistribution([1;0],1)
    end
        
    methods 
        function setState(this, state_)
            assert(isa(state_, 'VMFDistribution'));
            this.state = state_;
        end
        
        function est = getEstimate(this)
           est = this.state; 
        end
        
        function predictIdentity(this, sysNoise)
            % Predicts assuming identity system model, i.e.,
            % f(x(k+1)|x_k) = VMF(x(k+1) ; mu = Q_k * x_k, kappa_k^w)
            % where w is noise given by sysNoise.
            %
            % Parameters:
            %   sysNoise (VMFDistribution)
            %       distribution of noise
            assert(isa(sysNoise, 'VMFDistribution'));
            
            % todo: implement use of rotation Q_k
            
            this.state = this.state.convolve(sysNoise);
        end
        
        function predictNonlinear(this, a, sysNoise)
            % Predicts assuming nonlinear system model, i.e.,
            % f(x(k+1)|x_k) = VMF(x(k+1) ; mu = a_k(x_k), kappa_k^w)
            % where w is noise given by sysNoise.
            %
            % Parameters:
            %   a (function handle)
            %       function S^(d-1) -> S^(d-1)
            %   sysNoise (VMFDistribution)
            %       distribution of noise            
            assert(isa(sysNoise, 'VMFDistribution'));
            assert(isa(a, 'function_handle'));
            
            stateSamples = this.state.sampleDeterministic();
            
            for i=1:size(stateSamples,2)
                stateSamples(:,i) = a(stateSamples(:,i));
            end
            
            vmf = VMFDistribution.fit(stateSamples);
            this.state = vmf.convolve(sysNoise);
        end
        
        function predictNonlinearArbitraryNoise(this, a, noiseSamples, noiseWeights)
            % Predicts assuming nonlinear system model with arbitrary noise
            %
            % Parameters:
            %   a (function handle)
            %       function S^(d-1) x W -> S^(d-1) with noise space W
            %   noiseSamples (m x n)
            %       samples of the noise (where m = dim W)
            %   noiseWeights (1 x n)
            %       weight for each sample
            
            assert(isa(a, 'function_handle'));
            assert(size(noiseSamples,2) == size(noiseWeights,2));
            assert(size(noiseWeights,1) == 1);
            assert(all(noiseWeights>0));
            
            noiseWeights = noiseWeights/sum(noiseWeights); %normalize weights
            
            stateSamples = this.state.sampleDeterministic();
            k = 1;
            for i=1:size(stateSamples,2)
                for j=1:size(noiseSamples,2)
                    newSamples(:,k) = a(stateSamples(:,i), noiseSamples(:,j));
                    newWeights(k) = noiseWeights(j)/size(stateSamples,2);
                    k = k+1;
                end
            end
            
            this.state = VMFDistribution.fit(newSamples, newWeights);
        end
        
        function updateIdentity(this, measNoise, z)
            % Updates assuming identity measurement model, i.e.,
            % f(z(k)|x(k=) = VMF(z(k); mu(k) = Q_k x_k, kappa_k^v)
            % where v is noise given by measNoise.
            %
            % Parameters:
            %   measNoise (VMFDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)            
            
            % Mu of measNoise is discarded and the mean is set to z. This
            % operation is not ambiguous as measNoise is a
            % VMFDistribution and is thus zonal. To acknowledge that the mean 
            % is disregarded, the last entry of mu has to be set to 1 and 
            % all other to 0.
            
            assert(isa(measNoise, 'VMFDistribution'));
            assert(measNoise.mu(end)==1,'Set mu of measNoise to [0,0,...,1] to acknowledge that the mean is discarded');
            assert(size(z,1) == this.dim);
            assert(size(z,2) == 1);
            
            %todo how to handle rotation Q 
            
            measNoise.mu = z;
            
            this.state = this.state.multiply(measNoise);
        end
        
        function updateNonlinearProgressive(~)
            %todo later
            error('not implemented');
        end
        
        function mean = getEstimateMean(this)
            mean = this.state.mu;
        end
    end
end

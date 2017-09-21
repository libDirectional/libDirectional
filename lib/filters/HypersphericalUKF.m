classdef HypersphericalUKF < AbstractHypersphericalFilter
    
    properties
        state
        ukf
    end
    
    methods
        function this = HypersphericalUKF()
            % Constructor
            state_ = GaussianDistribution(1, 1);
            this.setState(state_)
            this.ukf = UKF();
        end
        
        function setState(this, state_)
            % Sets the current system state
            %
            % Parameters:
            %   state_ (GaussianDistribution)
            %       new state (1D Gaussian)
            if ~isa(state_,'GaussianDistribution')
                state_ = state_.toGaussian();
            end
           this.state = state_;
        end
        
        function predictNonlinear(this, f, gaussSys) 
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k) 
            % where w(k) is additive noise given by gaussSys.
            %
            % Parameters:
            %   f (function handle)
            %       function from R^(d-1) to R^(d-1)
            %   gaussSys (GaussianDistribution)
            %       distribution of additive noise (warning: mean is
            %       ignored)
            assert(isa(f,'function_handle'));
            if ~isa(gaussSys,'GaussianDistribution')
                gaussSys = gaussSys.toGaussian();
            end
            
            function y=g(x)
                x = x/norm(x);
                y = f(x);
                y = y/norm(y);
            end
            model = SysModelWrapper(@g);
            model.setNoise(Gaussian(zeros(size(gaussSys.mu)), gaussSys.C));
            
            this.ukf.setState(Gaussian(this.state.mu, this.state.C));
            this.ukf.predict(model);
            [this.state.mu, this.state.C] = this.ukf.getStateMeanAndCov();
            
            % normalize mean
            if norm(this.state.mu)~=0
                this.state.mu = this.state.mu/norm(this.state.mu);
            else
                error('mu was 0');
            end
        end
        
        function predictNonlinearArbitraryNoise(this, f, noiseSamples, noiseWeights)
            % Updates assuming nonlinear measurement model
            % z(k) = f(x(k), v_k)
            % where z is the measurement and v_k is the non-additive measurement noise
            % given by gaussMeas
            assert(isa(f, 'function_handle'));
            assert(size(noiseSamples,2) == size(noiseWeights,2));
            assert(size(noiseWeights,1) == 1);
            assert(all(noiseWeights>0));
            
            noiseWeights = noiseWeights/sum(noiseWeights); %normalize weights
            
            % get UKF samples
            ukfSampling = GaussianSamplingUKF();
            [stateSamples, stateWeights, numGaussianSamples] = ukfSampling.getSamples(Gaussian(this.state.mu, this.state.C));
            if isscalar(stateWeights)
                stateWeights = repmat(stateWeights, 1, numGaussianSamples);
            end
            k = 1;
            for i=1:size(stateSamples,2)
                for j=1:size(noiseSamples,2)
                    newSamples(:,k) = f(stateSamples(:,i), noiseSamples(:,j));
                    newWeights(k) = noiseWeights(j)*stateWeights(i);
                    k = k + 1;
                end
            end
            newWeights = newWeights/sum(newWeights);
                        
            % Compute predicted state mean and covariance
            [predictedStateMean, ...
             predictedStateCov] = Utils.getMeanAndCov(newSamples, newWeights);
            
            predictedStateMean = predictedStateMean/norm(predictedStateMean);
            this.state = GaussianDistribution(predictedStateMean,predictedStateCov);
        end
        
        function updateNonlinear(this, f, gaussMeas, z)
            % Updates assuming nonlinear measurement model
            % z(k) = f(x(k)) + v_k
            % where z is the measurement and v_k is the additive measurement noise
            % given by gaussMeas
            % 
            % Parameters:
            %   f (function handle)
            %       function from R^(d-1) to R^n, where R^n is
            %       the measurement space containing z
            %   gaussMeas (gaussianDistribution)
            %       distribution of additive noise (warning: mean is
            %       ignored)
            %   z (n x 1 vector)
            %       measurement            
            if ~isa(gaussMeas,'GaussianDistribution')
                gaussMeas = gaussMeas.toGaussian();
            end
            
            %UKF Update
            function y=g(x)
                x = x/norm(x);
                y = f(x);
            end
            model = MeasModelWrapper(@g);
            model.setNoise(Gaussian(zeros(size(gaussMeas.mu)), gaussMeas.C));
            
            this.ukf.setState(Gaussian(this.state.mu, this.state.C));
            this.ukf.update(model, z);
            [this.state.mu, this.state.C] = this.ukf.getStateMeanAndCov();
            
            if norm(this.state.mu)~=0 
                this.state.mu = this.state.mu/norm(this.state.mu);
            else
                error('mu was 0');
            end
        end
        
        function updateIdentity(this, gaussMeasNoise, z)
            this.updateNonlinear(@(x) x, gaussMeasNoise, z);
        end
        
        function predictIdentity(this, gaussSys)
            this.predictNonlinear(@(x) x, gaussSys);
        end

        function state = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   state (GaussianDistribution)
            %       current estimate
            state = this.state;
        end
        
        function estimateMean = getEstimateMean(this)
            % Return current estimate 
            %
            % Returns:
            %   estimateMean (n x 1 column vector)
            %       mean of current estimate 
            estimateMean = this.state.mu;
        end
    end
    
end


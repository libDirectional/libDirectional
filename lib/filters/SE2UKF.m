classdef SE2UKF < AbstractSE2Filter
    % Applies UKF to dual-quaternion-based Estimation on SE(2).
    %
    % A Stochastic Filter for Planar Rigid-Body Motions,
    % Igor Gilitschenski, Gerhard Kurz, and Uwe D. Hanebeck,
    % Proceedings of the 2015 IEEE International Conference on Multisensor
    % Fusion and Integration for Intelligent Systems (MFI),
    % San Diego, USA, 2015
    
    properties (Access = public)
        gauss % Current state estimate as Gaussian
    end
    
    methods
        function this = SE2UKF()
            % Constructs an SE2 UKF object.
        end
        
        function setState(this, state_)
            assert(isa(state_, 'GaussianDistribution'));
            assert(all(size(state_.C)==[4 4]), 'Wrong dimension.');
            assert(norm(state_.mu(1:2))==1, 'First two entries of estimate must be normalized');
            
            this.gauss = state_;
        end
        
        function this = predictIdentity(this, gaussSys)
            % Prediction step using a noisy identity model.
            %
            %  This prediction step assumes the system to evolve according to
            %        x_{t+1} = v [+] x_t ,
            %  where v is the system noise and x_t is the current system state.
            %  The model assumes the mean of the noise to be zero.
            
            assert(isa(gaussSys, 'GaussianDistribution'));
            assert(all(size(gaussSys.C)==[4 4]), 'System covariance matrix has wrong size');
            assert(all(size(gaussSys.mu)==[4,1]), 'Mean of system noise must be 4d column vector.');
            assert(norm(gaussSys.mu(1:2))==1, 'First two entries of the mean of the system noise must be normalized');
                        
            % Sampled estimates.
            systemStateSamples = [zeros(4,1) 2*chol(this.gauss.C)' -2*chol(this.gauss.C)'];
            systemStateSamples = systemStateSamples+repmat(this.gauss.mu,1,9);
            
            % Normalization of Measurement noise samples is part of the system
            % function.
            Norm=sqrt(sum(systemStateSamples(1:2,:).^2,1));
            systemStateSamples(1:2,:)=systemStateSamples(1:2,:)./repmat(Norm,2,1);
            
            % Compute measurement noise samples.
            systemNoiseSamples = [zeros(4,1) chol(4*gaussSys.C)' -chol(4*gaussSys.C)'];
            % Currently the system noise is considered to be zero mean (in the
            % manifold sense).
            systemNoiseSamples = systemNoiseSamples + repmat(gaussSys.mu, 1, 9);
            
            % Normalization of Measurement noise samples is part of the system
            % function.
            Norm=sqrt(sum(systemNoiseSamples(1:2,:).^2,1));
            systemNoiseSamples(1:2,:)=systemNoiseSamples(1:2,:)./repmat(Norm,2,1);
            
            
            % Compute predicted samples.
            % This step can be optimized by computing the correct covariance matrix as
            % we did in our Bingham Filter.
            systemPredictionSamples = ones(4,81);
            for i=1:9
                for j=1:9
                    systemPredictionSamples(:,(i-1)*9+j) = ...
                        SE2.dualQuaternionMultiply(systemStateSamples(:,i),systemNoiseSamples(:,j));
                end
            end
            
            %CP = systemStateSamples*diag(ones(1,9)/9)*systemStateSamples'+CW;
            CP = systemPredictionSamples*diag(ones(1,81)/81)*systemPredictionSamples';
            this.gauss.C = (CP+CP')/2; % Avoid numerical problems.
            this.gauss.mu = mean(systemStateSamples,2);
            
            % Renormalization.
            this.gauss.mu(1:2) = this.gauss.mu(1:2) / norm(this.gauss.mu(1:2));
        end
        
        
        function this = updateIdentity(this, gaussMeas, z)
            % Incorporates measurement into current estimate.
            %
            %  Parameters:
            %    z - Measurement given as a 4d dual quaternion restricted to
            %        SE(2)
            
            assert(all(size(gaussMeas.C)==[4 4]), 'Observation covariance matrix has wrong size');
            assert(all(size(gaussMeas.mu)==[4,1]), 'mu must be 4d column vector.');
            assert(norm(gaussMeas.mu(1:2))==1, 'First two entries of the mean of the observation noise must be normalized');
            
            % Take closer quaternion to improve estimate.
            if (norm(z-this.gauss.mu)>norm(-z-this.gauss.mu))
                z = -z;
            end
            
            % Compute covariance of augmented state and generate samples.
            CAUG = [this.gauss.C zeros(4,4); zeros(4,4) gaussMeas.C];
            xaug = [this.gauss.mu; gaussMeas.mu];
            augmentedSamplesPlain = [zeros(8,1), chol(8*CAUG)' -chol(8*CAUG)'];
            augmentedSamples = augmentedSamplesPlain + repmat(xaug,1,17);
            
            % Normalize state samples
            %Norm=sqrt(sum(augmentedSamples(1:2,:).^2,1));
            %augmentedSamples(1:2,:) = augmentedSamples(1:2,:)./repmat(Norm,2,1);
            
            % Normalize samples
            Norm=sqrt(sum(augmentedSamples(5:6,:).^2,1));
            normalizedNoiseSamples = augmentedSamples(5:6,:)./repmat(Norm,2,1);
            normalizedNoiseSamples = [ normalizedNoiseSamples; augmentedSamples(7:8,:)];
            %augmentedSamples(5:6,:) = augmentedSamples(5:6,:)./repmat(Norm,2,1);
            
            % Apply measurement function.
            measurementSamples = zeros(4,17);
            for i=1:17
                measurementSamples(:,i) = ...
                    SE2.dualQuaternionMultiply(augmentedSamples(1:4,i),normalizedNoiseSamples(:,i));
            end
            
            % Compute Covariance matrices.
            measurementSamplesDiv = measurementSamples-repmat(mean(measurementSamples,2), 1,17);
            C=augmentedSamples*measurementSamplesDiv';
            PY=cov(measurementSamplesDiv',1);
            PXY = C/17;
            
            
            K = PXY/PY;
            xeaug = xaug + K*(z-mean(measurementSamples,2));
            xe = xeaug(1:4);
            CEAUG = CAUG - K * PY * K';
            this.gauss.C = CEAUG(1:4,1:4);
            %eig(this.curCov)
            
            xe(1:2) = xe(1:2) / norm(xe(1:2));
            this.gauss.mu = xe; % Ensure quaternion to be on the unit sphere as in LaVoila (2003).
            
            % UKF without state augmentation
            % PY = CP+CV;
            % PXY = CP;
            % K = PXY/PY;
            % xe =xp + K*(z-xp);
            % CE = CP - K * PY * K';
            % xe = xe / norm(xe); % Ensure quaternion to be on the unit sphere as in LaVoila (2003).
        end
        
        function est = getEstimate(this)
            est = this.gauss ;
        end
        
        function est = getPointEstimate(~) %#ok<STOUT>
            error('Not yet implemented.');
        end
    end
    
end


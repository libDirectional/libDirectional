classdef SE2BinghamFilter < AbstractSE2Filter
    % Perfoms state estimation with a Bingham-like SE(2) distribution.
    %
    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
    % A New Probability Distribution for Simultaneous Representation of Uncertain Position and Orientation
    % Proceedings of the 17th International Conference on Information Fusion (Fusion 2014), Salamanca, Spain, July 2014.
    %
    % A Stochastic Filter for Planar Rigid-Body Motions,
    % Igor Gilitschenski, Gerhard Kurz, and Uwe D. Hanebeck,
    % Proceedings of the 2015 IEEE International Conference on Multisensor
    % Fusion and Integration for Intelligent Systems (MFI),
    % San Diego, USA, 2015 
    
    properties (Access = public)
        curEstimate SE2BinghamDistribution % SE2 Distribution object representing current estimate.
    end
       
    methods
        function this = SE2BinghamFilter()
            % Constructs an SE2 Filter Object.
            this.curEstimate = SE2BinghamDistribution( -eye(4,4));
        end
        
        function setState(this, distribution)
            % Sets the current system state
            %
            % Parameters:
            %   distribution (SE2BinghamDistribution)
            %       new state
            arguments
                this (1,1) SE2BinghamFilter
                distribution (1,1) SE2BinghamDistribution
            end
            this.curEstimate = distribution;
        end
        
        function this = predictIdentity(this, sysNoise)
            % Prediction step using a noisy identity model.
            %
            %  This prediction step assumes the system to evolve according to
            %        x_{t+1} = x_t [+] v,
            %  where v is the system noise and x_t is the current system state.
            %
            % Parameters:
            %   sysNoise (SE2BinghamDistribution)
            %       system noise v
            %
            arguments
                this (1,1) SE2BinghamFilter
                sysNoise (1,1) SE2BinghamDistribution
            end
            [eSamples, eWeights] = this.curEstimate.sampleDeterministic();
            [nSamples, nWeights] = sysNoise.sampleDeterministic();
            
            % Number of samples (these values may differ once more sampling
            % schemes are implemented).
            eNumSamples = length(eWeights);
            nNumSamples = length(nWeights);
            
            rSamples = zeros(4,eNumSamples*nNumSamples);
            rWeights = zeros(1,eNumSamples*nNumSamples);
            
            % Perform dual quaternion multiplication
            % This needs vectorization for better perfomance
            k = 1;
            for i=1:eNumSamples
                for j=1:nNumSamples
                    rSamples(:,k) = SE2.dualQuaternionMultiply(eSamples(:,i),nSamples(:,j));
                    
                    rWeights(k) = eWeights(i)*nWeights(j);
                    k = k+1;
                end
            end
            
            % Fit new samples.
            this.curEstimate = SE2BinghamDistribution.fit(rSamples); %todo: weights are ignored!
        end
        
        function this = predictNonlinearNonAdditive(this,a,sysNoise)
            % nonlinear prediction step  with system dynamics
            % described as x_{k+1} = a(x_k, w)
            % Input:
            %   a(function handle): system model
            %   sysNoise: system noise of arbitrary form
            arguments
                this (1,1) SE2BinghamFilter
                a (1,1) function_handle
                sysNoise (1,1) SE2BinghamDistribution % Otherwise, samples would need to be different
            end
            assert(nargin(a)==2, 'System function should take 2 input arguments for the prediction with nonadditive noise.');

            [eSamples, eWeights] = this.curEstimate.sampleDeterministic();
            [vSamples, vWeights] = sysNoise.sampleDeterministic();

            % number of samples 
            eNumSamples = numel(eWeights);
            vNumSamples = numel(vWeights);

            rSamples = zeros(4, eNumSamples*vNumSamples);
            rWeights = zeros(1, eNumSamples*vNumSamples);

            for i=1:eNumSamples
            	for j=1:vNumSamples
                    rSamples(:,j+vNumSamples*(i-1)) = a(eSamples(:,i),vSamples(:,j));
                	rWeights(:,j+vNumSamples*(i-1)) = eWeights(:,i)*vWeights(:,j);
            	end
            end

            % Fit new samples.
            this.curEstimate = SE2BinghamDistribution.fit(rSamples,rWeights);
        end
        
        function this = predictNonlinear(this, a, sysNoise, functionUsesAnglePos)
            arguments
                this (1,1) SE2BinghamFilter
                a (1,1) function_handle
                sysNoise (1,1) SE2BinghamDistribution % Otherwise, samples would need to be different
                functionUsesAnglePos (1,1) logical = true % Set to false if function takes quaternion instead
            end
            assert(nargin(a)==1, 'System function should take 1 input arguments for the prediction with additive noise.');
            function xkkQuaternion = axw(xk, wk)
                if functionUsesAnglePos
                    xkProp = a(AbstractSE2Distribution.dualQuaternionToAnglePos(xk));
                else
                    xkProp = AbstractSE2Distribution.dualQuaternionToAnglePos(a(xk));
                end
                xkk = xkProp + AbstractSE2Distribution.dualQuaternionToAnglePos(wk);
                xkk(1,:) = mod(xkk(1,:), 2*pi);
                xkkQuaternion = AbstractSE2Distribution.anglePosToDualQuaternion(xkk);
            end
            this.predictNonlinearNonAdditive(@axw, sysNoise);
        end
        
        function this = updateIdentity(this, measNoise, z)
            % Incorporates measurement into current estimate.
            %
            %  Parameters:
            %     measNoise (SE2BinghamDistribution)
            %       measurement noise
            %    z (4 x 1 vector)
            %       Measurement given as a 4d dual quaternion restricted to
            %        SE(2)
            arguments
                this (1,1) SE2BinghamFilter
                measNoise (1,1) SE2BinghamDistribution
                z (4,1) double
            end
            
            % Auxiliary quantities
            D = diag([1 -1 -1 -1]);
            zMatrix = SE2.dualQuaternionToMatrix(D*z);
            
            % New Parameter matrix.
            newC = D*zMatrix'*measNoise.C*zMatrix*D ...
                + this.curEstimate.C;
            
            % Fix symmetry
            newC = 0.5 * (newC + newC');
            
            % Measurement update
            this.curEstimate = SE2BinghamDistribution(newC);
        end
        
        function this = updateNonlinear(this, likelihood, z)
            % updateProgressive handles nonlinear cases, forward call to
            % that function
            this = updateProgressive(this, likelihood, z);
        end
        
        function this = updateProgressive(this, likelihood, z, tau)
            % update step that progressively fuse the measurement given the likelihood
            % Reference goes to:
            % Kailai Li, Gerhard Kurz, Lukas Bernreiter, Uwe D. Hanebeck,
            % Nonlinear Progressive Filtering for SE(2) Estimation,
            % Proceedings of the 21th International Conference on Information Fusion (Fusion 2018), 
            % Cambridge, United Kingdom, July 2018.
            % Input:
            %   likelihood(function handle): function from Zx(S^1xR^2) to
            %                                [0, infinity), where Z is the
            %                                measurement space containing z
            %   z(arbitrary): current measurement 
            %   tau(scalar): threshold of controlling progression step, should be
            %                between 0 and 1
            arguments
                this (1,1) SE2BinghamFilter
                likelihood (1,1) function_handle
                z (:,1) double % Measurement does not need to be in state space
                tau (1,1) double = 0.02
            end
            assert(nargin(likelihood)==2);

            lambda = 1;
            steps = 0;
            while lambda > 0
                [samples, weights] = this.curEstimate.sampleDeterministic();
                l = zeros(size(weights));
                for i=1:size(samples,2)
                    l(i) = likelihood(z,samples(:,i));  
                end        
                lmin = min(l);
                lmax = max(l);
                if lmax == 0
                    warning('Progressive update failed because likelihood is 0 everwhere');
                    return
                end
                currentLambda = min(log(tau)/log(lmin/lmax),lambda);
                if currentLambda <= 0
                    warning('Progressive update with given threshold impossible')
                    currentLambda = 0.001;
                end
                weightsNew = weights .* l.^currentLambda;
                weightsNew = weightsNew/sum(weightsNew);
                this.curEstimate = SE2BinghamDistribution.fit(samples,weightsNew);
                lambda = lambda - currentLambda;
                steps = steps + 1;
            end      
        end

        function est = getEstimate(this)
            arguments
                this (1,1) SE2BinghamFilter
            end
            % Return current estimate 
            %
            % Returns:
            %   est (SE2BinghamDistribution)
            %       current estimate
            est = this.curEstimate;
        end
        
        function est = getPointEstimate(this,convertToAnglePos)
            arguments
                this (1,1) SE2BinghamFilter
                convertToAnglePos (1,1) logical = true
            end
            est = this.curEstimate.mode();
            if convertToAnglePos
                [est(1),est(2:3)] = AbstractSE2Distribution.dualQuaternionToAnglePos(est);
                est(end) = [];
            end
        end
    end
end


classdef SE2BinghamFilter < handle
    % Perfoms state estimation with a Bingham-like SE(2) distribution.
    %
    % Igor Gilitschenski, Gerhard Kurz, Simon J. Julier, Uwe D. Hanebeck,
    % A New Probability Distribution for Simultaneous Representation of Uncertain Position and Orientation
    % Proceedings of the 17th International Conference on Information Fusion (Fusion 2014), Salamanca, Spain, July 2014.
    %
    % A Stochastic Filter for Planar Rigid-Body Motions,
    % Igor Gilitschenski, Gerhard Kurz, and Uwe D. Hanebeck,
    % Proceedings of the 2015 IEEE International Conference on Multisensor
    % Fusion and Integration for Intelligent Systems (MFI),
    % San Diego, USA, 2015 (submitted)    
    
    properties (Access = public)
        curEstimate % SE2 Distribution object representing current estimate.
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
            assert(isa(sysNoise, 'SE2BinghamDistribution'));
            
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
            this.curEstimate = SE2BinghamDistribution.fit(rSamples);
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
            assert(isa(measNoise, 'SE2BinghamDistribution'));
            assert(size(z,1) == 4);
            assert(size(z,2) == 1);
            
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
        
        function est = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   est (SE2BinghamDistribution)
            %       current estimate
            est = this.curEstimate;
        end
    end
end


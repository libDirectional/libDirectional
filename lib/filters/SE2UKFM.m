classdef SE2UKFM < AbstractSE2Filter
    properties
        state (1,1) struct
        C (3,3) double
        weights struct
    end
        
    methods
        function this = SE2UKFM(noiseDim, sigmaPointScaling)
            arguments
                noiseDim (1,1) double {mustBePositive,mustBeInteger}
                sigmaPointScaling (:,1) double = 10^-3*ones(3,1)% Should be between 10^-3 and 1
            end
            this.weights = ukf_set_weight(3, noiseDim, sigmaPointScaling); % sigmaPointScaling is alpha in the code by the authors
        end
        function setState(this, state_, C)
            arguments
                this (1,1) SE2UKFM
                state_
                C (3,3) double
            end
            if isa(state_,'struct')
                assert(isequal(size(state_.p),[2,1]));
                assert(isequal(size(state_.Rot),[2,2]));
                this.state = state_;
            else
                assert(isequal(size(state_),[3,1]));
                assert(isa(state_,'double'));
                rotMat = [cos(state_(1)),-sin(state_(1));sin(state_(1)),cos(state_(1))];
                this.state = struct('p',state_(2:3),'Rot',rotMat);
            end
            this.C = C;
        end
        
        function stateAndCov = getEstimate(this)
            arguments
                this (1,1) SE2UKFM
            end
            % We output the hybridMean as part of the state so we can use
            % .hybridMean of the state in the filter evaluation framework.
            stateAndCov = struct('state',this.state,'C',this.C,'hybridMean',this.getPointEstimate()); 
        end
        
        function mean = getPointEstimate(this)
            arguments
                this (1,1) SE2UKFM
            end
            mean = [atan2(this.state.Rot(2, 1),this.state.Rot(1, 1));this.state.p];
        end
        
        function predictIdentity(this, sysNoiseCov)
            arguments
                this (1,1) SE2UKFM
                sysNoiseCov (3,3) double
            end
            % Closely check the results for plausibility because this does not take an input!
            this.predictNonlinear(@(x, ~, ~, ~) x, sysNoiseCov);
        end
        
        function predictNonlinear(this, a, sysNoiseCov, u)
            arguments
                this (1,1) SE2UKFM
                a (1,1) function_handle
                sysNoiseCov (3,3) double
                u = []
            end
            Q = sysNoiseCov;
            cholQ = chol(Q);
            
            [this.state, this.C] = ukf_propagation(this.state, this.C, u, ...
                a, 1, @localization_phi, @localization_phi_inv, cholQ, this.weights);
        end
        
        function updatePositionMeasurement(this, measNoiseCov, z)
            arguments
                this (1,1) SE2UKFM
                measNoiseCov (:,:) double
                z (2,1) double
            end
            R = measNoiseCov;
            [this.state, this.C] = ukf_update(this.state, this.C, z, ...
                @(state)state.p, @localization_phi, R, this.weights);   
        end
        
        function updateIdentity(this, measNoiseCov, z)
            arguments
                this (1,1) SE2UKFM
                measNoiseCov (:,:) double
                z (3,1) double
            end
            R = measNoiseCov;
            [this.state, this.C] = ukf_update(this.state, this.C, z, ...
                @(state)[atan2(state.Rot(2, 1),state.Rot(1, 1));state.p],...
                @localization_phi, R, this.weights);   
        end
    end
end


classdef SecondOrderTaylorLinearGaussianFilter < LinearGaussianFilter
    % Abstract class for linear Gaussian filters that are based on the
    % second-order Taylor series approximation.
    %
    % SecondOrderTaylorLinearGaussianFilter Methods:
    %   SecondOrderTaylorLinearGaussianFilter - Class constructor.
    %   copy                                  - Copy a Filter instance.
    %   copyWithName                          - Copy a Filter instance and give the copy a new name/description.
    %   getName                               - Get the filter name/description.
    %   setColor                              - Set the filter color/plotting properties.
    %   getColor                              - Get the filter color/plotting properties.
    %   setState                              - Set the system state.
    %   getState                              - Get the system state.
    %   setStateMeanAndCov                    - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov                    - Get mean and covariance matrix of the system state.
    %   getStateDim                           - Get the dimension of the system state.
    %   predict                               - Perform a state prediction.
    %   update                                - Perform a measurement update.
    %   step                                  - Perform a combined state prediction and measurement update.
    %   setStateDecompDim                     - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                     - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing           - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing           - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing               - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing               - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold                - Set the measurement gating threshold.
    %   getMeasGatingThreshold                - Get the measurement gating threshold.
    
    % Literature:
    %   Michael Athans, Richard P. Wishner, and Anthony Bertolini,
    %   Suboptimal State Estimation for Continuous-Time Nonlinear Systems from Discrete Noisy Measurements,
    %   IEEE Transactions on Automatic Control, vol. 13, no. 5, pp. 504-514, Oct. 1968.
    
    % >> This function/class is part of the Nonlinear Estimation Toolbox
    %
    %    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
    %
    %    Copyright (C) 2015-2017  Jannik Steinbring <nonlinearestimation@gmail.com>
    %
    %    This program is free software: you can redistribute it and/or modify
    %    it under the terms of the GNU General Public License as published by
    %    the Free Software Foundation, either version 3 of the License, or
    %    (at your option) any later version.
    %
    %    This program is distributed in the hope that it will be useful,
    %    but WITHOUT ANY WARRANTY; without even the implied warranty of
    %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %    GNU General Public License for more details.
    %
    %    You should have received a copy of the GNU General Public License
    %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    methods
        function obj = SecondOrderTaylorLinearGaussianFilter(name)
            % Class constructor.
            %
            % Parameters:
            %   >> name (Char)
            %      An appropriate filter name / description of the implemented
            %      filter. The Filter subclass should set this during its
            %      construction to a meaningful default value (e.g., 'EKF'),
            %      or the user should specify an appropriate name (e.g.,
            %      'PF (10k Particles)').
            %
            % Returns:
            %   << obj (SecondOrderTaylorLinearGaussianFilter)
            %      A new SecondOrderTaylorLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@LinearGaussianFilter(name);
        end
    end
    
    methods (Sealed, Access = 'protected')
        % State prediction moment computation
        function [predictedStateMean, ...
                  predictedStateCov] = predictSysModel(obj, sysModel)
            [noiseMean, noiseCov, noiseCovSqrt] = sysModel.noise.getMeanAndCov();
            dimNoise = size(noiseMean, 1);
            
            % Compute system model derivatives around current state mean and noise mean
            [stateJacobian, ...
             noiseJacobian, ...
             stateHessians, ...
             noiseHessians] = sysModel.derivative(obj.stateMean, noiseMean);
            
            % Check computed derivatives
            obj.checkStateJacobian(stateJacobian, obj.dimState, obj.dimState);
            obj.checkStateHessians(stateHessians, obj.dimState, obj.dimState);
            
            obj.checkNoiseJacobian(noiseJacobian, obj.dimState, dimNoise);
            obj.checkNoiseHessians(noiseHessians, obj.dimState, dimNoise);
            
            [stateHessMean, ...
             stateHessCov, ...
             noiseHessMean, ...
             noiseHessCov] = obj.getHessianMomentsStateAndNoise(obj.dimState, stateHessians, obj.stateCov, ...
                                                                dimNoise, noiseHessians, noiseCov, obj.dimState);
            
            % Compute predicted state mean
            predictedStateMean = sysModel.systemEquation(obj.stateMean, noiseMean) + ...
                                 stateHessMean + noiseHessMean;
            
            % Compute predicted state covariance
            A = stateJacobian * obj.stateCovSqrt;
            B = noiseJacobian * noiseCovSqrt;
            
            predictedStateCov = A * A' + stateHessCov + B * B' + noiseHessCov;
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictAddNoiseSysModel(obj, sysModel)
            [noiseMean, noiseCov] = sysModel.noise.getMeanAndCov();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimNoise);
            
            % Compute system model derivatives around current state mean
            [stateJacobian, stateHessians] = sysModel.derivative(obj.stateMean);
            
            % Check computed derivatives
            obj.checkStateJacobian(stateJacobian, obj.dimState, obj.dimState);
            obj.checkStateHessians(stateHessians, obj.dimState, obj.dimState);
            
            [stateHessMean, ...
             stateHessCov] = obj.getHessianMomentsState(obj.dimState, stateHessians, ...
                                                        obj.stateCov, obj.dimState);
            
            % Compute predicted state mean
            predictedStateMean = sysModel.systemEquation(obj.stateMean) + stateHessMean + noiseMean;
            
            % Compute predicted state covariance
            A = stateJacobian * obj.stateCovSqrt;
            
            predictedStateCov = A * A' + stateHessCov + noiseCov;
        end
        
        % Helpers for the second-order Taylor-based linear measurement updates
        function [h, stateJacobian, ...
                  stateHessCov] = evaluateAddNoiseMeasModel(obj, measModel, dimMeas, ...
                                                            stateMean, stateCov)
            dimState = size(stateMean, 1);
            
            % Compute measurement model derivatives around current state mean
            [stateJacobian, stateHessians] = measModel.derivative(stateMean);
            
            % Check computed derivatives
            obj.checkStateJacobian(stateJacobian, dimMeas, dimState);
            obj.checkStateHessians(stateHessians, dimMeas, dimState);
            
            [stateHessMean, ...
             stateHessCov] = obj.getHessianMomentsState(dimState, stateHessians, stateCov, dimMeas);
            
            h = measModel.measurementEquation(stateMean) + stateHessMean;
        end
        
        function [h, stateJacobian, ...
                  noiseJacobian, ...
                  stateHessCov, ...
                  noiseHessCov] = evaluateMeasModel(obj, measModel, dimMeas, dimNoise, ...
                                                                  noiseMean, noiseCov, ...
                                                                  stateMean, stateCov)
            dimState = size(stateMean, 1);
            
            % Compute measurement model derivatives around current state mean and noise mean
            [stateJacobian, ...
             noiseJacobian, ...
             stateHessians, ...
             noiseHessians] = measModel.derivative(stateMean, noiseMean);
            
            % Check computed derivatives
            obj.checkStateJacobian(stateJacobian, dimMeas, dimState);
            obj.checkStateHessians(stateHessians, dimMeas, dimState);
            
            obj.checkNoiseJacobian(noiseJacobian, dimMeas, dimNoise);
            obj.checkNoiseHessians(noiseHessians, dimMeas, dimNoise);
            
            [stateHessMean, ...
             stateHessCov, ...
             noiseHessMean, ...
             noiseHessCov] = obj.getHessianMomentsStateAndNoise(dimState, stateHessians, stateCov, ...
                                                                dimNoise, noiseHessians, noiseCov, dimMeas);
            
            h = measModel.measurementEquation(stateMean, noiseMean) + ...
                stateHessMean + noiseHessMean;
        end
    end
    
    methods (Access = 'private')
        function [stateHessMean, ...
                  stateHessCov] = getHessianMomentsState(~, dimState, stateHessians, stateCov, dimOutput)
            stateHessProd = nan(dimState, dimState, dimState);
            
            stateHessMean = nan(dimOutput, 1);
            
            for i = 1:dimOutput
                stateHessProd(:, :, i) = stateHessians(:, :, i) * stateCov;
                stateHessMean(i)       = trace(stateHessProd(:, :, i));
            end
            
            stateHessCov = nan(dimOutput, dimOutput);
            
            for i = 1:dimOutput
                mat                = stateHessProd(:, :, i) .* stateHessProd(:, :, i)';
                stateHessCov(i, i) = sum(mat(:));        % = trace(stateHessProd(:, :, i)^2)
                
                for j = (i + 1):dimOutput
                    mat                = stateHessProd(:, :, i) .* stateHessProd(:, :, j)';
                    stateHessCov(i, j) = sum(mat(:));    % = trace(stateHessProd(:, :, i) * stateHessProd(:, :, j))
                    stateHessCov(j, i) = stateHessCov(i, j);
                end
            end
            
            stateHessMean = 0.5 * stateHessMean;
            stateHessCov  = 0.5 * stateHessCov;
        end
        
        function [stateHessMean, ...
                  stateHessCov, ...
                  noiseHessMean, ...
                  noiseHessCov] = getHessianMomentsStateAndNoise(~, dimState, stateHessians, stateCov, ...
                                                                 dimNoise, noiseHessians, noiseCov, dimOutput)
            stateHessProd = nan(dimState, dimState, dimState);
            noiseHessProd = nan(dimNoise, dimNoise, dimState);
            
            stateHessMean = nan(dimOutput, 1);
            noiseHessMean = nan(dimOutput, 1);
            
            for i = 1:dimOutput
                stateHessProd(:, :, i) = stateHessians(:, :, i) * stateCov;
                stateHessMean(i)       = trace(stateHessProd(:, :, i));
                
                noiseHessProd(:, :, i) = noiseHessians(:, :, i) * noiseCov;
                noiseHessMean(i)       = trace(noiseHessProd(:, :, i));
            end
            
            stateHessCov = nan(dimOutput, dimOutput);
            noiseHessCov = nan(dimOutput, dimOutput);
            
            for i = 1:dimOutput
                mat                = stateHessProd(:, :, i) .* stateHessProd(:, :, i)';
                stateHessCov(i, i) = sum(mat(:));           % = trace(stateHessProd(:, :, i)^2)
                
                mat                = noiseHessProd(:, :, i) .* noiseHessProd(:, :, i)';
                noiseHessCov(i, i) = sum(mat(:));           % = trace(noiseHessProd(:, :, i)^2)
                
                for j = (i + 1):dimOutput
                    mat                = stateHessProd(:, :, i) .* stateHessProd(:, :, j)';
                    stateHessCov(i, j) = sum(mat(:));       % = trace(stateHessProd(:, :, i) * stateHessProd(:, :, j))
                    stateHessCov(j, i) = stateHessCov(i, j);
                    
                    mat                = noiseHessProd(:, :, i) .* noiseHessProd(:, :, j)';
                    noiseHessCov(i, j) = sum(mat(:));       % = trace(noiseHessProd(:, :, i) * noiseHessProd(:, :, j))
                    noiseHessCov(j, i) = noiseHessCov(i, j);
                end
            end
            
            stateHessMean = 0.5 * stateHessMean;
            stateHessCov  = 0.5 * stateHessCov;
            
            noiseHessMean = 0.5 * noiseHessMean;
            noiseHessCov  = 0.5 * noiseHessCov;
        end
    end
end

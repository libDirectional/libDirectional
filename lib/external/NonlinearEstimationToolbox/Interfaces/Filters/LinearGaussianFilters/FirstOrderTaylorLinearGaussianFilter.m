
classdef FirstOrderTaylorLinearGaussianFilter < LinearGaussianFilter
    % Abstract base class for linear Gaussian filters that are based on the
    % first-order Taylor series approximation.
    %
    % FirstOrderTaylorLinearGaussianFilter Methods:
    %   FirstOrderTaylorLinearGaussianFilter - Class constructor.
    %   copy                                 - Copy a Filter instance.
    %   copyWithName                         - Copy a Filter instance and give the copy a new name/description.
    %   getName                              - Get the filter name/description.
    %   setColor                             - Set the filter color/plotting properties.
    %   getColor                             - Get the filter color/plotting properties.
    %   setState                             - Set the system state.
    %   getState                             - Get the system state.
    %   setStateMeanAndCov                   - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov                   - Get mean and covariance matrix of the system state.
    %   getStateDim                          - Get the dimension of the system state.
    %   predict                              - Perform a state prediction.
    %   update                               - Perform a measurement update.
    %   step                                 - Perform a combined state prediction and measurement update.
    %   setStateDecompDim                    - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                    - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing          - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing          - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing              - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing              - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold               - Set the measurement gating threshold.
    %   getMeasGatingThreshold               - Get the measurement gating threshold.
    
    % Literature:
    %   Dan Simon,
    %   Optimal State Estimation,
    %   Sections 13.2,
    %   1st ed. Wiley & Sons, 2006.
    
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
        function obj = FirstOrderTaylorLinearGaussianFilter(name)
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
            %   << obj (FirstOrderTaylorLinearGaussianFilter)
            %      A new FirstOrderTaylorLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@LinearGaussianFilter(name);
        end
    end
    
    methods (Sealed, Access = 'protected')
        % State prediction moment computation
        function [predictedStateMean, ...
                  predictedStateCov] = predictSysModel(obj, sysModel)
            [noiseMean, ~, noiseCovSqrt] = sysModel.noise.getMeanAndCov();
            dimNoise = size(noiseMean, 1);
            
            % Linearize system model around current state mean and noise mean
            [stateJacobian, ...
             noiseJacobian] = sysModel.derivative(obj.stateMean, noiseMean);
            
            % Check computed derivatives
            obj.checkStateJacobian(stateJacobian, obj.dimState, obj.dimState);
            obj.checkNoiseJacobian(noiseJacobian, obj.dimState, dimNoise);
            
            % Compute predicted state mean
            predictedStateMean = sysModel.systemEquation(obj.stateMean, noiseMean);
            
            % Compute predicted state covariance
            A = stateJacobian * obj.stateCovSqrt;
            B = noiseJacobian * noiseCovSqrt;
            
            predictedStateCov = A * A' + B * B';
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictAddNoiseSysModel(obj, sysModel)
            [noiseMean, noiseCov] = sysModel.noise.getMeanAndCov();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimNoise);
            
            % Linearize system model around current state mean
            stateJacobian = sysModel.derivative(obj.stateMean);
            
            % Check computed derivative
            obj.checkStateJacobian(stateJacobian, obj.dimState, obj.dimState);
            
            % Compute predicted state mean
            predictedStateMean = sysModel.systemEquation(obj.stateMean) + noiseMean;
            
            % Compute predicted state covariance
            A = stateJacobian * obj.stateCovSqrt;
            
            predictedStateCov = A * A' + noiseCov;
        end
        
        % Helpers for the first-order Taylor-based linear measurement updates
        function [h, stateJacobian] = evaluateAddNoiseMeasModel(obj, measModel, dimMeas, stateMean)
            dimState = size(stateMean, 1);
            
            % Linearize measurement model around current state mean
            stateJacobian = measModel.derivative(stateMean);
            
            % Check computed derivative
            obj.checkStateJacobian(stateJacobian, dimMeas, dimState);
            
            h = measModel.measurementEquation(stateMean);
        end
        
        function [h, stateJacobian, ...
                  noiseJacobian] = evaluateMeasModel(obj, measModel, dimMeas, dimNoise, ...
                                                     noiseMean, stateMean)
            dimState = size(stateMean, 1);
            
            % Linearize measurement model around current state mean and noise mean
            [stateJacobian, ...
             noiseJacobian] = measModel.derivative(stateMean, noiseMean);
            
            % Check computed derivatives
            obj.checkStateJacobian(stateJacobian, dimMeas, dimState);
            obj.checkNoiseJacobian(noiseJacobian, dimMeas, dimNoise);
            
            h = measModel.measurementEquation(stateMean, noiseMean);
        end
    end
end

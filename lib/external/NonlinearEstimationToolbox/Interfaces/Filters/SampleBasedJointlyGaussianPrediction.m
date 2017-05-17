
classdef SampleBasedJointlyGaussianPrediction < GaussianFilter
    % Abstract base class for Gaussian filters that use a sample-based prediction,
    % where the joint density of state and noise is assumed to be Gaussian.
    %
    % SampleBasedJointlyGaussianPrediction Methods:
    %   SampleBasedJointlyGaussianPrediction - Class constructor.
    %   copy                                 - Copy a Filter instance.
    %   copyWithName                         - Copy a Filter instance and give the copy a new name / description.
    %   getName                              - Get the filter name / description.
    %   setColor                             - Set the filter color / plotting properties.
    %   getColor                             - Get the current filter color / plotting properties.
    %   setState                             - Set the system state.
    %   getState                             - Get the current system state.
    %   getStateDim                          - Get the dimension of the current system state.
    %   predict                              - Perform a time update (prediction step).
    %   update                               - Perform a measurement update (filter step) using the given measurement(s).
    %   step                                 - Perform a combined time and measurement update.
    %   getPointEstimate                     - Get a point estimate of the current system state.
    %   setUseAnalyticSystemModel            - Enable or disable the use of analytic moment calculation during a prediction.
    %   getUseAnalyticSystemModel            - Get the current use of analytic moment calculation during a prediction.
    %   setStateDecompDim                    - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                    - Get the dimension of the unobservable part of the system state.
    
    % >> This function/class is part of the Nonlinear Estimation Toolbox
    %
    %    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
    %
    %    Copyright (C) 2016  Jannik Steinbring <jannik.steinbring@kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        http://isas.uka.de
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
        function obj = SampleBasedJointlyGaussianPrediction(name, samplingPrediction)
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
            %   >> samplingPrediction (Subclass of GaussianSampling)
            %      The Gaussian sampling method used to perform a state prediction.
            %
            % Returns:
            %   << obj (SampleBasedJointlyGaussianPrediction)
            %      A new SampleBasedJointlyGaussianPrediction instance.
            
            % Call superclass constructor
            obj = obj@GaussianFilter(name);
            
            if ~Checks.isClass(samplingPrediction, 'GaussianSampling')
                obj.error('InvalidGaussianSampling', ...
                          'samplingPrediction must be a subclass of GaussianSampling.');
            end
            
            obj.samplingPrediction = samplingPrediction;
        end
    end
    
    methods (Access = 'protected')
        function [predictedStateMean, ...
                  predictedStateCov] = predictedMomentsArbitraryNoise(obj, sysModel)
            [noiseMean, ~, noiseCovSqrt] = sysModel.noise.getMeanAndCovariance();
            
            % Generate state and noise samples
            [stateSamples, ...
             noiseSamples, ...
             weights, ...
             numSamples] = Utils.getStateNoiseSamples(obj.samplingPrediction, ...
                                                      obj.stateMean, obj.stateCovSqrt, ...
                                                      noiseMean, noiseCovSqrt);
            
            % Propagate samples through system equation
            predictedStates = sysModel.systemEquation(stateSamples, noiseSamples);
            
            % Check predicted state samples
            obj.checkPredictedStateSamples(predictedStates, numSamples);
            
            % Compute predicted state mean and covariance
            [predictedStateMean, ...
             predictedStateCov] = Utils.getMeanAndCov(predictedStates, weights);
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictedMomentsAdditiveNoise(obj, sysModel)
            [noiseMean, noiseCov] = sysModel.noise.getMeanAndCovariance();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimNoise);
            
            % Generate state samples
            [stateSamples, ...
             weights, ...
             numSamples] = Utils.getStateSamples(obj.samplingPrediction, ...
                                                 obj.stateMean, obj.stateCovSqrt);
            
            % Propagate samples through deterministic system equation
            predictedStates = sysModel.systemEquation(stateSamples);
            
            % Check predicted state samples
            obj.checkPredictedStateSamples(predictedStates, numSamples);
            
            [mean, cov] = Utils.getMeanAndCov(predictedStates, weights);
            
            % Compute predicted state mean
            predictedStateMean = mean + noiseMean;
            
            % Compute predicted state covariance
            predictedStateCov = cov + noiseCov;
        end
    end
    
    methods (Access = 'protected')
        function cpObj = copyElement(obj)
            cpObj = obj.copyElement@GaussianFilter();
            
            cpObj.samplingPrediction = obj.samplingPrediction.copy();
        end
    end
    
    properties (SetAccess = 'private', GetAccess = 'protected')
        % Gaussian sampling technique used for prediction.
        samplingPrediction;
    end
end


classdef SampleBasedGaussianFilter < GaussianFilter
    % Abstract base class for sample-based Gaussian filters.
    %
    % This type of filter is a special case of a Gaussian filter where prediction and filter
    % step are computed with the aid of samples (or particles).
    %
    % SampleBasedGaussianFilter Methods:
    %   SampleBasedGaussianFilter - Class constructor.
    %   getName                   - Get the filter name / description.
    %   setColor                  - Set the filter color / plotting properties.
    %   getColor                  - Get the current filter color / plotting properties.
    %   setState                  - Set the system state.
    %   getState                  - Get the current system state.
    %   getStateDim               - Get the dimension of the current system state.
    %   predict                   - Perform a time update (prediction step).
    %   update                    - Perform a measurement update (filter step) using the given measurement(s).
    %   step                      - Perform a combined time and measurement update.
    %   getPointEstimate          - Get a point estimate of the current system state.
    %   setUseAnalyticSystemModel - Enable or disable the use of analytic moment calculation during a prediction.
    %   getUseAnalyticSystemModel - Get the current use of analytic moment calculation during a prediction.
    
    % >> This function/class is part of the Nonlinear Estimation Toolbox
    %
    %    For more information, see https://bitbucket.org/nonlinearestimation/toolbox
    %
    %    Copyright (C) 2015  Jannik Steinbring <jannik.steinbring@kit.edu>
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
        function obj = SampleBasedGaussianFilter(name)
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
            %   << obj (SampleBasedGaussianFilter)
            %      A new SampleBasedGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@GaussianFilter(name);
            
            obj.setUseAnalyticSystemModel(false);
        end
        
        function setUseAnalyticSystemModel(obj, useAnalyticSysModel)
            % Enable or disable the use of analytic moment calculation during a prediction.
            %
            % If true, analytic moment calculation will be used for the
            % prediction if the given system model supports it (i.e., if a
            % system model is given, which inherits from
            % AnalyticSystemModel). Otherwise, a sample-based prediction
            % will be used (provided a system model supported by this filter
            % is given).
            %
            % By default, the use of analytic system models is disabled.
            %
            % Parameters:
            %   >> useAnalyticSysModel (Logical scalar)
            %      If true, analytic moment calculation will be used during
            %      the prediction. Otherwise, a sample-based prediction
            %      will be performed.
            
            if ~Checks.isFlag(useAnalyticSysModel)
                obj.error('InvalidFlag', ...
                          'useAnalyticSysModel must be a logical scalar.');
            end
            
            obj.useAnalyticSysModel = useAnalyticSysModel;
        end
        
        function useAnalyticSysModel = getUseAnalyticSystemModel(obj)
            % Get the current use of analytic moment calculation during a prediction.
            %
            % Returns:
            %   << useAnalyticSysModel (Logical scalar)
            %      If true, analytic moment calculation will be used during
            %      the prediction. Otherwise, a sample-based prediction
            %      will be performed.
            
            useAnalyticSysModel = obj.useAnalyticSysModel;
        end
    end
    
    methods (Abstract, Access = 'protected')
        predictArbitraryNoise(obj, sysModel);
        
        predictAdditiveNoise(obj, sysModel);
        
        predictMixedNoise(obj, sysModel);
    end
    
    methods (Access = 'protected')
        function performPrediction(obj, sysModel)
            if obj.useAnalyticSysModel && ...
               Checks.isClass(sysModel, 'AnalyticSystemModel')
                obj.predictAnalytic(sysModel);
            elseif Checks.isClass(sysModel, 'LinearSystemModel')
                obj.predictAnalytic(sysModel);
            elseif Checks.isClass(sysModel, 'SystemModel')
                obj.predictArbitraryNoise(sysModel);
            elseif Checks.isClass(sysModel, 'AdditiveNoiseSystemModel')
                obj.predictAdditiveNoise(sysModel);
            elseif Checks.isClass(sysModel, 'MixedNoiseSystemModel')
                obj.predictMixedNoise(sysModel);
            else
                obj.errorSysModel('AnalyticSystemModel (if enabled, see setUseAnalyticSystemModel())', ...
                                  'LinearSystemModel', ...
                                  'SystemModel', ...
                                  'AdditiveNoiseSystemModel', ...
                                  'MixedNoiseSystemModel');
            end
        end
        
        function predictJointGaussianArbitraryNoise(obj, sysModel, sampling)
            [predictedStateMean, ...
             predictedStateCov] = obj.predictArbitraryNoiseMoments(sysModel, sampling);
            
            obj.checkAndSavePrediction(predictedStateMean, predictedStateCov);
        end
        
        function predictJointGaussianAdditiveNoise(obj, sysModel, sampling)
            [noiseMean, noiseCov] = sysModel.noise.getMeanAndCovariance();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimNoise);
            
            % Generate state samples
            [stateSamples, ...
             weights, ...
             numSamples] = Utils.getStateSamples(sampling, ...
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
            
            obj.checkAndSavePrediction(predictedStateMean, predictedStateCov);
        end
        
        function predictJointGaussianMixedNoise(obj, sysModel, sampling)
            [addNoiseMean, addNoiseCov]  = sysModel.additiveNoise.getMeanAndCovariance();
            dimAddNoise = size(addNoiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimAddNoise);
            
            [mean, cov] = obj.predictArbitraryNoiseMoments(sysModel, sampling);
            
            % Compute predicted state mean
            predictedStateMean = mean + addNoiseMean;
            
            % Compute predicted state covariance
            predictedStateCov = cov + addNoiseCov;
            
            obj.checkAndSavePrediction(predictedStateMean, predictedStateCov);
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictArbitraryNoiseMoments(obj, sysModel, sampling)
            [noiseMean, ~, noiseCovSqrt] = sysModel.noise.getMeanAndCovariance();
            
            % Generate state and noise samples
            [stateSamples, ...
             noiseSamples, ...
             weights, ...
             numSamples] = Utils.getStateNoiseSamples(sampling, ...
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
    end
    
    properties (Access = 'private')
        % Flag that indicates the use of an AnalyticSystemMdoel.
        useAnalyticSysModel;
    end
end

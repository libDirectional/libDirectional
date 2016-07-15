
classdef LRKF < KF & SampleBasedGaussianFilter
    % Abstract base class for Linear Regression Kalman Filters (LRKFs)
    %
    % This type of filter implements a (nonlinear) Kalman filter with the aid of Gaussian
    % sampling techniques.
    %
    % LRKF Methods:
    %   LRKF                           - Class constructor.
    %   getName                        - Get the filter name / description.
    %   setColor                       - Set the filter color / plotting properties.
    %   getColor                       - Get the current filter color / plotting properties.
    %   setState                       - Set the system state.
    %   getState                       - Get the current system state.
    %   getStateDim                    - Get the dimension of the current system state.
    %   predict                        - Perform a time update (prediction step).
    %   update                         - Perform a measurement update (filter step) using the given measurement(s).
    %   step                           - Perform a combined time and measurement update.
    %   getPointEstimate               - Get a point estimate of the current system state.
    %   setMaxNumIterations            - Set the maximum number of iterations that will be performed during a measurement update.
    %   getMaxNumIterations            - Get the current maximum number of iterations that will be performed during a measurement update.
    %   setMeasValidationThreshold     - Set a threshold to perform a measurement validation (measurement acceptance/rejection).
    %   getMeasValidationThreshold     - Get the current measurement validation threshold.
    %   getLastUpdateData              - Get information from the last performed measurement update.
    %   setUseAnalyticSystemModel      - Enable or disable the use of analytic moment calculation during a prediction.
    %   getUseAnalyticSystemModel      - Get the current use of analytic moment calculation during a prediction.
    %   setUseAnalyticMeasurementModel - Enable or disable the use of analytic moment calculation during a filter step.
    %   getUseAnalyticMeasurementModel - Get the current use of analytic moment calculation during a filter step.
    
    % Literature:
    %   Jannik Steinbring and Uwe D. Hanebeck,
    %   LRKF Revisited: The Smart Sampling Kalman Filter (S²KF),
    %   Journal of Advances in Information Fusion, Vol. 9, No. 2, Dec 2014, pp. 106-123.
    %
    %   T. Lefebvre and H. Bruyninckx and J. De Schutter,
    %   Appendix A: The Linear Regression Kalman Filter,
    %   Nonlinear Kalman Filtering for Force-Controlled Robot Tasks
    %   vol. 19, Springer, 2005
    %
    %   T. Lefebvre and H. Bruyninckx and J. De Schutter,
    %   Kalman filters for non-linear systems: a comparison of performance,
    %   International Journal of Control vol. 77, pages 639-653, 2004
    
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
        function obj = LRKF(name, samplingPrediction, samplingUpdate)
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
            %      The Gaussian sampling method used by the LRKF to perform a state prediction.
            %
            %   >> samplingUpdate (Subclass of GaussianSampling)
            %      The Gaussian sampling method used by the LRKF to perform a measurement update.
            %      Default: the sampling for the prediction is used for the update, too.
            %
            % Returns:
            %   << obj (LRKF)
            %      A new LRKF instance.
            
            % Call superclass constructors
            obj = obj@KF(name);
            obj = obj@SampleBasedGaussianFilter(name);
            
            if ~Checks.isClass(samplingPrediction, 'GaussianSampling')
                obj.error('InvalidGaussianSampling', ...
                          'samplingPrediction must be a subclass of GaussianSampling.');
            end
            
            obj.samplingPrediction = samplingPrediction;
            
            if nargin == 3
                if ~Checks.isClass(samplingUpdate, 'GaussianSampling')
                    obj.error('InvalidGaussianSampling', ...
                              'samplingUpdate must be a subclass of GaussianSampling.');
                end
                
                obj.samplingUpdate = samplingUpdate;
            else
                % Use the same samplings for both prediction and update
                obj.samplingUpdate = samplingPrediction;
            end
            
            obj.setUseAnalyticMeasurementModel(false);
        end
        
        function setUseAnalyticMeasurementModel(obj, useAnalyticMeasModel)
            % Enable or disable the use of analytic moment calculation during a filter step.
            %
            % If true, analytic moment calculation will be used for the
            % measurement update if the given measurement model supports it,
            % (i.e., if a measurement model is given, which inherits from
            % AnalyticMeasurementModel). Otherwise, a sample-based
            % measurement update will be used (provided a measurement model
            % supported by this filter is given).
            %
            % By default, the use of analytic measurement models is disabled.
            %
            % Parameters:
            %   >> useAnalyticMeasModel (Logical scalar)
            %      If true, analytic moment calculation will be used during a
            %      measurement update. Otherwise, a sample-based measurement
            %      update will be performed.
            
            if ~Checks.isFlag(useAnalyticMeasModel)
                obj.error('InvalidFlag', ...
                          'useAnalyticMeasModel must be a logical scalar.');
            end
            
            obj.useAnalyticMeasModel = useAnalyticMeasModel;
        end
        
        function useAnalyticMeasModel = getUseAnalyticMeasurementModel(obj)
            % Get the current use of analytic moment calculation during a filter step.
            %
            % Returns:
            %   << useAnalyticMeasModel (Logical scalar)
            %      If true, analytic moment calculation will be used during a
            %      measurement update. Otherwise, a sample-based measurement
            %      update will be performed.
            
            useAnalyticMeasModel = obj.useAnalyticMeasModel;
        end
    end
    
    methods (Access = 'protected')
        function predictArbitraryNoise(obj, sysModel)
            obj.predictJointGaussianArbitraryNoise(sysModel, obj.samplingPrediction);
        end
        
        function predictAdditiveNoise(obj, sysModel)
            obj.predictJointGaussianAdditiveNoise(sysModel, obj.samplingPrediction);
        end
        
        function predictMixedNoise(obj, sysModel)
            obj.predictJointGaussianMixedNoise(sysModel, obj.samplingPrediction);
        end
        
        function performUpdate(obj, measModel, measurements)
            if obj.useAnalyticMeasModel && ...
               Checks.isClass(measModel, 'AnalyticMeasurementModel')
                obj.updateAnalytic(measModel, measurements);
            elseif Checks.isClass(measModel, 'LinearMeasurementModel')
                obj.updateAnalytic(measModel, measurements);
            elseif Checks.isClass(measModel, 'MeasurementModel')
                obj.updateArbitraryNoise(measModel, measurements);
            elseif Checks.isClass(measModel, 'AdditiveNoiseMeasurementModel')
                obj.updateAdditiveNoise(measModel, measurements);
            elseif Checks.isClass(measModel, 'MixedNoiseMeasurementModel')
                obj.updateMixedNoise(measModel, measurements);
            else
                obj.errorMeasModel('AnalyticMeasurementModel (if enabled, see setUseAnalyticMeasurementModel())', ...
                                   'LinearMeasurementModel', ...
                                   'MeasurementModel', ...
                                   'AdditiveNoiseMeausrementModel', ...
                                   'MixedNoiseMeasurementModel');
            end
        end
        
        function updateArbitraryNoise(obj, measModel, measurements)
            [dimMeas, numMeas]           = size(measurements);
            [noiseMean, ~, noiseCovSqrt] = measModel.noise.getMeanAndCovariance();
            dimNoise     = size(noiseMean, 1);
            noiseMean    = repmat(noiseMean, numMeas, 1);
            noiseCovSqrt = Utils.blockDiag(noiseCovSqrt, numMeas);
            
            % Perform state update
            obj.kalmanUpdate(measModel, measurements, @obj.momentFuncArbitraryNoise, ...
                             dimNoise, dimMeas, numMeas, noiseMean, noiseCovSqrt);
        end
        
        function updateAdditiveNoise(obj, measModel, measurements)
            [dimMeas, numMeas]    = size(measurements);
            [noiseMean, noiseCov] = measModel.noise.getMeanAndCovariance();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimNoise);
            
            % Perform state update
            obj.kalmanUpdate(measModel, measurements, @obj.momentFuncAdditiveNoise, ...
                             dimMeas, numMeas, noiseMean, noiseCov);
        end
        
        function updateMixedNoise(obj, measModel, measurements)
            [dimMeas, numMeas]           = size(measurements);
            [addNoiseMean, addNoiseCov]  = measModel.additiveNoise.getMeanAndCovariance();
            [noiseMean, ~, noiseCovSqrt] = measModel.noise.getMeanAndCovariance();
            dimAddNoise  = size(addNoiseMean, 1);
            dimNoise     = size(noiseMean, 1);
            noiseMean    = repmat(noiseMean, numMeas, 1);
            noiseCovSqrt = Utils.blockDiag(noiseCovSqrt, numMeas);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimAddNoise);
            
            % Perform state update
            obj.kalmanUpdate(measModel, measurements, @obj.momentFuncMixedNoise, ...
                             dimNoise, dimMeas, numMeas, addNoiseMean, addNoiseCov, noiseMean, noiseCovSqrt);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncArbitraryNoise(obj, measModel, iterNum, iterStateMean, iterStateCov, iterStateCovSqrt, ...
                                                                dimNoise, dimMeas, numMeas, noiseMean, noiseCovSqrt)
            % Generate state and noise samples
            [stateSamples, ...
             noiseSamples, ...
             weights, ...
             numSamples] = Utils.getStateNoiseSamples(obj.samplingUpdate, ...
                                                      iterStateMean, iterStateCovSqrt, ...
                                                      noiseMean, noiseCovSqrt);
            
            measSamples = nan(dimMeas * numMeas, numSamples);
            a = 1; c = 1;
            
            for i = 1:numMeas
                b = i * dimMeas;
                d = i * dimNoise;
                
                % Propagate samples through measurement equation
                meas = measModel.measurementEquation(stateSamples, noiseSamples(c:d, :));
                
                % Check computed measurements
                obj.checkComputedMeasurements(meas, dimMeas, numSamples);
                
                measSamples(a:b, :) = meas;
                
                a = b + 1;
                c = d + 1;
            end
            
            [measMean, measCov, ...
             stateMeasCrossCov] = Utils.getMeanCovAndCrossCov(iterStateMean, stateSamples, ...
                                                              measSamples, weights);
            
            if iterNum > 1
                [measMean, measCov, ...
                 stateMeasCrossCov] = obj.momentCorrection(iterStateMean, iterStateCov, ...
                                                           measMean, measCov, stateMeasCrossCov);
            end
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncAdditiveNoise(obj, measModel, iterNum, iterStateMean, iterStateCov, iterStateCovSqrt, ...
                                                               dimMeas, numMeas, noiseMean, noiseCov)
            % Generate state samples
            [stateSamples, ...
             weights, ...
             numSamples] = Utils.getStateSamples(obj.samplingUpdate, ...
                                                 iterStateMean, iterStateCovSqrt);
            
            % Propagate samples through deterministic measurement equation
            measSamples = measModel.measurementEquation(stateSamples);
            
            % Check computed measurements
            obj.checkComputedMeasurements(measSamples, dimMeas, numSamples);
            
            [mean, cov, crossCov] = Utils.getMeanCovAndCrossCov(iterStateMean, stateSamples, ...
                                                                measSamples, weights);
            
            % Compute measurement mean
            measMean = repmat(mean + noiseMean, numMeas, 1);
            
            % Compute measurement covariance
            measCov = Utils.baseBlockDiag(cov, noiseCov, numMeas);
            
            % Compute state measurement cross-covariance
            stateMeasCrossCov = repmat(crossCov, 1, numMeas);
            
            if iterNum > 1
                [measMean, measCov, ...
                 stateMeasCrossCov] = obj.momentCorrection(iterStateMean, iterStateCov, ...
                                                           measMean, measCov, stateMeasCrossCov);
            end
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncMixedNoise(obj, measModel, iterNum, iterStateMean, iterStateCov, iterStateCovSqrt, ...
                                                            dimNoise, dimMeas, numMeas, addNoiseMean, addNoiseCov, noiseMean, noiseCovSqrt)
            % Generate state and noise samples
            [stateSamples, ...
             noiseSamples, ...
             weights, ...
             numSamples] = Utils.getStateNoiseSamples(obj.samplingUpdate, ...
                                                      iterStateMean, iterStateCovSqrt, ...
                                                      noiseMean, noiseCovSqrt);
            
            measSamples = nan(dimMeas * numMeas, numSamples);
            a = 1; c = 1;
            
            for i = 1:numMeas
                b = i * dimMeas;
                d = i * dimNoise;
                
                % Propagate samples through measurement equation
                meas = measModel.measurementEquation(stateSamples, noiseSamples(c:d, :));
                
                % Check computed measurements
                obj.checkComputedMeasurements(meas, dimMeas, numSamples);
                
                measSamples(a:b, :) = meas;
                
                a = b + 1;
                c = d + 1;
            end
            
            [measMean, measCov, ...
             stateMeasCrossCov] = Utils.getMeanCovAndCrossCov(iterStateMean, stateSamples, ...
                                                              measSamples, weights);
            
            % Compute measurement mean
            measMean = measMean + repmat(addNoiseMean, numMeas, 1);
            
            % Compute measurement covariance
            measCov = measCov + Utils.blockDiag(addNoiseCov, numMeas);
            
            if iterNum > 1
                [measMean, measCov, ...
                 stateMeasCrossCov] = obj.momentCorrection(iterStateMean, iterStateCov, ...
                                                           measMean, measCov, stateMeasCrossCov);
            end
        end
    end
    
    properties (Access = 'private')
        % GaussianSampling used by the LRKF for the prediction step.
        samplingPrediction;
        
        % GaussianSampling used by the LRKF for the filter step.
        samplingUpdate;
        
        % Flag that indicates the use of an AnalyticMeasurementMdoel.
        useAnalyticMeasModel;
    end
end

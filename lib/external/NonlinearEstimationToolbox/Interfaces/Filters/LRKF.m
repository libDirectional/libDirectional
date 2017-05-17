
classdef LRKF < KF & SampleBasedJointlyGaussianPrediction
    % Abstract base class for Linear Regression Kalman Filters (LRKFs)
    %
    % This type of filter implements a (nonlinear) Kalman filter with the aid of Gaussian
    % sampling techniques.
    %
    % LRKF Methods:
    %   LRKF                           - Class constructor.
    %   copy                           - Copy a Filter instance.
    %   copyWithName                   - Copy a Filter instance and give the copy a new name / description.
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
    %   setUseAnalyticSystemModel      - Enable or disable the use of analytic moment calculation during a prediction.
    %   getUseAnalyticSystemModel      - Get the current use of analytic moment calculation during a prediction.
    %   setStateDecompDim              - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim              - Get the dimension of the unobservable part of the system state.
    %   setUseAnalyticMeasurementModel - Enable or disable the use of analytic moment calculation during a filter step.
    %   getUseAnalyticMeasurementModel - Get the current use of analytic moment calculation during a filter step.
    %   setMaxNumIterations            - Set the maximum number of iterations that will be performed during a measurement update.
    %   getMaxNumIterations            - Get the current maximum number of iterations that will be performed during a measurement update.
    %   setMeasValidationThreshold     - Set a threshold to perform a measurement validation (measurement acceptance/rejection).
    %   getMeasValidationThreshold     - Get the current measurement validation threshold.
    %   getLastUpdateData              - Get information from the last performed measurement update.
    
    % Literature:
    %   Ángel F. Garcı́a-Fernández, Lennart Svensson, Mark Morelande, and Simo Särkkä,
    %   Posterior Linearisation Filter: Principles and Implementation Using Sigma Points,
    %   IEEE Transactions on Signal Processing, Vol. 63, No. 20, Oct 2015, pp. 5561–5573.
    %
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
            obj = obj@SampleBasedJointlyGaussianPrediction(name, samplingPrediction);
            
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
        end
    end
    
    methods (Access = 'protected')
        function momentFunc = getMomentFuncArbitraryNoise(obj, measModel, measurements)
            [dimMeas, numMeas]           = size(measurements);
            [noiseMean, ~, noiseCovSqrt] = measModel.noise.getMeanAndCovariance();
            dimNoise     = size(noiseMean, 1);
            noiseMean    = repmat(noiseMean, numMeas, 1);
            noiseCovSqrt = Utils.blockDiag(noiseCovSqrt, numMeas);
            
            momentFunc = @(priorMean, priorCov, priorCovSqrt, iterNum, iterMean, iterCov, iterCovSqrt) ...
                         obj.momentFuncArbitraryNoise(priorMean, priorCov, priorCovSqrt, ...
                                                      iterNum, iterMean, iterCov, iterCovSqrt, ...
                                                      measModel, dimNoise, dimMeas, numMeas, noiseMean, noiseCovSqrt);
        end
        
        function momentFunc = getMomentFuncAdditiveNoise(obj, measModel, measurements)
            [dimMeas, numMeas]    = size(measurements);
            [noiseMean, noiseCov] = measModel.noise.getMeanAndCovariance();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimNoise);
            
            momentFunc = @(priorMean, priorCov, priorCovSqrt, iterNum, iterMean, iterCov, iterCovSqrt) ...
                         obj.momentFuncAdditiveNoise(priorMean, priorCov, priorCovSqrt, ...
                                                     iterNum, iterMean, iterCov, iterCovSqrt, ...
                                                     measModel, dimMeas, numMeas, noiseMean, noiseCov);
        end
        
        function momentFunc = getMomentFuncMixedNoise(obj, measModel, measurements)
            [dimMeas, numMeas]           = size(measurements);
            [addNoiseMean, addNoiseCov]  = measModel.additiveNoise.getMeanAndCovariance();
            [noiseMean, ~, noiseCovSqrt] = measModel.noise.getMeanAndCovariance();
            dimAddNoise  = size(addNoiseMean, 1);
            dimNoise     = size(noiseMean, 1);
            noiseMean    = repmat(noiseMean, numMeas, 1);
            noiseCovSqrt = Utils.blockDiag(noiseCovSqrt, numMeas);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimAddNoise);
            
            momentFunc = @(priorMean, priorCov, priorCovSqrt, iterNum, iterMean, iterCov, iterCovSqrt) ...
                         obj.momentFuncMixedNoise(priorMean, priorCov, priorCovSqrt, ...
                                                  iterNum, iterMean, iterCov, iterCovSqrt, ...
                                                  measModel, dimNoise, dimMeas, numMeas, addNoiseMean, ...
                                                  addNoiseCov, noiseMean, noiseCovSqrt);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncArbitraryNoise(obj, priorMean, priorCov, priorCovSqrt, ...
                                                                iterNum, iterMean, iterCov, iterCovSqrt, ...
                                                                measModel, dimNoise, dimMeas, numMeas, noiseMean, noiseCovSqrt)
            % Generate state and noise samples
            [stateSamples, ...
             noiseSamples, ...
             weights, ...
             numSamples] = Utils.getStateNoiseSamples(obj.samplingUpdate, ...
                                                      iterMean, iterCovSqrt, ...
                                                      noiseMean, noiseCovSqrt);
            
            measSamples = obj.getStackedMeasSamples(measModel, stateSamples, noiseSamples, ...
                                                    numSamples, dimMeas, numMeas, dimNoise);
            
            [measMean, measCov, ...
             stateMeasCrossCov] = Utils.getMeanCovAndCrossCov(iterMean, stateSamples, ...
                                                              measSamples, weights);
            
            if iterNum > 1
                [measMean, measCov, ...
                 stateMeasCrossCov] = KF.momentCorrection(priorMean, priorCov, priorCovSqrt, ...
                                                          iterMean, iterCov, iterCovSqrt, ...
                                                          measMean, measCov, stateMeasCrossCov);
            end
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncAdditiveNoise(obj, priorMean, priorCov, priorCovSqrt, ...
                                                               iterNum, iterMean, iterCov, iterCovSqrt, ...
                                                               measModel, dimMeas, numMeas, noiseMean, noiseCov)
            % Generate state samples
            [stateSamples, ...
             weights, ...
             numSamples] = Utils.getStateSamples(obj.samplingUpdate, ...
                                                 iterMean, iterCovSqrt);
            
            % Propagate samples through deterministic measurement equation
            measSamples = measModel.measurementEquation(stateSamples);
            
            % Check computed measurements
            obj.checkComputedMeasurements(measSamples, dimMeas, numSamples);
            
            [mean, cov, crossCov] = Utils.getMeanCovAndCrossCov(iterMean, stateSamples, ...
                                                                measSamples, weights);
            
            % Compute measurement mean
            measMean = repmat(mean + noiseMean, numMeas, 1);
            
            % Compute measurement covariance
            measCov = Utils.baseBlockDiag(cov, noiseCov, numMeas);
            
            % Compute state measurement cross-covariance
            stateMeasCrossCov = repmat(crossCov, 1, numMeas);
            
            if iterNum > 1
                [measMean, measCov, ...
                 stateMeasCrossCov] = KF.momentCorrection(priorMean, priorCov, priorCovSqrt, ...
                                                          iterMean, iterCov, iterCovSqrt, ...
                                                          measMean, measCov, stateMeasCrossCov);
            end
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov] = momentFuncMixedNoise(obj, priorMean, priorCov, priorCovSqrt, ...
                                                            iterNum, iterMean, iterCov, iterCovSqrt, ...
                                                            measModel, dimNoise, dimMeas, numMeas, addNoiseMean, addNoiseCov, noiseMean, noiseCovSqrt)
            % Generate state and noise samples
            [stateSamples, ...
             noiseSamples, ...
             weights, ...
             numSamples] = Utils.getStateNoiseSamples(obj.samplingUpdate, ...
                                                      iterMean, iterCovSqrt, ...
                                                      noiseMean, noiseCovSqrt);
            
            measSamples = obj.getStackedMeasSamples(measModel, stateSamples, noiseSamples, ...
                                                    numSamples, dimMeas, numMeas, dimNoise);
            
            [measMean, measCov, ...
             stateMeasCrossCov] = Utils.getMeanCovAndCrossCov(iterMean, stateSamples, ...
                                                              measSamples, weights);
            
            % Compute measurement mean
            measMean = measMean + repmat(addNoiseMean, numMeas, 1);
            
            % Compute measurement covariance
            measCov = measCov + Utils.blockDiag(addNoiseCov, numMeas);
            
            if iterNum > 1
                [measMean, measCov, ...
                 stateMeasCrossCov] = KF.momentCorrection(priorMean, priorCov, priorCovSqrt, ...
                                                          iterMean, iterCov, iterCovSqrt, ...
                                                          measMean, measCov, stateMeasCrossCov);
            end
        end
        
        function measSamples = getStackedMeasSamples(obj, measModel, stateSamples, noiseSamples, ...
                                                     numSamples, dimMeas, numMeas, dimNoise)
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
        end
    end
    
    methods (Access = 'protected')
        function cpObj = copyElement(obj)
            cpObj = obj.copyElement@SampleBasedJointlyGaussianPrediction();
            
            if obj.samplingPrediction == obj.samplingUpdate
                % Still use the same samplings for both prediction and update
                cpObj.samplingUpdate = cpObj.samplingPrediction;
            else
                cpObj.samplingUpdate = obj.samplingUpdate.copy();
            end
        end
    end
    
    properties (SetAccess = 'private', GetAccess = 'protected')
        % Gaussian sampling technique used for the measurement update.
        samplingUpdate;
    end
end

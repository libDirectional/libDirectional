
classdef SampleBasedRecursiveUpdateFilter < RecursiveUpdateFilter & SampleBasedLinearGaussianFilter
    % Abstract base class for sample-based recursive update filters.
    %
    % SampleBasedRecursiveUpdateFilter Methods:
    %   SampleBasedRecursiveUpdateFilter - Class constructor.
    %   copy                             - Copy a Filter instance.
    %   copyWithName                     - Copy a Filter instance and give the copy a new name/description.
    %   getName                          - Get the filter name/description.
    %   setColor                         - Set the filter color/plotting properties.
    %   getColor                         - Get the filter color/plotting properties.
    %   setState                         - Set the system state.
    %   getState                         - Get the system state.
    %   setStateMeanAndCov               - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov               - Get mean and covariance matrix of the system state.
    %   getStateDim                      - Get the dimension of the system state.
    %   predict                          - Perform a state prediction.
    %   update                           - Perform a measurement update.
    %   step                             - Perform a combined state prediction and measurement update.
    %   setStateDecompDim                - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing      - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing      - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing          - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing          - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold           - Set the measurement gating threshold.
    %   getMeasGatingThreshold           - Get the measurement gating threshold.
    %   setNumRecursionSteps             - Set the number of recursion steps that are performed by a measurement update.
    %   getNumRecursionSteps             - Get the number of recursion steps that are performed by a measurement update.
    
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
    
    % Literature:
    %   Yulong Huang, Yonggang Zhang, Ning Li, and Lin Zhao,
    %   Design of Sigma-Point Kalman Filter with Recursive Updated Measurement,
    %   Circuits, Systems, and Signal Processing, pp. 1-16, Aug. 2015.
    
    methods
        function obj = SampleBasedRecursiveUpdateFilter(name)
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
            %   << obj (SampleBasedRecursiveUpdateFilter)
            %      A new SampleBasedRecursiveUpdateFilter instance.
            
            % Call superclass constructors
            obj = obj@RecursiveUpdateFilter(name);
            obj = obj@SampleBasedLinearGaussianFilter(name);
        end
    end
    
    methods (Sealed, Access = 'protected')
        function setupMeasModel(obj, measModel, dimMeas)
            [noiseMean, noiseCov, noiseCovSqrt] = measModel.noise.getMeanAndCov();
            
            obj.momentFuncHandle = @(stateMean, stateCovSqrt) ...
                                   obj.momentFuncMeasModel(measModel, dimMeas, ...
                                                           noiseMean, noiseCovSqrt, ...
                                                           stateMean, stateCovSqrt);
            
            obj.momentFuncCorrHandle = @(stateMean, stateCov, ...
                                         stateCovSqrt, stateNoiseCrossCov) ...
                                       obj.momentFuncCorrMeasModel(measModel, dimMeas, ...
                                                                   noiseMean, noiseCov, ...
                                                                   stateMean, stateCov, stateCovSqrt, ...
                                                                   stateNoiseCrossCov);
        end
        
        function setupAddNoiseMeasModel(obj, measModel, dimMeas)
            [addNoiseMean, addNoiseCov] = measModel.noise.getMeanAndCov();
            dimAddNoise = size(addNoiseMean, 1);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimAddNoise);
            
            obj.momentFuncHandle = @(stateMean, stateCovSqrt) ...
                                   obj.momentFuncAddNoiseMeasModel(measModel, dimMeas, ...
                                                                   addNoiseMean, addNoiseCov, ...
                                                                   stateMean, stateCovSqrt);
            
            obj.momentFuncCorrHandle = @(stateMean, stateCov, ...
                                         stateCovSqrt, stateNoiseCrossCov) ...
                                       obj.momentFuncCorrAddNoiseMeasModel(measModel, dimMeas, ...
                                                                           addNoiseMean, addNoiseCov, ...
                                                                           stateMean, stateCov, stateCovSqrt, ...
                                                                           stateNoiseCrossCov);
        end
        
        function setupMixedNoiseMeasModel(obj, measModel, dimMeas)
            [noiseMean, noiseCov, noiseCovSqrt] = measModel.noise.getMeanAndCov();
            [addNoiseMean, addNoiseCov]  = measModel.additiveNoise.getMeanAndCov();
            dimAddNoise = size(addNoiseMean, 1);
            
            obj.checkAdditiveMeasNoise(dimMeas, dimAddNoise);
            
            obj.momentFuncHandle = @(stateMean, stateCovSqrt) ...
                                   obj.momentFuncMixedNoiseMeasModel(measModel, dimMeas, ...
                                                                     noiseMean, noiseCovSqrt, ...
                                                                     addNoiseMean, addNoiseCov, ...
                                                                     stateMean, stateCovSqrt);
         	
            obj.momentFuncCorrHandle = @(stateMean, stateCov, ...
                                         stateCovSqrt, stateNoiseCrossCov) ...
                                       obj.momentFuncCorrMixedNoiseMeasModel(measModel, dimMeas, ...
                                                                             noiseMean, noiseCov, ...
                                                                             addNoiseMean, addNoiseCov, ...
                                                                             stateMean, stateCov, stateCovSqrt, ...
                                                                             stateNoiseCrossCov);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoiseCrossCov] = getMeasMoments(obj, priorStateMean, ~, priorStateCovSqrt)
            [measMean, measCov, ...
             stateMeasCrossCov, ...
             measNoiseCrossCov] = obj.momentFuncHandle(priorStateMean, priorStateCovSqrt);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoiseCrossCov] = getMeasMomentsCorr(obj, updatedStateMean, updatedStateCov, updatedStateCovSqrt, ...
                                                          updatedStateNoiseCrossCov)
            [measMean, measCov, ...
             stateMeasCrossCov, ...
             measNoiseCrossCov] = obj.momentFuncCorrHandle(updatedStateMean, updatedStateCov, updatedStateCovSqrt, ...
                                                           updatedStateNoiseCrossCov);
        end
    end
    
    methods (Access = 'private')
        % Measurement update moment computation with uncorrelated system state and measurement noise
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoiseCrossCov] = momentFuncMeasModel(obj, measModel, dimMeas, ...
                                                           noiseMean, noiseCovSqrt, ...
                                                           stateMean, stateCovSqrt)
            % No correlations between system state and measurement noise:
            % E[(x - E[x])*(v - E[v])'] = 0
            
            [hSamples, weights, ...
             zeroMeanStateSamples, ...
             zeroMeanNoiseSamples] = obj.evaluateMeasModelUncorr(measModel, dimMeas, ...
                                                                 stateMean, stateCovSqrt, ...
                                                                 noiseMean, noiseCovSqrt);
            
            % Compute moments:
            %  - Measurement mean:
            %    E[y] = E[h(x,v)]
            %  - Measurement covariance matrix:
            %    E[(y - E[y])*(y - E[y])'] = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %  - State--measurement cross-covariance matrix:
            %    E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %  - Measurement--noise cross-covariance matrix:
            %    E[(y - E[y])*(v - E[v])'] = E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            [measMean, measCov, ...
             stateMeasCrossCov, ...
             measNoiseCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                            zeroMeanStateSamples, ...
                                                            zeroMeanNoiseSamples);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoiseCrossCov] = momentFuncAddNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                   addNoiseMean, addNoiseCov, ...
                                                                   stateMean, stateCovSqrt)
            % No correlations between system state and measurement noise:
            % E[(x - E[x])*(v - E[v])'] = 0
            
            [hSamples, weights, ...
             zeroMeanStateSamples] = obj.evaluateAddNoiseMeasModel(measModel, dimMeas, ...
                                                                   stateMean, stateCovSqrt);
          	
            % Compute moments:
            %  - E[h(x)]
            %  - E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %  - E[(x - E[x])*(h(x) - E[h(x)])']
            [hMean, hCov, ...
             stateHCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                         zeroMeanStateSamples);
            
            % Measurement mean:
            % E[y] = E[h(x)] + E[v]
            measMean = hMean + addNoiseMean;
            
            % Measurement covariance matrix:
            % E[(y - E[y])*(y - E[y])'] = E[(h(x) + v - E[h(x)] - E[v])*(h(x) + v - E[h(x)] - E[v])']
            %                           = E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %                           + E[(v - E[v])*(v - E[v])']
            %                           + E[(h(x) - E[h(x)])*(v - E[v])']
            %                           + E[(v - E[v])*(h(x) - E[h(x)])']
            %                           = E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %                           + E[(v - E[v])*(v - E[v])']
            measCov = hCov + addNoiseCov;
            
            % State--measurement cross-covariance matrix
            % E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x) + v - E[h(x)] - E[v])']
            %                           = E[(x - E[x])*(h(x) - E[h(x)])']
            %                           + E[(x - E[x])*(v - E[v])']
            %                           = E[(x - E[x])*(h(x) - E[h(x)])']
            stateMeasCrossCov = stateHCrossCov;
            
            % Measurement--noise cross-covariance matrix
            % E[(y - E[y])*(v - E[v])'] = E[(h(x) + v - E[h(x)] - E[v])*(v - E[v])']
            %                           = E[(h(x) - E[h(x)])*(v - E[v])']
            %                           + E[(v - E[v])*(v - E[v])']
            %                           = E[(v - E[v])*(v - E[v])']
            measNoiseCrossCov = addNoiseCov;
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoisesCrossCov] = momentFuncMixedNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                      noiseMean, noiseCovSqrt, ...
                                                                      addNoiseMean, addNoiseCov, ...
                                                                      stateMean, stateCovSqrt)
            % No correlations between system state and measurement noise:
            % E[(x - E[x])*(v - E[v])'] = 0
            % E[(x - E[x])*(r - E[r])'] = 0
            
            [hSamples, weights, ...
             zeroMeanStateSamples, ...
             zeroMeanNoiseSamples] = obj.evaluateMeasModelUncorr(measModel, dimMeas, ...
                                                                 stateMean, stateCovSqrt, ...
                                                                 noiseMean, noiseCovSqrt);
            
            % Compute moments:
            %  - E[h(x,v)]
            %  - E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %  - E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %  - E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            [hMean, hCov, ...
             stateHCrossCov, ...
             hNoiseCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                         zeroMeanStateSamples, ...
                                                         zeroMeanNoiseSamples);
            
            % Measurement mean:
            % E[y] = E[h(x,v)] + E[r]
            measMean = hMean + addNoiseMean;
            
            % Measurement covariance matrix:
            % E[(y - E[y])*(y - E[y])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(h(x,v) + r - E[h(x,v)] - E[r])']
            %                           = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %                           + E[(r - E[r])*(r - E[r])']
            %                           + E[(h(x,v) - E[h(x,v)])*(r - E[r])']
            %                           + E[(r - E[r])*(h(x,v) - E[h(x,v)])']
            %                           = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %                           + E[(r - E[r])*(r - E[r])']
            measCov = hCov + addNoiseCov;
            
            % State--measurement cross-covariance matrix
            % E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x,v) + r - E[h(x,v)] - E[r])']
            %                           = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %                           + E[(x - E[x])*(r - E[r])']
            %                           = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            stateMeasCrossCov = stateHCrossCov;
            
            % Measurement--noise cross-covariance matrix
            % E[(y - E[y])*(v - E[v])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(v - E[v])']
            %                           = E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            %                           + E[(r - E[r])*(v - E[v])']
            %                           = E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            %
            % E[(y - E[y])*(r - E[r])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(r - E[r])']
            %                           = E[(h(x,v) - E[h(x,v)])*(r - E[r])']
            %                           + E[(r - E[r])*(r - E[r])']
            %                           = E[(r - E[r])*(r - E[r])']
            measNoisesCrossCov = [hNoiseCrossCov addNoiseCov];
        end
        
        % Measurement update moment computation with correlated system state and measurement noise
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoiseCrossCov] = momentFuncCorrMeasModel(obj, measModel, dimMeas, ...
                                                               noiseMean, noiseCov, ...
                                                               stateMean, stateCov, ~, ...
                                                               stateNoiseCrossCov)
            % Correlations between system state and measurement noise:
            % E[(x - E[x])*(v - E[v])'] != 0
            
            [hSamples, weights, ...
             zeroMeanStateSamples, ...
             zeroMeanNoiseSamples] = obj.evaluateMeasModelCorr(measModel, dimMeas, ...
                                                               stateMean, stateCov, ...
                                                               noiseMean, noiseCov, ...
                                                               stateNoiseCrossCov);
            
            % Compute moments:
            %  - Measurement mean:
            %    E[y] = E[h(x,v)]
            %  - Measurement covariance matrix:
            %    E[(y - E[y])*(y - E[y])'] = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %  - State--measurement cross-covariance matrix:
            %    E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %  - Measurement--noise cross-covariance matrix:
            %    E[(y - E[y])*(v - E[v])'] = E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            [measMean, measCov, ...
             stateMeasCrossCov, ...
             measNoiseCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                            zeroMeanStateSamples, ...
                                                            zeroMeanNoiseSamples);
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoiseCrossCov] = momentFuncCorrAddNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                       addNoiseMean, addNoiseCov, ...
                                                                       stateMean, stateCov, stateCovSqrt, ...
                                                                       stateAddNoiseCrossCov)
            % Correlations between system state and measurement noise:
            % E[(x - E[x])*(v - E[v])'] != 0
            
            [hSamples, weights, ...
             zeroMeanStateSamples] = obj.evaluateAddNoiseMeasModel(measModel, dimMeas, ...
                                                                   stateMean, stateCovSqrt);
          	
            % Compute moments:
            %  - E[h(x)]
            %  - E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %  - E[(x - E[x])*(h(x) - E[h(x)])']
            [hMean, hCov, ...
             stateHCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                         zeroMeanStateSamples);
            
            % Measurement mean:
            % E[y] = E[h(x)] + E[v]
            measMean = hMean + addNoiseMean;
            
            % E[(x - E[x])*(v - E[v])'] != 0 => E[(h(x) - E[h(x)])*(v - E[v])'] != 0
            H = stateHCrossCov' / stateCov;
            hAddNoiseCrossCov = H * stateAddNoiseCrossCov;
            
            % Measurement covariance matrix:
            % E[(y - E[y])*(y - E[y])'] = E[(h(x) + v - E[h(x)] - E[v])*(h(x) + v - E[h(x)] - E[v])']
            %                           = E[(h(x) - E[h(x)])*(h(x) - E[h(x)])']
            %                           + E[(v - E[v])*(v - E[v])']
            %                           + E[(h(x) - E[h(x)])*(v - E[v])']
            %                           + E[(v - E[v])*(h(x) - E[h(x)])']
            measCov = hCov + addNoiseCov + (hAddNoiseCrossCov + hAddNoiseCrossCov');
            
            % State--measurement cross-covariance matrix
            % E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x) + v - E[h(x)] - E[v])']
            %                           = E[(x - E[x])*(h(x) - E[h(x)])']
            %                           + E[(x - E[x])*(v - E[v])']
            stateMeasCrossCov = stateHCrossCov + stateAddNoiseCrossCov;
            
            % Measurement--noise cross-covariance matrix
            % E[(y - E[y])*(v - E[v])'] = E[(h(x) + v - E[h(x)] - E[v])*(v - E[v])']
            %                           = E[(h(x) - E[h(x)])*(v - E[v])']
            %                           + E[(v - E[v])*(v - E[v])']
            measNoiseCrossCov = hAddNoiseCrossCov + addNoiseCov;
        end
        
        function [measMean, measCov, ...
                  stateMeasCrossCov, ...
                  measNoisesCrossCov] = momentFuncCorrMixedNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                          noiseMean, noiseCov, ...
                                                                          addNoiseMean, addNoiseCov, ...
                                                                          stateMean, stateCov, ~, ...
                                                                          stateNoisesCrossCov)
            % Correlations between system state and measurement noise:
            % E[(x - E[x])*(v - E[v])'] != 0
            % E[(x - E[x])*(r - E[r])'] != 0
            
            dimState              = size(stateMean, 1);
            dimNoise              = size(noiseMean, 1);
            stateNoiseCrossCov    = stateNoisesCrossCov(:, 1:dimNoise);
            stateAddNoiseCrossCov = stateNoisesCrossCov(:, dimNoise+1:end);
            
            [hSamples, weights, ...
             zeroMeanStateSamples, ...
             zeroMeanNoiseSamples, ...
             stateNoiseCov] = obj.evaluateMeasModelCorr(measModel, dimMeas, ...
                                                        stateMean, stateCov, ...
                                                        noiseMean, noiseCov, ...
                                                        stateNoiseCrossCov);
            
            % Compute moments:
            %  - E[h(x,v)]
            %  - E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %  - E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %  - E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            [hMean, hCov, ...
             stateHCrossCov, ...
             hNoiseCrossCov] = obj.getMeasModelMomements(hSamples, weights, ...
                                                         zeroMeanStateSamples, ...
                                                         zeroMeanNoiseSamples);
            
            % Measurement mean:
            % E[y] = E[h(x,v)] + E[r]
            measMean = hMean + addNoiseMean;
            
            % E[(x - E[x])*(r - E[r])'] != 0 => E[(h(x,v) - E[h(x,v)])*(r - E[r])'] != 0
            hStateNoiseCrossCov = [stateHCrossCov' hNoiseCrossCov];
            H = hStateNoiseCrossCov / stateNoiseCov;
            % Note that E[(v - E[v])*(r - E[r])'] = 0
            hAddNoiseCrossCov = H(:, 1:dimState) * stateAddNoiseCrossCov;
            
            % Measurement covariance matrix:
            % E[(y - E[y])*(y - E[y])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(h(x,v) + r - E[h(x,v)] - E[r])']
            %                           = E[(h(x,v) - E[h(x,v)])*(h(x,v) - E[h(x,v)])']
            %                           + E[(r - E[r])*(r - E[r])']
            %                           + E[(h(x,v) - E[h(x,v)])*(r - E[r])']
            %                           + E[(r - E[r])*(h(x,v) - E[h(x,v)])']
            measCov = hCov + addNoiseCov + (hAddNoiseCrossCov + hAddNoiseCrossCov');
            
            % State--measurement cross-covariance matrix
            % E[(x - E[x])*(y - E[y])'] = E[(x - E[x])*(h(x,v) + r - E[h(x,v)] - E[r])']
            %                           = E[(x - E[x])*(h(x,v) - E[h(x,v)])']
            %                           + E[(x - E[x])*(r - E[r])']
            stateMeasCrossCov = stateHCrossCov + stateAddNoiseCrossCov;
            
            % Measurement--noise cross-covariance matrix
            % E[(y - E[y])*(v - E[v])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(v - E[v])']
            %                           = E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            %                           + E[(r - E[r])*(v - E[v])']
            %                           = E[(h(x,v) - E[h(x,v)])*(v - E[v])']
            %
            % E[(y - E[y])*(r - E[r])'] = E[(h(x,v) + r - E[h(x,v)] - E[r])*(r - E[r])']
            %                           = E[(h(x,v) - E[h(x,v)])*(r - E[r])']
            %                           + E[(r - E[r])*(r - E[r])']
            measNoisesCrossCov = [hNoiseCrossCov hAddNoiseCrossCov + addNoiseCov];
        end
    end
    
    properties (Access = 'private')
        % Function handle to the currently used moment computation method
        % for uncorrelated system state and measurement noise.
        momentFuncHandle;
        
        % Function handle to the currently used moment computation method
        % for correlated system state and measurement noise.
        momentFuncCorrHandle;
    end
end

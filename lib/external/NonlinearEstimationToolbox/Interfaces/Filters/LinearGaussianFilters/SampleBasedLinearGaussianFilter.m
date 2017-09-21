
classdef SampleBasedLinearGaussianFilter < LinearGaussianFilter
    % Abstract base class for sample-based linear Gaussian filters.
    %
    % SampleBasedLinearGaussianFilter Methods:
    %   SampleBasedLinearGaussianFilter - Class constructor.
    %   copy                            - Copy a Filter instance.
    %   copyWithName                    - Copy a Filter instance and give the copy a new name/description.
    %   getName                         - Get the filter name/description.
    %   setColor                        - Set the filter color/plotting properties.
    %   getColor                        - Get the filter color/plotting properties.
    %   setState                        - Set the system state.
    %   getState                        - Get the system state.
    %   setStateMeanAndCov              - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov              - Get mean and covariance matrix of the system state.
    %   getStateDim                     - Get the dimension of the system state.
    %   predict                         - Perform a state prediction.
    %   update                          - Perform a measurement update.
    %   step                            - Perform a combined state prediction and measurement update.
    %   setStateDecompDim               - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim               - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing     - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing     - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing         - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing         - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold          - Set the measurement gating threshold.
    %   getMeasGatingThreshold          - Get the measurement gating threshold.
    
    % Literature:
    %   Michael Roth, Gustaf Hendeby, and Fredrik Gustafsson,
    %   Nonlinear Kalman Filters Explained: A Tutorial on Moment Computations and Sigma Point Methods,
    %   Journal of Advances in Information Fusion, vol. 11, no. 1, pp. 47-70, Jun. 2016.
    %
    %   Jannik Steinbring and Uwe D. Hanebeck,
    %   LRKF Revisited: The Smart Sampling Kalman Filter (S²KF),
    %   Journal of Advances in Information Fusion, vol. 9, no. 2, pp. 106-123, Dec. 2014.
    
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
        function obj = SampleBasedLinearGaussianFilter(name)
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
            %   << obj (SampleBasedLinearGaussianFilter)
            %      A new SampleBasedLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@LinearGaussianFilter(name);
        end
    end
    
    methods (Sealed, Access = 'protected')
        % State prediction moment computation
        function [predictedStateMean, ...
                  predictedStateCov] = predictSysModel(obj, sysModel)
            [noiseMean, ~, noiseCovSqrt] = sysModel.noise.getMeanAndCov();
            dimNoise    = size(noiseMean, 1);
            dimAugState = obj.dimState + dimNoise;
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = obj.getStdNormalSamplesPrediction(dimAugState);
            
            % Generate state samples
            zeroMeanStateSamples = obj.stateCovSqrt * stdNormalSamples(1:obj.dimState, :);
            stateSamples         = bsxfun(@plus, zeroMeanStateSamples, obj.stateMean);
            
            % Generate noise samples
            zeroMeanNoiseSamples = noiseCovSqrt * stdNormalSamples(obj.dimState+1:end, :);
            noiseSamples         = bsxfun(@plus, zeroMeanNoiseSamples, noiseMean);
            
            % Propagate samples through system equation a(x, w)
            aSamples = sysModel.systemEquation(stateSamples, noiseSamples);
            
            % Check predicted state samples
            obj.checkPredictedStateSamples(aSamples, numSamples);
            
            % Compute moments:
            %  - Predicted state mean:
            %    E[xp] = E[a(x,w)]
            %  - Predicted state covariance matrix:
            %    E[(xp - E[xp])*(xp - E[xp])'] = E[(a(x,w) - E[a(x,w)])*(a(x,w) - E[a(x,w)])']
            [predictedStateMean, ...
             predictedStateCov] = Utils.getMeanAndCov(aSamples, weights);
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictAddNoiseSysModel(obj, sysModel)
            [noiseMean, noiseCov] = sysModel.noise.getMeanAndCov();
            dimNoise = size(noiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimNoise);
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = obj.getStdNormalSamplesPrediction(obj.dimState);
            
            % Generate state samples
            zeroMeanStateSamples = obj.stateCovSqrt * stdNormalSamples;
            stateSamples         = bsxfun(@plus, zeroMeanStateSamples, obj.stateMean);
            
            % Propagate samples through system equation a(x)
            aSamples = sysModel.systemEquation(stateSamples);
            
            % Check predicted state samples
            obj.checkPredictedStateSamples(aSamples, numSamples);
            
            % Compute moments:
            %  - E[a(x)]
            %  - E[(a(x) - E[a(x)])*(a(x) - E[a(x)])']
            [mean, cov] = Utils.getMeanAndCov(aSamples, weights);
            
            % Predicted state mean:
            % E[xp] = E[a(x)] + E[w]
            predictedStateMean = mean + noiseMean;
            
            % Predicted state covariance matrix:
            % E[(xp - E[xp])*(xp - E[xp])'] = E[(a(x) + w - E[a(x)] - E[w])*(a(x) + w - E[a(x)] - E[w])']
            %                               = E[(a(x) - E[a(x)])*(a(x) - E[a(x)])']
            %                               + E[(w - E[w])(w - E[w])']
            predictedStateCov = cov + noiseCov;
        end
        
        % Helpers for the sample-based linear measurement updates
        function [hSamples, weights, ...
                  zeroMeanStateSamples] = evaluateAddNoiseMeasModel(obj, measModel, dimMeas, ...
                                                                    stateMean, stateCovSqrt)
            dimState = size(stateMean, 1);
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = obj.getStdNormalSamplesUpdate(dimState);
            
            % Generate state samples
            zeroMeanStateSamples = stateCovSqrt * stdNormalSamples;
            stateSamples         = bsxfun(@plus, zeroMeanStateSamples, stateMean);
            
            % Propagate samples through measurement equation h(x)
            hSamples = measModel.measurementEquation(stateSamples);
            
            % Check computed measurements
            obj.checkComputedMeasurements(hSamples, dimMeas, numSamples);
        end
        
        function [hSamples, weights, ...
                  zeroMeanStateSamples, ...
                  zeroMeanNoiseSamples] = evaluateMeasModelUncorr(obj, measModel, dimMeas, ...
                                                                  stateMean, stateCovSqrt, ...
                                                                  noiseMean, noiseCovSqrt)
            dimState    = size(stateMean, 1);
            dimNoise    = size(noiseMean, 1);
            dimAugState = dimState + dimNoise;
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = obj.getStdNormalSamplesUpdate(dimAugState);
            
            % Generate state samples
            zeroMeanStateSamples = stateCovSqrt * stdNormalSamples(1:dimState, :);
            stateSamples         = bsxfun(@plus, zeroMeanStateSamples, stateMean);
            
            % Generate noise samples
            zeroMeanNoiseSamples = noiseCovSqrt * stdNormalSamples(dimState+1:end, :);
            noiseSamples         = bsxfun(@plus, zeroMeanNoiseSamples, noiseMean);
            
            % Propagate samples through measurement equation h(x,v)
            hSamples = measModel.measurementEquation(stateSamples, noiseSamples);
            
            % Check computed measurements
            obj.checkComputedMeasurements(hSamples, dimMeas, numSamples);
        end
        
        function [hSamples, weights, ...
                  zeroMeanStateSamples, ...
                  zeroMeanNoiseSamples, ...
                  stateNoiseCov] = evaluateMeasModelCorr(obj, measModel, dimMeas, ...
                                                         stateMean, stateCov, ...
                                                         noiseMean, noiseCov, ...
                                                         stateNoiseCrossCov)
            dimState    = size(stateMean, 1);
            dimNoise    = size(noiseMean, 1);
            dimAugState = dimState + dimNoise;
            
            stateNoiseCov = [stateCov            stateNoiseCrossCov
                             stateNoiseCrossCov' noiseCov          ];
            
            stateNoiseCovSqrt = obj.checkCovUpdate(stateNoiseCov, 'System state and measurement noise joint');
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = obj.getStdNormalSamplesUpdate(dimAugState);
            
            % Generate state samples
            zeroMeanStateSamples = stateNoiseCovSqrt(1:dimState, :) * stdNormalSamples;
            stateSamples         = bsxfun(@plus, zeroMeanStateSamples, stateMean);
            
            % Generate noise samples
            zeroMeanNoiseSamples = stateNoiseCovSqrt(dimState+1:end, :) * stdNormalSamples;
            noiseSamples         = bsxfun(@plus, zeroMeanNoiseSamples, noiseMean);
            
            % Propagate samples through measurement equation h(x,v)
            hSamples = measModel.measurementEquation(stateSamples, noiseSamples);
            
            % Check computed measurements
            obj.checkComputedMeasurements(hSamples, dimMeas, numSamples);
        end
        
        function [hMean, hCov, ...
                  stateHCrossCov, ...
                  hNoiseCrossCov] = getMeasModelMomements(~, hSamples, weights, ...
                                                          zeroMeanStateSamples, ...
                                                          zeroMeanNoiseSamples)
            if numel(weights) == 1
                % Equally weighted samples
                
                hMean = sum(hSamples, 2) * weights;
                
                zeroMeanHSamples = bsxfun(@minus, hSamples, hMean);
                
                hCov = (zeroMeanHSamples * zeroMeanHSamples') * weights;
                
                stateHCrossCov = (zeroMeanStateSamples * zeroMeanHSamples') * weights;
                
                if nargin == 5
                    hNoiseCrossCov = (zeroMeanHSamples * zeroMeanNoiseSamples') * weights;
                end
            else
                % Samples are not equally weighted
                
                hMean = hSamples * weights';
                
                zeroMeanHSamples = bsxfun(@minus, hSamples, hMean);
                
                % Weights can be negative => we have to treat them separately
                
                % Positive weights
                idx = weights >= 0;
                
                sqrtWeights                  = sqrt(weights(idx));
                weightedZeroMeanHSamples     = bsxfun(@times, zeroMeanHSamples(:, idx), sqrtWeights);
                weightedZeroMeanStateSamples = bsxfun(@times, zeroMeanStateSamples(:, idx), sqrtWeights);
                
                hCov = weightedZeroMeanHSamples * weightedZeroMeanHSamples';
                
                stateHCrossCov = weightedZeroMeanStateSamples * weightedZeroMeanHSamples';
                
                if nargin == 5
                    weightedZeroMeanNoiseSamples = bsxfun(@times, zeroMeanNoiseSamples(:, idx), sqrtWeights);
                    
                    hNoiseCrossCov = weightedZeroMeanHSamples * weightedZeroMeanNoiseSamples';
                end
                
                % Negative weights
                if ~all(idx)
                    idx = ~idx;
                    
                    sqrtWeights                  = sqrt(abs(weights(idx)));
                    weightedZeroMeanHSamples     = bsxfun(@times, zeroMeanHSamples(:, idx), sqrtWeights);
                    weightedZeroMeanStateSamples = bsxfun(@times, zeroMeanStateSamples(:, idx), sqrtWeights);
                    
                    hCov = hCov - weightedZeroMeanHSamples * weightedZeroMeanHSamples';
                    
                    stateHCrossCov = stateHCrossCov - weightedZeroMeanStateSamples * weightedZeroMeanHSamples';
                    
                    if nargin == 5
                        weightedZeroMeanNoiseSamples = bsxfun(@times, zeroMeanNoiseSamples(:, idx), sqrtWeights);
                        
                        hNoiseCrossCov = hNoiseCrossCov - weightedZeroMeanHSamples * weightedZeroMeanNoiseSamples';
                    end
                end
            end
        end
    end
    
    methods (Abstract, Access = 'protected')
        [stdNormalSamples, ...
         weights, numSamples] = getStdNormalSamplesPrediction(obj, dim);
        
        [stdNormalSamples, ...
         weights, numSamples] = getStdNormalSamplesUpdate(obj, dim);
    end
end

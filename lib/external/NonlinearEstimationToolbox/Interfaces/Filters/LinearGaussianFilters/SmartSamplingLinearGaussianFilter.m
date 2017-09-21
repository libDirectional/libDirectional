
classdef SmartSamplingLinearGaussianFilter < SampleBasedLinearGaussianFilter
    % Abstract base class for linear Gaussian filters that are based on the
    % point-symmetric LCD-based Gaussian sampling technique.
    %
    % SmartSamplingLinearGaussianFilter Methods:
    %   SmartSamplingLinearGaussianFilter - Class constructor.
    %   copy                              - Copy a Filter instance.
    %   copyWithName                      - Copy a Filter instance and give the copy a new name/description.
    %   getName                           - Get the filter name/description.
    %   setColor                          - Set the filter color/plotting properties.
    %   getColor                          - Get the filter color/plotting properties.
    %   setState                          - Set the system state.
    %   getState                          - Get the system state.
    %   setStateMeanAndCov                - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov                - Get mean and covariance matrix of the system state.
    %   getStateDim                       - Get the dimension of the system state.
    %   predict                           - Perform a state prediction.
    %   update                            - Perform a measurement update.
    %   step                              - Perform a combined state prediction and measurement update.
    %   setStateDecompDim                 - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                 - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing       - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing       - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing           - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing           - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold            - Set the measurement gating threshold.
    %   getMeasGatingThreshold            - Get the measurement gating threshold.
    %   setNumSamples                     - Set absolute numbers of samples used for state prediction and measurement update.
    %   setNumSamplesByFactors            - Set linear factors to determine the number of samples used for state prediction and measurement update.
    %   getNumSamplesConfigPrediction     - Get the number of samples configuration used for the state prediction.
    %   getNumSamplesConfigUpdate         - Get the number of samples configuration used for the measurement update.
    
    % Literature:
    %   Jannik Steinbring, Martin Pander, and Uwe D. Hanebeck,
    %   The Smart Sampling Kalman Filter with Symmetric Samples
    %   Journal of Advances in Information Fusion, vol. 11, no. 1, pp. 71-90, Jun. 2016.
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
    
    methods (Sealed)
        function obj = SmartSamplingLinearGaussianFilter(name)
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
            %   << obj (SmartSamplingLinearGaussianFilter)
            %      A new SmartSamplingLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@SampleBasedLinearGaussianFilter(name);
            
            obj.samplingPrediction = GaussianSamplingLCD();
            obj.samplingUpdate     = GaussianSamplingLCD();
            
            % The point-symmetric LCD-based sampling is used.
            obj.samplingPrediction.setSymmetricMode(true);
            obj.samplingUpdate.setSymmetricMode(true);
            
            % By default, determine the number of samples for state prediction
            % and measurement update by using a linear factor of 10.
            obj.samplingPrediction.setNumSamplesByFactor(10);
            obj.samplingUpdate.setNumSamplesByFactor(10);
        end
        
        function setNumSamples(obj, numSamplesPrediction, numSamplesUpdate)
            % Set absolute numbers of samples used for state prediction and measurement update.
            %
            % This also overwrites a possible previous setting, where the number of
            % samples are determined by a linear factor (see setNumSamplesByFactor()).
            %
            % By default, a linear factor 10 is used for both state prediction
            % and measurement update.
            %
            % Parameters:
            %   >> numSamplesPrediction (Positive scalar or empty matrix)
            %      The new absolute number of samples used for the state prediction.
            %      Pass an empty matrix to keep the current configuration for the
            %      state prediction.
            %
            %   >> numSamplesUpdate (Positive scalar or empty matrix)
            %      The new absolute number of samples used for the measurement update.
            %      Pass an empty matrix to keep the current configuration for the
            %      measurement update. If nothing is passed, the absolute number of
            %      samples specified for the state prediction is also used for the
            %      measurement update.
            
            if ~isempty(numSamplesPrediction)
                obj.samplingPrediction.setNumSamples(numSamplesPrediction);
            end
            
            if nargin == 3
                if ~isempty(numSamplesUpdate)
                    obj.samplingUpdate.setNumSamples(numSamplesUpdate);
                end
            elseif ~isempty(numSamplesPrediction)
                obj.samplingUpdate.setNumSamples(numSamplesPrediction);
            end
        end
        
        function setNumSamplesByFactors(obj, factorPrediction, factorUpdate)
            % Set linear factors to determine the number of samples used for state prediction and measurement update.
            %
            % The actual number of samples will be computed according to
            %
            %    Number of samples = factor * dimension + 1 - mod(factor * dimension, 2)
            %
            % i.e., always an odd number of samples is used.
            %
            % This also overwrites a possible previous setting, where the number of
            % samples are determined in an absolute way (see setNumSamples()).
            %
            % By default, a linear factor of 10 is used for both state prediction
            % and measurement update.
            %
            % Parameters:
            %   >> factorPrediction (Positive scalar or empty matrix)
            %      The new linear factor to determine the number of samples used for
            %      the state prediction. Pass an empty matrix to keep the current
            %      configuration for the state prediction.
            %
            %   >> factorUpdate (Positive scalar or empty matrix)
            %      The new linear factor to determine the number of samples used for
            %      the measurement update. Pass an empty matrix to keep the current
            %      configuration for the measurement update. If nothing is passed, the
            %      linear factor specified for the state prediction is also used for
            %      the measurement update.
            
            if ~isempty(factorPrediction)
                obj.samplingPrediction.setNumSamplesByFactor(factorPrediction);
            end
            
            if nargin == 3
                if ~isempty(factorUpdate)
                    obj.samplingUpdate.setNumSamplesByFactor(factorUpdate);
                end
            elseif ~isempty(factorPrediction)
                obj.samplingUpdate.setNumSamplesByFactor(factorPrediction);
            end
        end
        
        function [numSamplesAbs, ...
                  numSamplesFactor] = getNumSamplesConfigPrediction(obj)
            % Get the number of samples configuration used for the state prediction.
            %
            % Returns:
            %   << numSamplesAbs (Positive scalar or empty matrix)
            %      Equals the absolute number of samples if set.
            %      Otherwise, an empty matrix.
            %
            %   << numSamplesFactor (Positive scalar or empty matrix)
            %      Equals the sample factor if set.
            %      Otherwise, an empty matrix.
            
            [numSamplesAbs, ...
             numSamplesFactor] = obj.samplingPrediction.getNumSamplesConfig();
        end
        
        function [numSamplesAbs, ...
                  numSamplesFactor] = getNumSamplesConfigUpdate(obj)
            % Get the number of samples configuration used for the measurement update.
            %
            % Returns:
            %   << numSamplesAbs (Positive scalar or empty matrix)
            %      Equals the absolute number of samples if set.
            %      Otherwise, an empty matrix.
            %
            %   << numSamplesFactor (Positive scalar or empty matrix)
            %      Equals the sample factor if set.
            %      Otherwise, an empty matrix.
            
            [numSamplesAbs, ...
             numSamplesFactor] = obj.samplingUpdate.getNumSamplesConfig();
        end
    end
    
    methods (Sealed, Access = 'protected')
        function [stdNormalSamples, ...
                  weights, numSamples] = getStdNormalSamplesPrediction(obj, dim)
            [stdNormalSamples, ...
             weights, numSamples] = obj.samplingPrediction.getStdNormalSamples(dim);
        end
        
        function [stdNormalSamples, ...
                  weights, numSamples] = getStdNormalSamplesUpdate(obj, dim)
            [stdNormalSamples, ...
             weights, numSamples] = obj.samplingUpdate.getStdNormalSamples(dim);
        end
    end
    
    methods (Access = 'protected')
        % Copy Gaussian sampling techniques correctly
        function cpObj = copyElement(obj)
            cpObj = obj.copyElement@SampleBasedLinearGaussianFilter();
            
            cpObj.samplingPrediction = obj.samplingPrediction.copy();
            cpObj.samplingUpdate     = obj.samplingUpdate.copy();
        end
    end
    
    properties (Access = 'private')
        % Gaussian sampling technique used for the state prediction.
        samplingPrediction;
        
        % Gaussian sampling technique used for the measurement update.
        samplingUpdate;
    end
end

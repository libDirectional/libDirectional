
classdef RandomizedUnscentedLinearGaussianFilter < SampleBasedLinearGaussianFilter
    % Abstract base class for linear Gaussian filters that are based on the
    % randomized unscented Gaussian sampling technique.
    %
    % RandomizedUnscentedLinearGaussianFilter Methods:
    %   RandomizedUnscentedLinearGaussianFilter - Class constructor.
    %   copy                                    - Copy a Filter instance.
    %   copyWithName                            - Copy a Filter instance and give the copy a new name/description.
    %   getName                                 - Get the filter name/description.
    %   setColor                                - Set the filter color/plotting properties.
    %   getColor                                - Get the filter color/plotting properties.
    %   setState                                - Set the system state.
    %   getState                                - Get the system state.
    %   setStateMeanAndCov                      - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov                      - Get mean and covariance matrix of the system state.
    %   getStateDim                             - Get the dimension of the system state.
    %   predict                                 - Perform a state prediction.
    %   update                                  - Perform a measurement update.
    %   step                                    - Perform a combined state prediction and measurement update.
    %   setStateDecompDim                       - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim                       - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing             - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing             - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing                 - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing                 - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold                  - Set the measurement gating threshold.
    %   getMeasGatingThreshold                  - Get the measurement gating threshold.
    %   setNumSamplesFactors                    - Set the linear factors to determine the number of samples used for state prediction and measurement update.
    %   getNumSamplesFactors                    - Get the linear factors to determine the number of samples used for state prediction and measurement upate.
    
    % Literature:
    %   Jindrich Dunik, Ondrej Straka, and Miroslav Simandl,
    %   The Development of a Randomised Unscented Kalman Filter,
    %   Proceedings of the 18th IFAC World Congress, Milano, Italy, pp. 8-13, Aug. 2011.
    
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
        function obj = RandomizedUnscentedLinearGaussianFilter(name)
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
            %   << obj (RandomizedUnscentedLinearGaussianFilter)
            %      A new RandomizedUnscentedLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@SampleBasedLinearGaussianFilter(name);
            
            obj.samplingPrediction = GaussianSamplingRUKF();
            obj.samplingUpdate     = GaussianSamplingRUKF();
            
            % By default, 5 iterations are used for both
            % state prediction and measurement update.
            obj.samplingPrediction.setNumIterations(5);
            obj.samplingUpdate.setNumIterations(5);
        end
        
        function setNumSamplesFactors(obj, factorPrediction, factorUpdate)
            % Set the linear factors to determine the number of samples used for state prediction and measurement update.
            %
            % The actual number of samples will be computed according to
            %
            %    Number of samples = factor * 2 * dimension + 1
            %
            % By default, a linear factor of 5 is used for both state
            % prediction and measurement upate.
            %
            % Parameters:
            %   >> factorPrediction (Positive scalar)
            %      The new linear factor to determine the number of samples
            %      used for the state prediction.
            %
            %   >> factorUpdate (Positive scalar)
            %      The new linear factor to determine the number of samples
            %      used for the measurement update.
            %      If nothing is passed, the linear factor specified for the
            %      state prediction is also used for the measurement update.
            
            obj.samplingPrediction.setNumIterations(factorPrediction);
            
            if nargin == 3
                obj.samplingUpdate.setNumIterations(factorUpdate);
            else
                obj.samplingUpdate.setNumIterations(factorPrediction);
            end
        end
        
        function [factorPrediction, ...
                  factorUpdate] = getNumSamplesFactors(obj)
            % Get the linear factors to determine the number of samples used for state prediction and measurement upate.
            %
            % Returns:
            %   << factorPrediction (Positive scalar)
            %      The linear factor to determine the number of samples
            %      used for the state prediction.
            %
            %   << factorUpdate (Positive scalar)
            %      The linear factor to determine the number of samples
            %      used for the measurement update.
            
            factorPrediction = obj.samplingPrediction.getNumIterations();
            factorUpdate     = obj.samplingUpdate.getNumIterations();
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

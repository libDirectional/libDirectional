
classdef UnscentedLinearGaussianFilter < SampleBasedLinearGaussianFilter
    % Abstract base class for linear Gaussian filters that are based on the
    % unscented Gaussian sampling technique.
    %
    % UnscentedLinearGaussianFilter Methods:
    %   UnscentedLinearGaussianFilter - Class constructor.
    %   copy                          - Copy a Filter instance.
    %   copyWithName                  - Copy a Filter instance and give the copy a new name/description.
    %   getName                       - Get the filter name/description.
    %   setColor                      - Set the filter color/plotting properties.
    %   getColor                      - Get the filter color/plotting properties.
    %   setState                      - Set the system state.
    %   getState                      - Get the system state.
    %   setStateMeanAndCov            - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov            - Get mean and covariance matrix of the system state.
    %   getStateDim                   - Get the dimension of the system state.
    %   predict                       - Perform a state prediction.
    %   update                        - Perform a measurement update.
    %   step                          - Perform a combined state prediction and measurement update.
    %   setStateDecompDim             - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim             - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing   - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing   - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing       - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing       - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold        - Set the measurement gating threshold.
    %   getMeasGatingThreshold        - Get the measurement gating threshold.
    %   setSampleScalings             - Set the sample scaling factors used for state prediction and measurement update.
    %   getSampleScalings             - Get the sample scaling factors used for state prediction and measurement update.
    
    % Literature:
    %   Simon J. Julier and Jeffrey K. Uhlmann,
    %   Unscented Filtering and Nonlinear Estimation,
    %   Proceedings of the IEEE, vol. 92 no. 3, pp. 401-422, 2004.
    
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
        function obj = UnscentedLinearGaussianFilter(name)
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
            %   << obj (UnscentedLinearGaussianFilter)
            %      A new UnscentedLinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@SampleBasedLinearGaussianFilter(name);
            
            obj.samplingPrediction = GaussianSamplingUKF();
            obj.samplingUpdate     = GaussianSamplingUKF();
            
            % By default, all samples are equally weighted for
            % state prediction and measurement update.
            obj.samplingPrediction.setSampleScaling(0.5);
            obj.samplingUpdate.setSampleScaling(0.5);
        end
        
        function setSampleScalings(obj, scalingPrediction, scalingUpdate)
            % Set the sample scaling factors used for state prediction and measurement update.
            %
            % For example, a scaling factor of 0.5 results in an equal sample
            % weight for all samples, a factor of 1 results in a double
            % weighted sample located at the state space origin, and a factor
            % of 0 results in a zero weight for the sample located at the state
            % space origin.
            %
            % Note: a valid sampling requires a scaling factor larger than -N,
            % where N denotes the requested dimension of the samples.
            %
            % By default, the sample scaling factor is set to 0.5 for both
            % state prediction and measurement update.
            %
            % Parameters:
            %   >> scalingPrediction (Scalar)
            %      The new sample scaling factor used for the state prediction.
            %
            %   >> scalingUpdate (Scalar)
            %      The new sample scaling factor used for the measurement update.
            %      If nothing is passed, the sample scaling factor specified for
            %      the state prediction is also used for the measurement update.
            
            obj.samplingPrediction.setSampleScaling(scalingPrediction);
            
            if nargin == 3
                obj.samplingUpdate.setSampleScaling(scalingUpdate);
            else
                obj.samplingUpdate.setSampleScaling(scalingPrediction);
            end
        end
        
        function [scalingPrediction, ...
                  scalingUpdate] = getSampleScalings(obj)
            % Get the sample scaling factors used for state prediction and measurement update.
            %
            % Returns:
            %   << scalingPrediction (Scalar)
            %      The sample scaling factor for the prediction.
            %
            %   << scalingUpdate (Scalar)
            %      The sample scaling factor for the update.
            
            scalingPrediction = obj.samplingPrediction.getSampleScaling();
            scalingUpdate     = obj.samplingUpdate.getSampleScaling();
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


classdef LinearGaussianFilter < GaussianFilter
    % Abstract base class for linear Gaussian filters.
    %
    % LinearGaussianFilter Methods:
    %   LinearGaussianFilter        - Class constructor.
    %   copy                        - Copy a Filter instance.
    %   copyWithName                - Copy a Filter instance and give the copy a new name/description.
    %   getName                     - Get the filter name/description.
    %   setColor                    - Set the filter color/plotting properties.
    %   getColor                    - Get the filter color/plotting properties.
    %   setState                    - Set the system state.
    %   getState                    - Get the system state.
    %   setStateMeanAndCov          - Set the system state by means of mean and covariance matrix.
    %   getStateMeanAndCov          - Get mean and covariance matrix of the system state.
    %   getStateDim                 - Get the dimension of the system state.
    %   predict                     - Perform a state prediction.
    %   update                      - Perform a measurement update.
    %   step                        - Perform a combined state prediction and measurement update.
    %   setStateDecompDim           - Set the dimension of the unobservable part of the system state.
    %   getStateDecompDim           - Get the dimension of the unobservable part of the system state.
    %   setPredictionPostProcessing - Set a post-processing method for the state prediction.
    %   getPredictionPostProcessing - Get the post-processing method for the state prediction.
    %   setUpdatePostProcessing     - Set a post-processing method for the measurement update.
    %   getUpdatePostProcessing     - Get the post-processing method for the measurement update.
    %   setMeasGatingThreshold      - Set the measurement gating threshold.
    %   getMeasGatingThreshold      - Get the measurement gating threshold.
    
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
        function obj = LinearGaussianFilter(name)
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
            %   << obj (LinearGaussianFilter)
            %      A new LinearGaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@GaussianFilter(name);
            
            % By default, measurement gating is disabled
            obj.measGatingThreshold = 1;
        end
        
        function setMeasGatingThreshold(obj, threshold)
            % Set the measurement gating threshold.
            %
            % A value of zero means that only measurements, which exactly match the expected
            % measurement mean, will be accepted. In contrast, a value of one means that all
            % measurements will be accepted (i.e., no gating will be performed at all).
            %
            % By default, the gating threshold is set to 1 (i.e., no gating).
            %
            % The gating is based on the Normalized Innovation Squared (NIS) as defined in
            %   Yaakov Bar-Shalom, X. Rong Li, and Thiagalingam Kirubarajan,
            %   Estimation with Applications to Tracking and Navigation,
            %   Wiley-Interscience, 2001, page 236
            %
            % Parameters:
            %   >> threshold (Scalar in [0, 1])
            %      The new measurement gating threshold.
            
            if ~Checks.isScalarIn(threshold, 0, 1)
                obj.error('InvalidMeasurementGatingThreshold', ...
                          'threshold must be a scalar in [0, 1].');
            end
            
            obj.measGatingThreshold = threshold;
        end
        
        function threshold = getMeasGatingThreshold(obj)
            % Get the measurement gating threshold.
            %
            % Returns:
            %   << threshold (Scalar in [0, 1])
            %      The measurement gating threshold.
            
            threshold = obj.measGatingThreshold;
        end
    end
    
    methods (Sealed, Access = 'protected')
        function [updatedStateMean, ...
                  updatedStateCov] = performUpdateObservable(obj, measModel, measurement, ...
                                                             priorStateMean, priorStateCov, priorStateCovSqrt)
            obj.checkMeasurementVector(measurement);
            
            dimMeas = size(measurement, 1);
            
            if Checks.isClass(measModel, 'LinearMeasurementModel')
                [updatedStateMean, ...
                 updatedStateCov] = obj.updateLinearMeasModel(measModel, measurement, ...
                                                              priorStateMean, priorStateCov, priorStateCovSqrt);
            else
                if Checks.isClass(measModel, 'MeasurementModel')
                    obj.setupMeasModel(measModel, dimMeas);
                elseif Checks.isClass(measModel, 'AdditiveNoiseMeasurementModel')
                    obj.setupAddNoiseMeasModel(measModel, dimMeas);
                elseif Checks.isClass(measModel, 'MixedNoiseMeasurementModel')
                    obj.setupMixedNoiseMeasModel(measModel, dimMeas);
                else
                    obj.errorMeasModel('LinearMeasurementModel', ...
                                       'MeasurementModel', ...
                                       'AdditiveNoiseMeausrementModel', ...
                                       'MixedNoiseMeasurementModel');
                end
                
                [updatedStateMean, ...
                 updatedStateCov] = obj.updateNonlinear(measurement, ...
                                                        priorStateMean, priorStateCov, priorStateCovSqrt);
            end
        end
        
        function [updatedStateMean, ...
                  updatedStateCov] = updateLinearMeasModel(obj, measModel, measurement, ...
                                                           priorStateMean, priorStateCov, priorStateCovSqrt)
            [measMean, measCov, ...
             stateMeasCrossCov] = measModel.analyticMoments(priorStateMean, priorStateCov, priorStateCovSqrt);
            
            try
                if obj.isMeasGatingEnabled()
                    dimMeas = size(measurement, 1);
                    
                    [updatedStateMean, ...
                     updatedStateCov, ...
                     sqMeasMahalDist] = Utils.kalmanUpdate(priorStateMean, priorStateCov, measurement, ...
                                                           measMean, measCov, stateMeasCrossCov);
                    
                    obj.measurementGating(dimMeas, sqMeasMahalDist);
                else
                    [updatedStateMean, ...
                     updatedStateCov] = Utils.kalmanUpdate(priorStateMean, priorStateCov, measurement, ...
                                                           measMean, measCov, stateMeasCrossCov);
                end
            catch ex
                obj.ignoreMeas(ex.message);
            end
        end
        
        function isEnabled = isMeasGatingEnabled(obj)
            isEnabled = obj.measGatingThreshold ~= 1;
        end
        
        function measurementGating(obj, dimMeas, sqMeasMahalDist)
            % Normalize Mahalanobis distance
            normalizedValue = chi2cdf(sqMeasMahalDist, dimMeas);
            
            % Check for threshold exceedance
            if normalizedValue > obj.measGatingThreshold
                % Discard measurement
               error('Measurement exceeds gating threshold (value: %f).', ...
                     normalizedValue);
            end
        end
    end
    
    methods (Abstract, Access = 'protected')
        setupMeasModel(obj, measModel, dimMeas);
        
        setupAddNoiseMeasModel(obj, measModel, dimMeas);
        
        setupMixedNoiseMeasModel(obj, measModel, dimMeas);
        
        [updatedStateMean, ...
         updatedStateCov] = updateNonlinear(obj, measurement, ...
                                            priorStateMean, priorStateCov, priorStateCovSqrt);
    end
    
    properties (Access = 'private')
        % Measurement gating threshold in [0, 1].
        measGatingThreshold;
    end
end


classdef RecursiveUpdateFilter < LinearGaussianFilter
    % Abstract base class for recursive update filters.
    %
    % RecursiveUpdateFilter Methods:
    %   RecursiveUpdateFilter       - Class constructor.
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
    %   setNumRecursionSteps        - Set the number of recursion steps that are performed by a measurement update.
    %   getNumRecursionSteps        - Get the number of recursion steps that are performed by a measurement update.
    
    % Literature:
    %   Renato Zanetti,
    %   Recursive Update Filtering for Nonlinear Estimation,
    %   IEEE Transactions on Automatic Control, vol. 57, no. 6, pp. 1481-1490, Jun. 2012.
    
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
        function obj = RecursiveUpdateFilter(name)
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
            %   << obj (RecursiveUpdateFilter)
            %      A new RecursiveUpdateFilter instance.
            
            % Call superclass constructor
            obj = obj@LinearGaussianFilter(name);
            
            % By default, 10 recursion steps are performed by a measurement update.
            obj.numRecursionSteps = 10;
        end
        
        function setNumRecursionSteps(obj, numRecursionSteps)
            % Set the number of recursion steps that are performed by a measurement update.
            %
            % By default, 10 recursion steps are performed by a measurement update.
            %
            % Parameters:
            %   >> numRecursionSteps (Positive scalar)
            %      The new number of recursion steps that are performed by a measurement update.
            
            if ~Checks.isPosScalar(numRecursionSteps)
                obj.error('InvalidNumberOfRecursionSteps', ...
                          'numRecursionSteps must be a positive scalar.');
            end
            
            obj.numRecursionSteps = ceil(numRecursionSteps);
        end
        
        function numRecursionSteps = getNumRecursionSteps(obj)
            % Get the number of recursion steps that are performed by a measurement update.
            %
            % Returns:
            %   << numRecursionSteps (Scalar)
            %      The number of recursion steps that are performed by a measurement update.
            
            numRecursionSteps = obj.numRecursionSteps;
        end
    end
    
    methods (Sealed, Access = 'protected')
        function [updatedStateMean, ...
                  updatedStateCov] = updateNonlinear(obj, measurement, ...
                                                     priorStateMean, priorStateCov, priorStateCovSqrt)
            % First recursion step
            [measMean, measCov, ...
             stateMeasCrossCov, ...
             measNoiseCrossCov] = obj.getMeasMoments(priorStateMean, priorStateCov, priorStateCovSqrt);
            
            measCovSqrt = obj.checkCovUpdate(measCov, 'Measurement');
            
            r = 1 / obj.numRecursionSteps;
            A = stateMeasCrossCov / measCovSqrt';
            K = r * (A / measCovSqrt);
            
            innovation                = measurement - measMean;
            updatedStateMean          = priorStateMean + K * innovation;
            updatedStateCov           = priorStateCov + ((r - 2) * r) * (A * A');
            updatedStateNoiseCrossCov = -K * measNoiseCrossCov;
            
            if obj.isMeasGatingEnabled()
                % Perform measurement gating before any measurement information is processed
                dimMeas = size(measurement, 1);
                
                t = measCovSqrt \ innovation;
                sqMeasMahalDist = t' * t;
                
                try
                    obj.measurementGating(dimMeas, sqMeasMahalDist);
                catch ex
                    obj.ignoreMeas(ex.message);
                end
            end
            
            % Next recursion steps
            for i = 2:obj.numRecursionSteps
                % Check if intermediate state covariance matrix is valid
                updatedStateCovSqrt = obj.checkCovUpdate(updatedStateCov, 'Intermediate state');
                
                [measMean, measCov, ...
                 stateMeasCrossCov, ...
                 measNoiseCrossCov] = obj.getMeasMomentsCorr(updatedStateMean, updatedStateCov, updatedStateCovSqrt, ...
                                                             updatedStateNoiseCrossCov);
                
                measCovSqrt = obj.checkCovUpdate(measCov, 'Measurement');
                
                r = 1 / (obj.numRecursionSteps + 1 - i);
                A = stateMeasCrossCov / measCovSqrt';
                K = r * (A / measCovSqrt);
                
                innovation                 = measurement - measMean;
                updatedStateMean          = updatedStateMean + K * innovation;
                updatedStateCov           = updatedStateCov + ((r - 2) * r) * (A * A');
                updatedStateNoiseCrossCov = updatedStateNoiseCrossCov - K * measNoiseCrossCov;
            end
        end
    end
    
    methods (Abstract, Access = 'protected')
        [measMean, measCov, ...
         stateMeasCrossCov, ...
         measNoiseCrossCov] = getMeasMoments(obj, priorStateMean, priorStateCov, priorStateCovSqrt);
        
        [measMean, measCov, ...
         stateMeasCrossCov, ...
         measNoiseCrossCov] = getMeasMomentsCorr(obj, updatedStateMean, updatedStateCov, updatedStateCovSqrt, ...
                                                 updatedStateNoiseCrossCov);
    end
    
    properties (Access = 'private')
        % The number of recursion steps that are performed by a measurement update.
        numRecursionSteps;
    end
end

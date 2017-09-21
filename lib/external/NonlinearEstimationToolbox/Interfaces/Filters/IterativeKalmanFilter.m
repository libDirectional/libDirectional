
classdef IterativeKalmanFilter < LinearGaussianFilter
    % Abstract base class for iterative Kalman filters.
    %
    % IterativeKalmanFilter Methods:
    %   IterativeKalmanFilter       - Class constructor.
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
    %   setMaxNumIterations         - Set the maximum number of iterations that will be performed by a measurement update.
    %   getMaxNumIterations         - Get the maximum number of iterations that will be performed by a measurement update.
    %   getNumIterations            - Get number of iterations performed by the last measurement update.
    %   setConvergenceCheck         - Set a convergence check to determine if no further iterations are required.
    %   getConvergenceCheck         - Get the convergence check.
    
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
        function obj = IterativeKalmanFilter(name)
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
            %   << obj (IterativeKalmanFilter)
            %      A new IterativeKalmanFilter instance.
            
            % Call superclass constructor
            obj = obj@LinearGaussianFilter(name);
            
            % By default, only one iteration is performed during a measurement
            % update, i.e., the usual non-iterated Kalman filter update.
            obj.maxNumIterations = 1;
            
            obj.numIterations = 0;
            
            % By default, no check for convergence is performed
            obj.convergenceCheck = [];
        end
        
        function setMaxNumIterations(obj, maxNumIterations)
            % Set the maximum number of iterations that will be performed by a measurement update.
            %
            % By default, only a single iteration is performed by a measurement update,
            % i.e., the usual non-iterated Kalman filter update.
            %
            % Parameters:
            %   >> maxNumIterations (Positive scalar)
            %      The new maximum number of iterations that will be performed by a measurement update.
            
            if ~Checks.isPosScalar(maxNumIterations)
                obj.error('InvalidNumberOfIterations', ...
                          'maxNumIterations must be a positive scalar.');
            end
            
            obj.maxNumIterations = ceil(maxNumIterations);
        end
        
        function maxNumIterations = getMaxNumIterations(obj)
            % Get the maximum number of iterations that will be performed by a measurement update.
            %
            % Returns:
            %   << maxNumIterations (Positive scalar)
            %      The maximum number of iterations that will be performed by a measurement update.
            
            maxNumIterations = obj.maxNumIterations;
        end
        
        function numIterations = getNumIterations(obj)
            % Get number of iterations performed by the last measurement update.
            %
            % Returns:
            %   << numIterations (Scalar)
            %      The number of iterations performed by the last measurement update.
            
            numIterations = obj.numIterations;
        end
        
        function setConvergenceCheck(obj, convergenceCheck)
            % Set a convergence check to determine if no further iterations are required.
            %
            % The convergence check is executed after each iteration and has
            % to be of the form
            %
            %   convReached = convergenceCheck(lastStateMean, lastStateCov, lastStateCovSqrt, ...
            %                                  updatedStateMean, updatedStateCov, updatedStateCovSqrt)
            %
            % For example, a convergence check can be based on the Kulback-Leibler
            % divergence (see Utils.getGaussianKLD()) or the L2 distance
            % (see Utils.getGaussianL2Distance()) between the state estimates of
            % the last and current iteration.
            %
            % Remove a set convergence check by passing an empty matrix.
            %
            % Parameters:
            %   >> convergenceCheck (Function handle or empty matrix)
            %      Function handle that implements the convergence check or an
            %      empty matrix to remove a previously set convergence check.
            
            if ~Checks.isClass(convergenceCheck, 'function_handle') && ...
               ~isempty(convergenceCheck)
                obj.error('InvalidConvergenceCheck', ...
                          'convergenceCheck must be a function handle or an empty matrix.');
            end
            
            obj.convergenceCheck = convergenceCheck;
        end
        
        function convergenceCheck = getConvergenceCheck(obj)
            % Get the convergence check.
            %
            % Returns:
            %   << convergenceCheck (Function handle or empty matrix)
            %      Function handle that implements the convergence check or an
            %      empty matrix if no convergence check is set.
            
            convergenceCheck = obj.convergenceCheck;
        end
    end
    
    methods (Sealed, Access = 'protected')
        function [updatedStateMean, ...
                  updatedStateCov] = updateNonlinear(obj, measurement, ...
                                                     priorStateMean, priorStateCov, priorStateCovSqrt)
            % Initialize iteration counter
            obj.numIterations = 1;
            
            % Perform first measurement update
            [measMean, measCov, ...
             stateMeasCrossCov] = obj.getMeasMoments(priorStateMean, priorStateCov, priorStateCovSqrt);
            
            try
                if obj.isMeasGatingEnabled()
                    % Perform measurement gating before any measurement information is processed
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
            
            % More updates to perform (iterative update)?
            if obj.maxNumIterations > 1
                % Initially, the last estimate is the prior one.
                lastStateMean    = priorStateMean;
                lastStateCov     = priorStateCov;
                lastStateCovSqrt = priorStateCovSqrt;
                
                while obj.numIterations < obj.maxNumIterations
                    % Check if intermediate state covariance matrix is valid
                    updatedStateCovSqrt = obj.checkCovUpdate(updatedStateCov, 'Intermediate state');
                    
                    % Check for convergence
                    if obj.isConvergenceReached(lastStateMean, lastStateCov, lastStateCovSqrt, ...
                                                updatedStateMean, updatedStateCov, updatedStateCovSqrt)
                        % No more iteration is required according to a
                        % possibly set convergence check
                        return;
                    end
                    
                    % Save estimate from last iteration
                    lastStateMean    = updatedStateMean;
                    lastStateCov     = updatedStateCov;
                    lastStateCovSqrt = updatedStateCovSqrt;
                    
                    % Increment iteration counter
                    obj.numIterations = obj.numIterations + 1;
                    
                    % % Perform measurement update for the current iteration
                    [measMean, measCov, ...
                     stateMeasCrossCov] = obj.getMeasMomentsIteration(priorStateMean, priorStateCov, priorStateCovSqrt, ...
                                                                      updatedStateMean, updatedStateCov, updatedStateCovSqrt);
                    
                    try
                        [updatedStateMean, ...
                         updatedStateCov] = Utils.kalmanUpdate(priorStateMean, priorStateCov, measurement, ...
                                                               measMean, measCov, stateMeasCrossCov);
                    catch ex
                        obj.ignoreMeas(ex.message);
                    end
                end
            end
        end
    end
    
    methods (Abstract, Access = 'protected')
        [measMean, measCov, ...
         stateMeasCrossCov] = getMeasMoments(obj, priorStateMean, priorStateCov, priorStateCovSqrt);
        
        [measMean, measCov, ...
         stateMeasCrossCov] = getMeasMomentsIteration(obj, priorStateMean, priorStateCov, priorStateCovSqrt, ...
                                                      updatedStateMean, updatedStateCov, updatedStateCovSqrt);
    end
    
    methods (Access = 'private')
        function convReached = isConvergenceReached(obj, lastStateMean, lastStateCov, lastStateCovSqrt, ...
                                                    updatedStateMean, updatedStateCov, updatedStateCovSqrt)
            if ~isempty(obj.convergenceCheck)
                convReached = obj.convergenceCheck(lastStateMean, lastStateCov, lastStateCovSqrt, ...
                                                   updatedStateMean, updatedStateCov, updatedStateCovSqrt);
                
                if ~Checks.isFlag(convReached)
                    obj.error('InvalidConvergenceResult', ...
                              'Check for convergence must return a boolean scalar.');
                end
            else
                convReached = false;
            end
        end
    end
    
    properties (Access = 'private')
        % The maximum number of iterations that will be performed by a measurement update.
        maxNumIterations;
        
        % The number of iterations performed by the last measurement update.
        numIterations;
        
        % Function handle to check for convergence.
        convergenceCheck;
    end
end

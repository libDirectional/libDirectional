
classdef GaussianFilter < Filter
    % Abstract base class for Gaussian filters.
    %
    % GaussianFilter Methods:
    %   GaussianFilter              - Class constructor.
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
    
    % Literature:
    %   Jannik Steinbring, Antonio Zea, Uwe D. Hanebeck,
    %   Semi-Analytic Progressive Gaussian Filtering,
    %   Proceedings of the 2016 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems (MFI), Baden-Baden, Germany, Sep. 2016.
    %
    %   Tine Lefebvre, Herman Bruyninckx, Joris De Schutter,
    %   Nonlinear Kalman Filtering for Force-Controlled Robot Tasks,
    %   Appendix E: Partial Observation with the Kalman Filter,
    %   ser. Springer Tracts in Advanced Robotics. Berlin Heidelberg: Springer, vol. 19, 2005.
    
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
        function obj = GaussianFilter(name)
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
            %   << obj (GaussianFilter)
            %      A new GaussianFilter instance.
            
            % Call superclass constructor
            obj = obj@Filter(name);
            
            % By default, it is assumed that the entire system state is
            % required by the measurement model/likelihood function.
            obj.stateDecompDim = 0;
            
            % By default, no post-processing is enabled
            obj.predictionPostProcessing = [];
            obj.updatePostProcessing     = [];
        end
        
        function state = getState(obj)
            state = Gaussian(obj.stateMean, obj.stateCov);
        end
        
        function [stateMean, stateCov, stateCovSqrt] = getStateMeanAndCov(obj)
            stateMean    = obj.stateMean;
            stateCov     = obj.stateCov;
            stateCovSqrt = obj.stateCovSqrt;
        end
        
        function setStateDecompDim(obj, dim)
            % Set the dimension of the unobservable part of the system state.
            %
            % Consider a measurement model/likelihood function that only
            % depends on the subspace o of the entire system state [o, u]',
            % called the observable part. Due to the Gaussian state
            % estimate, the filter step can be divided into two parts.
            % First, the filter step of the respective filter is used to
            % only update the state estimate of the observable part o.
            % Second, the updated estimate of the entire state can now be
            % computed in closed-form based on the updated estimate for the
            % subspace o and the prior state estimate.
            %
            % For example, consider a 5D system state [a, b, c, d, e]'. If
            % the dimension of the unobservable part is set to two, the
            % observable part is o = [a, b, c]' and the unobservable part
            % is u = [d, e]'. That is, it is assumed that the unobservable
            % part comprises the last dimensions of the system state.
            %
            % A value of zero means that the entire system state is
            % required by the measurement model/likelihood function (i.e.,
            % the usual case).
            %
            % By default, the dimension of the unobservable part is set to
            % zero.
            %
            % Parameters:
            %   >> dim (Non-negative scalar)
            %      The new dimension of the unobservable part of the system state.
            
            if ~Checks.isNonNegativeScalar(dim)
                obj.error('InvalidDimension', ...
                          'dim must be a non-negative scalar.');
            end
            
            obj.stateDecompDim = ceil(dim);
        end
        
        function dim = getStateDecompDim(obj)
            % Get the dimension of the unobservable part of the system state.
            %
            % Returns:
            %   << dim (Non-negative scalar)
            %      Dimension of the unobservable part of the system state.
            
            dim = obj.stateDecompDim;
        end
        
        function setPredictionPostProcessing(obj, predictionPostProcessing)
            % Set a post-processing method for the state prediction.
            %
            % The post-processing method is executed after each state prediction
            % and has to be of the form
            %
            %   [postProcessedStateMean, ...
            %    postProcessedStateCov] = predictionPostProcessing(predictedStateMean, ...
            %                                                      predictedStateCov, ...
            %                                                      predictedStateCovSqrt)
            %
            % This can be used, for example, to implement a constrained state estimation.
            %
            % Remove a set post-processing method by passing an empty matrix.
            %
            % Parameters:
            %   >> predictionPostProcessing (Function handle or empty matrix)
            %      Function handle that implements the post-processing or
            %      an empty matrix to remove a previously set post-processing
            %      method.
            
            if ~Checks.isClass(predictionPostProcessing, 'function_handle') && ...
               ~isempty(predictionPostProcessing)
                obj.error('InvalidPostPrecessingMethod', ...
                          'predictionPostProcessing must be a function handle or an empty matrix.');
            end
            
            obj.predictionPostProcessing = predictionPostProcessing;
        end
        
        function predictionPostProcessing = getPredictionPostProcessing(obj)
            % Get the post-processing method for the state prediction.
            %
            % Returns:
            %   << predictionPostProcessing (Function handle or empty matrix)
            %      Function handle that implements the post-processing or
            %      an empty matrix if no post-processing method is set.
            
            predictionPostProcessing = obj.predictionPostProcessing;
        end
        
        function setUpdatePostProcessing(obj, updatePostProcessing)
            % Set a post-processing method for the measurement update.
            %
            % The post-processing method is executed after each measurement update
            % and has to be of the form
            %
            %   [postProcessedStateMean, ...
            %    postProcessedStateCov] = updatePostProcessing(updatedStateMean, ...
            %                                                  updatedStateCov, ...
            %                                                  updatedStateCovSqrt)
            %
            % This can be used, for example, to implement a constrained state estimation.
            %
            % Remove a set post-processing method by passing an empty matrix.
            %
            % Parameters:
            %   >> updatePostProcessing (Function handle or empty matrix)
            %      Function handle that implements the post-processing or
            %      an empty matrix to remove a previously set post-processing
            %      method.
            
            if ~Checks.isClass(updatePostProcessing, 'function_handle') && ...
               ~isempty(updatePostProcessing)
                obj.error('InvalidPostPrecessingMethod', ...
                          'updatePostProcessing must be a function handle or an empty matrix.');
            end
            
            obj.updatePostProcessing = updatePostProcessing;
        end
        
        function updatePostProcessing = getUpdatePostProcessing(obj)
            % Get the post-processing method for the measurement update.
            %
            % Returns:
            %   << updatePostProcessing (Function handle or empty matrix)
            %      Function handle that implements the post-processing or
            %      an empty matrix if no post-processing method is set.
            
            updatePostProcessing = obj.updatePostProcessing;
        end
    end
    
    methods (Sealed, Access = 'protected')
        function performSetState(obj, state)
            [obj.stateMean, obj.stateCov, obj.stateCovSqrt] = state.getMeanAndCov();
        end
        
        function performSetStateMeanAndCov(obj, stateMean, stateCov, stateCovSqrt)
            obj.stateMean    = stateMean;
            obj.stateCov     = stateCov;
            obj.stateCovSqrt = stateCovSqrt;
        end
        
        function performPrediction(obj, sysModel)
            if Checks.isClass(sysModel, 'LinearSystemModel')
                [predictedStateMean, ...
                 predictedStateCov] = obj.predictLinearSysModel(sysModel);
            elseif Checks.isClass(sysModel, 'SystemModel')
                [predictedStateMean, ...
                 predictedStateCov] = obj.predictSysModel(sysModel);
            elseif Checks.isClass(sysModel, 'AdditiveNoiseSystemModel')
                [predictedStateMean, ...
                 predictedStateCov] = obj.predictAddNoiseSysModel(sysModel);
            elseif Checks.isClass(sysModel, 'MixedNoiseSystemModel')
                [predictedStateMean, ...
                 predictedStateCov] = obj.predictMixedNoiseSysModel(sysModel);
            else
                obj.errorSysModel('LinearSystemModel', ...
                                  'SystemModel', ...
                                  'AdditiveNoiseSystemModel', ...
                                  'MixedNoiseSystemModel');
            end
            
            obj.checkAndSavePrediction(predictedStateMean, predictedStateCov);
        end
        
        function performUpdate(obj, measModel, measurement)
            observableStateDim = obj.getObservableStateDim();
            
            % Use decomposed state update?
            if observableStateDim < obj.dimState
                % Extract observable part of the system state
                idx     = 1:observableStateDim;
                mean    = obj.stateMean(idx);
                cov     = obj.stateCov(idx, idx);
                covSqrt = obj.stateCovSqrt(idx, idx);
                
                % Update observable state variables
                [updatedMean, ...
                 updatedCov] = obj.performUpdateObservable(measModel, measurement, ...
                                                           mean, cov, covSqrt);
                
                % Check if updated observable state covariance is valid
                updatedCovSqrt = obj.checkCovUpdate(updatedCov, 'Updated observable state');
                
                % Update entire system state
                [updatedStateMean, ...
                 updatedStateCov] = Utils.decomposedStateUpdate(obj.stateMean, obj.stateCov, obj.stateCovSqrt, ...
                                                                updatedMean, updatedCov, updatedCovSqrt);
            else
                % Update entire system state
                [updatedStateMean, ...
                 updatedStateCov] = obj.performUpdateObservable(measModel, measurement, ...
                                                                obj.stateMean, obj.stateCov, obj.stateCovSqrt);
            end
            
            obj.checkAndSaveUpdate(updatedStateMean, updatedStateCov);
        end
        
        function checkAndSavePrediction(obj, predictedStateMean, predictedStateCov)
            % Check if predicted state covariance is valid
            predictedStateCovSqrt = obj.checkCovPrediction(predictedStateCov, ...
                                                           'Predicted state');
            
            if ~isempty(obj.predictionPostProcessing)
                % Perform post-processing
                [predictedStateMean, ...
                 predictedStateCov] = obj.predictionPostProcessing(predictedStateMean, ...
                                                                   predictedStateCov, ...
                                                                   predictedStateCovSqrt);
                
                if ~Checks.isColVec(predictedStateMean, obj.dimState)
                    obj.error('InvalidStateMean', ...
                              ['Post-processed predicted state mean must be a ' ...
                               'column vector of dimension %d.'], ...
                              obj.dimState);
                end
                
                if ~Checks.isSquareMat(predictedStateCov, obj.dimState)
                    obj.error('InvalidStateCovariance', ...
                              ['Post-processed predicted state covariance must be a ' ...
                               'square matrix of dimension %dx%d.'], ...
                              obj.dimState, obj.dimState);
                end
                
                % Check if post-processed state covariance is valid
                predictedStateCovSqrt = obj.checkCovPrediction(predictedStateCov, ...
                                                               'Post-processed predicted state');
            end
            
            % Save predicted state estimate
            obj.stateMean    = predictedStateMean;
            obj.stateCov     = predictedStateCov;
            obj.stateCovSqrt = predictedStateCovSqrt;
        end
        
        function checkAndSaveUpdate(obj, updatedStateMean, updatedStateCov)
            % Check if updated state covariance is valid
            updatedStateCovSqrt = obj.checkCovUpdate(updatedStateCov, ...
                                                     'Updated state');
            
            if ~isempty(obj.updatePostProcessing)
                % Perform post-processing
                [updatedStateMean, ...
                 updatedStateCov] = obj.updatePostProcessing(updatedStateMean, ...
                                                             updatedStateCov, ...
                                                             updatedStateCovSqrt);
                
                if ~Checks.isColVec(updatedStateMean, obj.dimState)
                    obj.error('InvalidStateMean', ...
                              ['Post-processed updated state mean must be a ' ...
                               'column vector of dimension %d.'], ...
                              obj.dimState);
                end
                
                if ~Checks.isSquareMat(updatedStateCov, obj.dimState)
                    obj.error('InvalidStateCovariance', ...
                              ['Post-processed updated state covariance must be a ' ...
                               'square matrix of dimension %dx%d.'], ...
                              obj.dimState, obj.dimState);
                end
                
                % Check if post-processed state covariance is valid
                updatedStateCovSqrt = obj.checkCovUpdate(updatedStateCov, ...
                                                         'Post-processed updated state');
            end
            
            % Save updated state estimate
            obj.stateMean    = updatedStateMean;
            obj.stateCov     = updatedStateCov;
            obj.stateCovSqrt = updatedStateCovSqrt;
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictLinearSysModel(obj, sysModel)
            [predictedStateMean, ...
             predictedStateCov] = sysModel.analyticMoments(obj.stateMean, ...
                                                           obj.stateCov, ...
                                                           obj.stateCovSqrt);
        end
        
        function [predictedStateMean, ...
                  predictedStateCov] = predictMixedNoiseSysModel(obj, sysModel)
            [addNoiseMean, addNoiseCov]  = sysModel.additiveNoise.getMeanAndCov();
            dimAddNoise = size(addNoiseMean, 1);
            
            obj.checkAdditiveSysNoise(dimAddNoise);
            
            [mean, cov] = obj.predictSysModel(sysModel);
            
            % Compute predicted state mean
            predictedStateMean = mean + addNoiseMean;
            
            % Compute predicted state covariance
            predictedStateCov = cov + addNoiseCov;
        end
    end
    
    methods (Abstract, Access = 'protected')
        [predictedStateMean, ...
         predictedStateCov] = predictSysModel(obj, sysModel);
        
        [predictedStateMean, ...
         predictedStateCov] = predictAddNoiseSysModel(obj, sysModel);
        
        [updatedStateMean, ...
         updatedStateCov] = performUpdateObservable(obj, measModel, measurement, ...
                                                    priorStateMean, priorStateCov, priorStateCovSqrt);
    end
    
    methods (Access = 'private')
        function observableStateDim = getObservableStateDim(obj)
            observableStateDim = obj.dimState - obj.stateDecompDim;
            
            % At least one variable of the system state must be "observable"
            if observableStateDim <= 0
                obj.error('InvalidUnobservableStateDimension', ...
                          'At least one variable of the system state must be "observable".');
            end
        end
    end
    
    properties (GetAccess = 'protected', SetAccess = 'private')
        % The system mean vector.
        stateMean;
        
        % The system covariance matrix.
        stateCov;
        
        % Square root of the system covariance matrix.
        stateCovSqrt;
    end
    
    properties (Access = 'private')
        % Dimension of the unobservable part of the system state.
        stateDecompDim;
        
        % Function handle for state prediction post-processing
        predictionPostProcessing;
        
        % Function handle for measurement update post-processing
        updatePostProcessing;
    end
end


classdef GaussianFilter < Filter
    % Abstract base class for Gaussian filters.
    %
    % This type of filter approximates the estimate of the system state as a Gaussian distribution.
    %
    % GaussianFilter Methods:
    %   GaussianFilter   - Class constructor.
    %   getName          - Get the filter name / description.
    %   setColor         - Set the filter color / plotting properties.
    %   getColor         - Get the current filter color / plotting properties.
    %   setState         - Set the system state.
    %   getState         - Get the current system state.
    %   getStateDim      - Get the dimension of the current system state.
    %   predict          - Perform a time update (prediction step).
    %   update           - Perform a measurement update (filter step) using the given measurement(s).
    %   step             - Perform a combined time and measurement update.
    %   getPointEstimate - Get a point estimate of the current system state.
    
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
        end
        
        function setState(obj, state)
            if ~Checks.isClass(state, 'Distribution')
                obj.error('UnsupportedSystemState', ...
                          'state must be a subclass of Distribution.');
            end
            
            obj.dimState = state.getDimension();
            
            [obj.stateMean, obj.stateCov, obj.stateCovSqrt] = state.getMeanAndCovariance();
        end
        
        function state = getState(obj)
            state = Gaussian(obj.stateMean, obj.stateCov);
        end
        
        function [pointEstimate, uncertainty] = getPointEstimate(obj)
            % Get a point estimate of the current system state.
            %
            % Returns:
            %   << pointEstimate (Column vector)
            %      The current state mean.
            %
            %   << uncertainty (Positive definite matrix)
            %      The current state covariance.
            
            pointEstimate = obj.stateMean;
            uncertainty   = obj.stateCov;
        end
    end
    
    methods (Access = 'protected')
        function predictAnalytic(obj, sysModel)
            % Compute predicted state moments
            [predictedStateMean, ...
             predictedStateCov] = sysModel.analyticPredictedMoments(obj.stateMean, ...
                                                                    obj.stateCov);
            
            % Check predicted moments
            covSqrt = obj.checkPredictedMoments(predictedStateMean, predictedStateCov);
            
            % Save new state estimate
            obj.stateMean    = predictedStateMean;
            obj.stateCov     = predictedStateCov;
            obj.stateCovSqrt = covSqrt;
        end
        
        function checkAndSavePrediction(obj, predictedStateMean, predictedStateCov)
            % Check predicted state covariance is valid
            [isPosDef, covSqrt] = Checks.isCov(predictedStateCov);
            
            if ~isPosDef
                obj.warnIgnorePrediction('Predicted state covariance is not positive definite.');
                return;
            end
            
            % Save new state estimate
            obj.stateMean    = predictedStateMean;
            obj.stateCov     = predictedStateCov;
            obj.stateCovSqrt = covSqrt;
        end
    end
    
    methods (Access = 'private')
        function covSqrt = checkPredictedMoments(obj, mean, covariance)
            if ~Checks.isColVec(mean, obj.dimState)
                obj.error('InvalidPredictedStateMean', ...
                          ['Predicted state mean must be a ' ...
                           'column vector of dimension %d.'], ...
                          obj.dimState);
            end
            
            [isPosDef, covSqrt] = Checks.isCov(covariance, obj.dimState);
            
            if ~isPosDef
                obj.error('InvalidPredictedStateCovariance', ...
                          ['Predicted state covariance must be a ' ...
                           'positive definite matrix of dimension %dx%d.'], ...
                           obj.dimState, obj.dimState);
            end
        end
    end
    
    properties (Access = 'protected')
        % Current system mean vector.
        stateMean;
        
        % Current system covariance matrix.
        stateCov;
        
        % Sqrt of the current system covariance matrix.
        stateCovSqrt;
    end
end

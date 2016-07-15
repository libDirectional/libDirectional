
classdef UKF < LRKF
    % The Unscented Kalman Filter (UKF).
    %
    % UKF Methods:
    %   UKF                            - Class constructor.
    %   getName                        - Get the filter name / description.
    %   setColor                       - Set the filter color / plotting properties.
    %   getColor                       - Get the current filter color / plotting properties.
    %   setState                       - Set the system state.
    %   getState                       - Get the current system state.
    %   getStateDim                    - Get the dimension of the current system state.
    %   predict                        - Perform a time update (prediction step).
    %   update                         - Perform a measurement update (filter step) using the given measurement(s).
    %   step                           - Perform a combined time and measurement update.
    %   getPointEstimate               - Get a point estimate of the current system state.
    %   setMaxNumIterations            - Set the maximum number of iterations that will be performed during a measurement update.
    %   getMaxNumIterations            - Get the current maximum number of iterations that will be performed during a measurement update.
    %   setMeasValidationThreshold     - Set a threshold to perform a measurement validation (measurement acceptance/rejection).
    %   getMeasValidationThreshold     - Get the current measurement validation threshold.
    %   getLastUpdateData              - Get information from the last performed measurement update.
    %   setUseAnalyticSystemModel      - Enable or disable the use of analytic moment calculation during a prediction.
    %   getUseAnalyticSystemModel      - Get the current use of analytic moment calculation during a prediction.
    %   setUseAnalyticMeasurementModel - Enable or disable the use of analytic moment calculation during a filter step.
    %   getUseAnalyticMeasurementModel - Get the current use of analytic moment calculation during a filter step.
    %   setSampleScaling               - Set the sample scaling factor.
    %   getSampleScaling               - Get the current sample scaling factor.
    
    % Literature:
    %   Simon J. Julier and Jeffrey K. Uhlmann,
    %   Unscented Filtering and Nonlinear Estimation,
    %   Proceedings of the IEEE volume 92 No. 3, pages 401-422, 2004
    
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
        function obj = UKF(name)
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
            %      Default name: 'UKF'.
            %
            % Returns:
            %   << obj (UKF)
            %      A new UKF instance.
            
            if nargin < 1
                name = 'UKF';
            end
            
            sampling = GaussianSamplingUKF();
            
            % Call superclass constructor
            obj = obj@LRKF(name, sampling);
            
            obj.ukfSampling = sampling;
            
            obj.setSampleScaling(0.5);
        end
        
        function setSampleScaling(obj, scaling)
            % Set the sample scaling factor.
            % 
            % For example, a scaling factor of 0.5 results in an equal sample
            % weight for all samples, a factor of 1 results in a double
            % weighted sample located at the state space origin, and a factor
            % of 0 results in a zero weight for the sample located at the state
            % space origin.
            %
            % By default, the sample scaling factor is set to 0.5.
            %
            % Parameters:
            %   >> scaling (Non-negative scalar)
            %      The new sample scaling factor.
            
            obj.ukfSampling.setSampleScaling(scaling);
        end
        
        function scaling = getSampleScaling(obj)
            % Get the current sample scaling factor.
            % 
            % Returns:
            %   << scaling (Non-negative scalar)
            %      The current sample scaling factor.
            
            scaling = obj.ukfSampling.getSampleScaling();
        end
    end
    
    properties (Access = 'private')
        % Gaussian sampling used for prediction and update.
        ukfSampling;
    end
end

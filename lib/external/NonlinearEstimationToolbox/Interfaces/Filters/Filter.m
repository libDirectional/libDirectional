
classdef Filter < handle & matlab.mixin.Copyable
    % Abstract base class for a general (nonlinear) filter.
    %
    % Filter Methods:
    %   Filter           - Class constructor.
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
        function obj = Filter(name)
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
            %   << obj (Filter)
            %      A new Filter instance.
            
            if ~ischar(name)
                error('Filter:InvalidFilterName', ...
                      'name must be a char.');
            end
            
            % Set filter name
            obj.name = name;
            
            % Set default filter color
            obj.color = 'b';
            
            % Initially, there is no valid state.
            obj.dimState = 0;
        end
        
        function filterName = getName(obj)
            % Get the filter name / description.
            %
            % Returns:
            %   << filterName (Char)
            %      The filter name / description.
            
            filterName = obj.name;
        end
        
        function setColor(obj, color)
            % Set the filter color / plotting properties.
            %
            % Assign a color (e.g., 'r') or a color with additional line
            % plotting settings (e.g., { '--', 'LineWidth', 2, 'Color',
            % [0 0.5 0] }) to the filter.
            %
            % Inteded / useful for plotting purposes.
            %
            % Parameters:
            %   >> color (Arbitrary data)
            %      The new filter color / plotting properties.
            
            obj.color = color;
        end
        
        function color = getColor(obj)
            % Get the current filter color / plotting properties.
            %
            % Returns:
            %   << color (Arbitrary data)
            %      The current filter color / plotting properties.
            
            color = obj.color;
        end
        
        function dim = getStateDim(obj)
            % Get the dimension of the current system state.
            %
            % Returns:
            %   << dim (Scalar)
            %      The dimension of the current system state.
            
            dim = obj.dimState;
        end
        
        function runtime = predict(obj, sysModel)
            % Perform a time update (prediction step).
            %
            % Parameters:
            %   >> sysModel (Arbitrary class (filter dependent))
            %      System model that provides the mapping between the prior system
            %      state and the predicted state (i.e., the system state's temporal evolution).
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the prediction step.
            
            if nargout == 1
                s = tic;
                obj.performPrediction(sysModel);
                runtime = toc(s);
            else
                obj.performPrediction(sysModel);
            end
        end
        
        function runtime = update(obj, measModel, measurements)
            % Perform a measurement update (filter step) using the given measurement(s).
            %
            % Parameters:
            %   >> measModel (Arbitrary class (filter dependent))
            %      Measurement model that provides the mapping between measurements and the system state.
            %
            %   >> measurements (Matrix)
            %      Column-wise arranged measurement vectors, where each column represents an individual
            %      measurement. In case of two or more measurements (i.e., two or more columns), the
            %      filter assumes that the measurements originate from the same measurement model and
            %      i.i.d. measurement noise. For example, in case of a measurement model h(x, v) and two
            %      measurements m1 and m2 the filter assumes
            %
            %          m1 = h(x, v) and m2 = h(x, v) .
            %
            %      The advantage is that one has to set the measurement noise v only for one
            %      measurement, no matter how many measurements will be provided in one filter step.
            %      That is, the measurement noise is assumed to be i.i.d. for all measurements.
            %      However, in case of non-i.i.d. measurement noise 
            %
            %          m1 = h(x, v1) and m2 = h(x, v2)
            %
            %      (e.g., existing correlations between noise for different measurements or in general
            %      different noise for different measurements) one has to explicitly stack the
            %      measurement noise to v = [v1; v2] and pass the measurements m1 and m2 as a stacked
            %      measurement vector m = [m1; m2].
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the measurement update.
            
            obj.checkMeasurements(measurements);
            
            if nargout == 1
                s = tic;
                obj.performUpdate(measModel, measurements);
                runtime = toc(s);
            else
                obj.performUpdate(measModel, measurements);
            end
        end
        
        function runtime = step(obj, sysModel, measModel, measurements)
            % Perform a combined time and measurement update.
            %
            % By default, this is equal to execute predict() followed by an update().
            %
            % However, each filter can overwrite this behavior with a custom implementation
            % of a combined time and measurement update in order to improve estimation quality
            % (e.g., done by the Auxiliary SIR-PF).
            %
            % Parameters:
            %   >> sysModel (Arbitrary class (filter dependent))
            %      System model that provides the mapping between the prior system
            %      state and the predicted state (i.e., the system state's temporal evolution).
            %
            %   >> measModel (Arbitrary class (filter dependent))
            %      Measurement model that provides the mapping between measurements and the system state.
            %
            %   >> measurements (Matrix)
            %      Column-wise arranged measurement vectors, where each column represents an individual
            %      measurement. In case of two or more measurements (i.e., two or more columns), the
            %      filter assumes that the measurements originate from the same measurement model and
            %      i.i.d. measurement noise. For example, in case of a measurement model h(x, v) and two
            %      measurements m1 and m2 the filter assumes
            %
            %          m1 = h(x, v) and m2 = h(x, v) .
            %
            %      The advantage is that one has to set the measurement noise v only for one
            %      measurement, no matter how many measurements will be provided in one filter step.
            %      That is, the measurement noise is assumed to be i.i.d. for all measurements.
            %      However, in case of non-i.i.d. measurement noise
            %
            %          m1 = h(x, v1) and m2 = h(x, v2)
            %
            %      (e.g., existing correlations between noise for different measurements or in general
            %      different noise for different measurements) one has to explicitly stack the
            %      measurement noise to v = [v1; v2] and pass the measurements m1 and m2 as a stacked
            %      measurement vector m = [m1; m2].
            %
            % Returns:
            %   << runtime (Scalar)
            %      Time needed to perform the combined time and measurement update.
            
            obj.checkMeasurements(measurements);
            
            if nargout == 1
                s = tic;
                obj.performStep(sysModel, measModel, measurements);
                runtime = toc(s);
            else
                obj.performStep(sysModel, measModel, measurements);
            end
        end
    end
    
    methods (Abstract)
        % Set the system state.
        %
        % This function is mainly used to set an initial system state, as
        % it is intended that the Filter is responsible for modifying the
        % system state by exploiting system and measurement models.
        %
        % Parameters:
        %   >> state (Subclass of Distribution)
        %      The new system state.
        setState(obj, state);
        
        % Get the current system state.
        %
        % Returns:
        %   << state (Subclass of Distribution)
        %      The current system state.
        state = getState(obj);
        
        % Get a point estimate of the current system state.
        %
        % Returns:
        %   << pointEstimate (Column vector)
        %      Point estimate of the current system state.
        %
        %   << uncertainty (Positive definite matrix)
        %      Uncertainty of the current system state point estimate (e.g.,
        %      a covariance matrix).
        [pointEstimate, uncertainty] = getPointEstimate(obj);
    end
    
    methods (Abstract, Access = 'protected')
        performPrediction(obj, sysModel);
        
        performUpdate(obj, measModel, measurements);
    end
    
    methods (Access = 'protected')
        function performStep(obj, sysModel, measModel, measurements)
            obj.performPrediction(sysModel);
            obj.performUpdate(measModel, measurements);
        end
        
        function checkMeasurements(obj, measurements)
            if ~Checks.isMat(measurements)
                obj.error('InvalidMeasurements', ...
                          'measurements must be a matrix.');
            end
        end
        
        function checkPredictedStateSamples(obj, samples, numSamples)
            if ~Checks.isMat(samples, obj.dimState, numSamples)
                obj.error('InvalidPredictedStateSamples', ...
                          ['Predicted state samples have to be stored as a ' ...
                           'matrix of dimension %dx%d.'], ...
                          obj.dimState, numSamples);
            end
            
            areFinite = isfinite(samples);
            
            if ~all(areFinite(:))
                obj.error('InvalidPredictedStateSamples', ...
                          'At least one predicted state sample contains NaN or Inf.');
            end
        end
        
        function checkComputedMeasurements(obj, measurements, dimMeas, numMeas)
            if ~Checks.isMat(measurements, dimMeas, numMeas)
                obj.error('InvalidMeasurements', ...
                          ['Computed measurements have to be stored as a ' ...
                           'matrix of dimension %dx%d.'], ...
                          dimMeas, numMeas);
            end
            
            areFinite = isfinite(measurements);
            
            if ~all(areFinite(:))
                obj.error('InvalidMeasurements', ...
                          'At least one measurement contains NaN or Inf.');
            end
        end
        
        function checkLogLikelihoodEvaluations(obj, values, numEvaluations)
            if ~Checks.isRowVec(values, numEvaluations)
                obj.error('InvalidLogLikelihoodEvaluations', ...
                          ['Logarithmic likelihood evaluations have to be stored as a ' ...
                           'row vector of length %d.'], ...
                          numEvaluations);
            end
            
            if any(isnan(values)) || any(values > 0 & isinf(values))
                obj.error('InvalidLogLikelihoodEvaluations', ...
                          'At least one logarithmic likelihood evaluation is NaN or +Inf.');
            end
        end
        
        function checkAdditiveSysNoise(obj, dimNoise)
            if dimNoise ~= obj.dimState
                obj.error('InvalidAdditiveSystemNoise', ...
                          'System state and additive system noise with incompatible dimensions.');
            end
        end
        
        function checkAdditiveMeasNoise(obj, dimMeas, dimNoise)
            if dimNoise ~= dimMeas
                obj.error('InvalidAdditiveMeasurementNoise', ...
                          'Measurement and additive measurement noise with incompatible dimensions.');
            end
        end
        
        function warning(obj, id, msg, varargin)
            msg = sprintf(msg, varargin{:});
            
            warning(sprintf('Filter:%s', id), ...
                    'From "%s":\n%s', obj.name, msg);
        end
        
        function warnIgnorePrediction(obj, reason, varargin)
            reason = sprintf(reason, varargin{:});
            
            obj.warning('IgnoringPrediction', ...
                        '%s\nIgnoring prediction and leaving state estimate unchanged.', reason);
        end
        
        function warnIgnoreMeas(obj, reason, varargin)
            reason = sprintf(reason, varargin{:});
            
            obj.warning('IgnoringMeasurement', ...
                        '%s\nIgnoring measurement and leaving state estimate unchanged.', reason);
        end
        
        function error(obj, id, msg, varargin)
            msg = sprintf(msg, varargin{:});
            
            error(sprintf('Filter:%s', id), ...
                  'From "%s":\n%s', obj.name, msg);
        end
        
        function errorSysModel(obj, varargin)
            numModels = numel(varargin);
            
            msg = 'Supported system models:';
            
            for i = 1:numModels
                msg = sprintf('%s\n * %s', msg, varargin{i});
            end
            
            obj.error('UnsupportedSystemModel', msg);
        end
        
        function errorMeasModel(obj, varargin)
            numModels = numel(varargin);
            
            msg = 'Supported measurement models:';
            
            for i = 1:numModels
                msg = sprintf('%s\n * %s', msg, varargin{i});
            end
            
            obj.error('UnsupportedMeasurementModel', msg);
        end
    end
    
    properties (Access = 'protected')
        % Dimension of the current system state.
        dimState;
    end
    
    properties (Access = 'private')
        % The filter name / description.
        name;
        
        % The filter color / plotting properties.
        color;
    end
end

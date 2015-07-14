
classdef AdditiveNoiseMeasurementModel < Likelihood
    % Abstract base class for measurement models corrupted by additive noise.
    %
    % AdditiveNoiseMeasurementModel Methods:
    %   setNoise            - Set the measurement noise.
    %   measurementEquation - The measurement equation.
    %   logLikelihood       - Evaluate the logarithmic likelihood function of the implemented measurement equation.
    %   derivative          - Compute the derivative of the implemented measurement equation.
    %   simulate            - Simulate one ore more measurements for a given system state.
    
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
        function setNoise(obj, noise)
            % Set the measurement noise.
            %
            % Parameters:
            %   >> noise (Subclass of Distribution or cell array containing subclasses of Distribution)
            %      The new measurement noise.
            
            if Checks.isClass(noise, 'Distribution')
                obj.noise = noise;
            else
                obj.noise = JointDistribution(noise);
            end
        end
        
        function logValues = logLikelihood(obj, stateSamples, measurements)
            % Evaluate the logarithmic likelihood function of the implemented measurement equation.
            %
            % Parameters:
            %   >> stateSamples (Matrix)
            %      L column-wise arranged state samples.
            %
            %   >> measurements (Matrix)
            %      Column-wise arranged measurement vectors, where each column represents an
            %      individual measurement.
            %
            % Returns:
            %   << logValues (Row vector)
            %      L column-wise arranged logarithmic likelihood function values.
            
            [dimMeas, numMeas] = size(measurements);
            dimNoise   = obj.noise.getDimension();
            numSamples = size(stateSamples, 2); 
            
            if dimMeas ~= dimNoise
                error('AdditiveNoiseMeasurementModel:InvalidMeasurementNoise', ...
                      'Measurement and additive measurement noise with different dimensions.');
            end
            
            % Evaluate deterministic measurement equation
            deterministicEvals = obj.measurementEquation(stateSamples);
            
            % Check computed deterministic measurements
            if ~Checks.isMat(deterministicEvals, dimMeas, numSamples)
                error('AdditiveNoiseMeasurementModel:InvalidMeasurements', ...
                      ['Computed deterministic measurements have to be ' ...
                       'stored as a matrix of dimension %dx%d.'], ...
                       dimMeas, numSamples);
            end
            
            logValues = zeros(1, numSamples);
            
            for i = 1:numMeas
                values = bsxfun(@minus, measurements(:, i), deterministicEvals);
                
                logV = obj.noise.logPdf(values);
                
                logValues = logValues + logV;
            end
        end
        
        function stateJacobian = derivative(obj, nominalState)
            % Compute the derivative of the implemented measurement equation.
            %
            % By default, the Jacobians are computed using a difference quotient.
            %
            % Mainly used by the EKF.
            %
            % Parameters:
            %   >> nominalState (Column vector)
            %      The nominal system state vector to linearize the measurement equation.
            %
            % Returns:
            %   << stateJacobian (Square matrix)
            %      The Jacobian of the state variables.
            
            stateJacobian = Utils.diffQuotientState(@(s) obj.measurementEquation(s), ...
                                                    nominalState);
        end
        
        function measurements = simulate(obj, state, numMeasurements)
            % Simulate one ore more measurements for a given system state.
            %
            % Parameters:
            %   >> state (Column vector)
            %      The system state.
            %
            %   >> numMeasurements (Positive scalar)
            %      The number of measurements to simulate.
            %      Default: One measurement will be simulated.
            %
            % Returns:
            %   << measurements (Matrix)
            %      Column-wise arranged simulated measurements.
            
            if ~Checks.isColVec(state)
                error('AdditiveNoiseMeasurementModel:InvalidSystemState', ...
                      'state must be a column vector.');
            end
            
            if nargin < 3
                numMeasurements = 1;
            else
                if ~Checks.isPosScalar(numMeasurements)
                    error('AdditiveNoiseMeasurementModel:InvalidSystemState', ...
                          'numMeasurements must be a positive scalar.');
                end
                
                numMeasurements = ceil(numMeasurements);
            end
            
            noiseSamples = obj.noise.drawRndSamples(numMeasurements);
            
            measurements = obj.measurementEquation(state);
            measurements = repmat(measurements, 1, numMeasurements) + noiseSamples;
        end
    end
    
    methods (Abstract)
        % The measurement equation.
        %
        % Parameters:
        %   >> stateSamples (Matrix)
        %      L column-wise arranged state samples.
        %
        % Returns:
        %   << measurements (Matrix)
        %      L column-wise arranged measurement samples.
        measurements = measurementEquation(obj, stateSamples);
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        % The measurement noise.
        noise;
    end
end

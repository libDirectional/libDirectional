
classdef MixedNoiseMeasurementModel < handle
    % Abstract base class for system models corrupted by a mixture additive and arbitrary noise.
    %
    % MixedNoiseMeasurementModel Methods:
    %   setNoise            - Set the measurement noise.
    %   setAdditiveNoise    - Set the pure additive measurement noise.
    %   measurementEquation - The measurement equation.
    %   derivative          - Compute the first-order and second-order derivatives of the implemented measurement equation.
    %   simulate            - Simulate a measurement for the given system state.
    
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
    
    methods
        function setNoise(obj, noise)
            % Set the measurement noise.
            %
            % Parameters:
            %   >> noise (Subclass of Distribution)
            %      The new measurement noise.
            
            if Checks.isClass(noise, 'Distribution')
                obj.noise = noise;
            else
                error('MixedNoiseMeasurementModel:InvalidNoise', ...
                      'noise must be a subclass of Distribution.');
            end
        end
        
        function setAdditiveNoise(obj, noise)
            % Set the pure additive measurement noise.
            %
            % Parameters:
            %   >> noise (Subclass of Distribution)
            %      The new pure additive measurement noise.
            
            if Checks.isClass(noise, 'Distribution')
                obj.additiveNoise = noise;
            else
                error('MixedNoiseMeasurementModel:InvalidNoise', ...
                      'noise must be a subclass of Distribution.');
            end
        end
        
        function [stateJacobian, noiseJacobian, ...
                  stateHessians, noiseHessians] = derivative(obj, nominalState, nominalNoise)
            % Compute the first-order and second-order derivatives of the implemented measurement equation.
            %
            % By default, the derivatives are computed using difference quotients.
            %
            % Mainly used by EKF and EKF2.
            %
            % Parameters:
            %   >> nominalState (Column vector)
            %      The nominal system state vector.
            %
            %   >> nominalNoise (Column vector)
            %      The nominal system noise vector.
            %
            % Returns:
            %   << stateJacobian (Square matrix)
            %      The Jacobian of the state variables.
            %
            %   << noiseJacobian (Matrix)
            %      The Jacobian of the noise variables.
            %
            %   << stateHessians (3D matrix)
            %      The Hessians of the state variables.
            %
            %   << noiseHessians (3D matrix)
            %      The Hessians of the noise variables.
            
            if nargout <= 2
                [stateJacobian, noiseJacobian] = Utils.diffQuotientStateAndNoise(@obj.measurementEquation, ...
                                                                                 nominalState, nominalNoise);
            else
                [stateJacobian, noiseJacobian, ...
                 stateHessians, noiseHessians] = Utils.diffQuotientStateAndNoise(@obj.measurementEquation, ...
                                                                                 nominalState, nominalNoise);
            end
        end
        
        function measurement = simulate(obj, state)
            % Simulate a measurement for the given system state.
            %
            % Parameters:
            %   >> state (Column vector)
            %      The system state.
            %
            % Returns:
            %   << measurement (Column vector)
            %      The simulated measurement.
            
            if ~Checks.isColVec(state)
                error('MixedNoiseMeasurementModel:InvalidSystemState', ...
                      'state must be a column vector.');
            end
            
            addNoiseSamples = obj.additiveNoise.drawRndSamples(1);
            noiseSamples    = obj.noise.drawRndSamples(1);
            
            measurement = obj.measurementEquation(state, noiseSamples) + addNoiseSamples;
        end
    end
    
    methods (Abstract)
        % The measurement equation.
        %
        % Parameters:
        %   >> stateSamples (Matrix)
        %      L column-wise arranged state samples.
        %
        %   >> noiseSamples (Matrix)
        %      L column-wise arranged measurement noise samples.
        %
        % Returns:
        %   << measurements (Matrix)
        %      L column-wise arranged measurement samples.
        measurements = measurementEquation(obj, stateSamples, noiseSamples);
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        % The measurement noise.
        noise;
        
        % The pure additive measurement noise.
        additiveNoise;
    end
end

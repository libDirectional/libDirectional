
classdef AdditiveNoiseSystemModel < handle
    % Abstract base class for system models corrupted by additive noise.
    %
    % AdditiveNoiseSystemModel Methods:
    %   setNoise       - Set the system noise.
    %   systemEquation - The system equation.
    %   derivative     - Compute the first-order and second-order derivatives of the implemented system equation.
    %   simulate       - Simulate the temporal evolution for a given system state.
    
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
            % Set the system noise.
            %
            % Parameters:
            %   >> noise (Subclass of Distribution or cell array containing subclasses of Distribution)
            %      The new system noise.
            
            if Checks.isClass(noise, 'Distribution')
                obj.noise = noise;
            else
                obj.noise = JointDistribution(noise);
            end
        end
        
        function [stateJacobian, stateHessians] = derivative(obj, nominalState)
            % Compute the first-order and second-order derivatives of the implemented system equation.
            %
            % By default, the derivatives are computed using difference quotients.
            %
            % Mainly used by EKF and EKF2.
            %
            % Parameters:
            %   >> nominalState (Column vector)
            %      The nominal system state vector.
            %
            % Returns:
            %   << stateJacobian (Square matrix)
            %      The Jacobian of the state variables.
            %
            %   << stateHessians (3D matrix)
            %      The Hessians of the state variables.
            
            if nargout == 1
                stateJacobian = Utils.diffQuotientState(@obj.systemEquation, nominalState);
            else
                [stateJacobian, ...
                 stateHessians] = Utils.diffQuotientState(@obj.systemEquation, nominalState);
            end
        end
        
        function predictedState = simulate(obj, state)
            % Simulate the temporal evolution for a given system state.
            %
            % Parameters:
            %   >> state (Column vector)
            %      The system state to predict.
            %
            % Returns:
            %   << predictedState (Column vector)
            %      The simulated temporal system state evolution.
            
            if ~Checks.isColVec(state)
                error('AdditiveNoiseSystemModel:InvalidSystemState', ...
                      'state must be a column vector.');
            end
            
            noiseSample = obj.noise.drawRndSamples(1);
            
            predictedState = obj.systemEquation(state) + noiseSample;
        end
    end
    
    methods (Abstract)
        % The system equation.
        %
        % Parameters:
        %   >> stateSamples (Matrix)
        %      L column-wise arranged state samples.
        %
        % Returns:
        %   << predictedStates (Matrix)
        %      L column-wise arranged predicted state samples.
        predictedStates = systemEquation(obj, stateSamples);
    end
    
    properties (SetAccess = 'private', GetAccess = 'public')
        % The system noise.
        noise;
    end
end

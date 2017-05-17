
classdef GaussianSamplingUKF < GaussianSampling
    % Implements the UKF Gaussian sampling technique.
    %
    % GaussianSamplingUKF Methods:
    %   GaussianSamplingUKF - Class constructor.
    %   copy                - Copy a GaussianSampling instance.
    %   getStdNormalSamples - Get a set of samples approximating a standard normal distribution.
    %   getSamples          - Get a set of samples approximating a Gaussian distribution.
    %   getSampleScaling    - Get the current sample scaling factor.
    
    % Literature:
    %   Simon J. Julier and Jeffrey K. Uhlmann,
    %   Unscented Filtering and Nonlinear Estimation,
    %   Proceedings of the IEEE, Vol. 92 No. 3, 2004, pp. 401-422.
    
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
        function obj = GaussianSamplingUKF()
            % Class constructor.
            %
            % Returns:
            %   << obj (GaussianSamplingUKF)
            %      A new GaussianSamplingUKF instance.
            
            % Default scaling
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
            % Note: a valid sampling requires a scaling factor larger than -N,
            % where N denotes the requested dimension of the samples.
            %
            % By default, the sample scaling factor is set to 0.5.
            %
            % Parameters:
            %   >> scaling (Scalar)
            %      The new sample scaling factor.
            
            if ~Checks.isScalar(scaling)
                error('GaussianSamplingUKF:InvalidScaling', ...
                      'scaling must be a scalar.');
            end
            
            obj.scaling = scaling;
        end
        
        function scaling = getSampleScaling(obj)
            % Get the current sample scaling factor.
            % 
            % Returns:
            %   << scaling (Scalar)
            %      The current sample scaling factor.
            
            scaling = obj.scaling;
        end
        
        function [samples, weights, numSamples] = getStdNormalSamples(obj, dimension)
            if ~Checks.isPosScalar(dimension)
                error('GaussianSamplingUKF:InvalidDimension', ...
                      'dimension must be a positive scalar.');
            end
            
            if dimension + obj.scaling <= 0
                error('GaussianSamplingUKF:InvalidScaling', ...
                      'The sample scaling factor is too small for the requested dimension.');
            end
            
            numSamples = 2 * dimension + 1;
            
            mat = sqrt(dimension + obj.scaling) * speye(dimension);
            
            samples = [zeros(dimension, 1) -mat mat];
            
            weights = (1 / (2 * (dimension + obj.scaling))) * ...
                      [2 * obj.scaling ones(1, numSamples - 1)];
        end
    end
    
    properties (Access = 'private')
        % Sample scaling factor.
        scaling;
    end
end

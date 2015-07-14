
classdef Distribution < handle
    % Abstract base class for probability distributions.
    %
    % Distribution Methods:
    %   getDimension         - Get the dimension of the distribution.
    %   getMeanAndCovariance - Get mean and covariance of the distribution.
    %   drawRndSamples       - Draw random samples from the distribution.
    %   logPdf               - Evaluate the logarithmic probability density function (pdf) of the distribution.
    
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
    %    aldong with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    methods (Abstract)
        % Get the dimension of the distribution.
        %
        % Returns:
        %   << dim (Positive scalar)
        %      The dimension of the distribution.
        dim = getDimension(obj);
        
        % Get mean and covariance of the distribution.
        %
        % Returns:
        %   << mean (Vector)
        %      Mean of the distribution.
        %
        %   << covariance (Positive definite matrix)
        %      Covariance of the distribution.
        %
        %   << covSqrt (Square matrix)
        %      Square root of the distribution's covariance.
        [mean, covariance, covSqrt] = getMeanAndCovariance(obj);
        
        % Draw random samples from the distribution.
        %
        % Parameters:
        %   >> numSamples (Positive scalar)
        %      Number of samples to draw from the distribution.
        %
        % Returns:
        %   << rndSamples (Matrix)
        %      Column-wise arranged random samples.
        rndSamples = drawRndSamples(obj, numSamples);
        
        % Evaluate the logarithmic probability density function (pdf) of the distribution.
        %
        % Parameters:
        %   >> values (Matrix)
        %      Column-wise arranged vectors of the distribution's dimension.
        %
        % Returns:
        %   << logValues (Row vector)
        %      Column-wise arranged corresponding log-likelihood evaluations.
        logValues = logPdf(obj, values);
    end
    
    methods (Access = 'protected')
        function checkValues(obj, values)
            dim = obj.getDimension();
            
            if ~Checks.isFixedRowMat(values, dim)
                error('Distribution:InvalidValues', ...
                      'values must be a matrix of dimension %dxN.', dim);
            end
        end
    end
end

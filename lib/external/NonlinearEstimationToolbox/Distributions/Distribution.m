
classdef Distribution < handle & matlab.mixin.Copyable
    % Abstract base class for probability distributions.
    %
    % Distribution Methods:
    %   copy           - Copy a distribution instance.
    %   getDim         - Get the dimension of the distribution.
    %   getMeanAndCov  - Get mean and covariance matrix of the distribution.
    %   drawRndSamples - Draw random samples from the distribution.
    %   logPdf         - Evaluate the logarithmic probability density function (PDF) of the distribution.
    
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
    
    methods (Abstract)
        % Get the dimension of the distribution.
        %
        % Returns:
        %   << dim (Positive scalar)
        %      The dimension of the distribution.
        dim = getDim(obj);
        
        % Get mean and covariance matrix of the distribution.
        %
        % Returns:
        %   << mean (Column vector)
        %      Mean of the distribution.
        %
        %   << cov (Positive definite matrix)
        %      Covariance matrix of the distribution.
        %
        %   << covSqrt (Square matrix)
        %      Lower Cholesky decomposition of the distribution's covariance matrix.
        [mean, cov, covSqrt] = getMeanAndCov(obj);
        
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
        
        % Evaluate the logarithmic probability density function (PDF) of the distribution.
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
    
    methods (Sealed, Access = 'protected')
        function checkValues(obj, values)
            dim = obj.getDim();
            
            if ~Checks.isFixedRowMat(values, dim)
                error('Distribution:InvalidValues', ...
                      'values must be a matrix of dimension %dxN.', dim);
            end
        end
    end
end

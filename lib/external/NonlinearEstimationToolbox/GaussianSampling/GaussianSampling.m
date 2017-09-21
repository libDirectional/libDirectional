
classdef GaussianSampling < handle & matlab.mixin.Copyable
    % Abstract base class for multivariate Gaussian sampling methods.
    %
    % GaussianSampling Methods:
    %   copy                - Copy a GaussianSampling instance.
    %   getStdNormalSamples - Get a set of samples approximating a standard normal distribution.
    %   getSamples          - Get a set of samples approximating a Gaussian distribution.
    
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
        function [samples, weights, numSamples] = getSamples(obj, gaussian)
            % Get a set of samples approximating a Gaussian distribution.
            %
            % The number of samples, their positions, and their weights are determined (and
            % controlled) by the respective Gaussian sampling technique.
            %
            % Parameters:
            %   >> gaussian (Gaussian)
            %      Gaussian distribution to be approximated.
            %
            % Returns:
            %   << stateSamples (Matrix)
            %      Column-wise arranged sample positions approximating the Gaussian distribution.
            %
            %   << weights (Row vector)
            %      Either column-wise arranged corresponding sample weights
            %      or single scalar weight in case of equally weighted samples.
            %
            %   << numSamples (Positive scalar)
            %      Number of samples approximating the Gaussian distribution.
            
            if ~Checks.isClass(gaussian, 'Gaussian')
                error('GaussianSampling:InvalidGaussian', ...
                      'gaussian must be a Gaussian distribution.');
            end
            
            [mean, ~, covSqrt] = gaussian.getMeanAndCov();
            dim = gaussian.getDim();
            
            % Get standard normal approximation
            [stdNormalSamples, weights, numSamples] = obj.getStdNormalSamples(dim);
            
            % Generate samples
            samples = covSqrt * stdNormalSamples;
            samples = bsxfun(@plus, samples, mean);
        end
    end
    
    methods (Abstract)
        % Get a set of samples approximating a standard normal distribution.
        %
        % The number of samples, their positions, and their weights are determined (and
        % controlled) by the respective Gaussian sampling technique.
        %
        % Parameters:
        %   >> dimension (Positive scalar)
        %      Dimension of the standard normal distribution.
        %
        % Returns:
        %   << samples (Matrix)
        %      The column-wise arranged sample positions.
        %
        %   << weights (Row vector)
        %      Either column-wise arranged corresponding sample weights
        %      or single scalar weight in case of equally weighted samples.
        %
        %   << numSamples (Positive scalar)
        %      The number of samples.
        [samples, weights, numSamples] = getStdNormalSamples(obj, dimension);
    end
end

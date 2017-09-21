
classdef DiracMixture < Distribution
    % This class represents a multivariate Dirac mixture (sample) distribution.
    %
    % DiracMixture Methods:
    %   DiracMixture     - Class constructor.
    %   copy             - Copy a distribution instance.
    %   set              - Set the parameters of the Dirac mixture distribution.
    %   getDim           - Get the dimension of the distribution.
    %   getMeanAndCov    - Get mean and covariance matrix of the distribution.
    %   drawRndSamples   - Draw random samples from the distribution.
    %   logPdf           - Evaluate the logarithmic probability density function (PDF) of the distribution.
    %   getNumComponents - Get the number of Dirac mixture components.
    %   getComponents    - Get the Dirac mixture components.
    
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
    
    methods (Sealed)
        function obj = DiracMixture(samples, weights)
            % Class constructor.
            %
            % The default constructor results an uninitialized Dirac mixture
            % distribution of zero dimension consisting of zero compoentns.
            %
            % Parameters
            %   >> samples (Matrix)
            %      Column-wise arranged sample positions.
            %
            %   >> weights (Row vector)
            %      Row-wise arranged weights of the samples.
            %      If no weights are passed, all samples are assumed to be equally weighted.
            
            if nargin == 2
                obj.set(samples, weights)
            elseif nargin == 1
                obj.set(samples);
            else
                % Default distribution information
                obj.dim        = 0;
                obj.numComps   = 0;
                obj.samples    = [];
                obj.weights    = [];
                obj.cumWeights = [];
                obj.mean       = [];
                obj.cov        = [];
                obj.covSqrt    = [];
            end
        end
        
        function set(obj, samples, weights)
            % Set the parameters of the Dirac mixture distribution.
            %
            % Parameters
            %   >> samples (Matrix)
            %      Column-wise arranged sample positions.
            %
            %   >> weights (Row vector)
            %      Row-wise arranged weights of the samples.
            %      If no weights are passed, all samples are assumed to be equally weighted.
            
            try
                if ~Checks.isMat(samples)
                    error('DiracMixture:InvalidSamples', ...
                          'samples must be a matrix.');
                end
                
                [obj.dim, obj.numComps] = size(samples);
                
                obj.samples = samples;
                
                if nargin == 3
                    if ~Checks.isNonNegativeRowVec(weights, obj.numComps)
                        error('DiracMixture:InvalidWeights', ...
                              ['weights must be a row vector of dimension' ...
                               '%d containing only non-negative values.'], ...
                              obj.numComps);
                    end
                    
                    % Normalize sample weights
                    sumWeights = sum(weights);
                    
                    if sumWeights <= 0
                        error('DiracMixture:InvalidWeights', ...
                              'Sum of sample weights is not positive.');
                    end
                    
                    obj.weights = weights / sumWeights;
                else
                    % Equally weighted samples
                    obj.weights = repmat(1 / obj.numComps, 1, obj.numComps);
                end
            catch ex
                % Reset all distribution information
                obj.dim        = 0;
                obj.numComps   = 0;
                obj.samples    = [];
                obj.weights    = [];
                obj.cumWeights = [];
                obj.mean       = [];
                obj.cov        = [];
                obj.covSqrt    = [];
                
                ex.rethrow();
            end
        end
        
        function dim = getDim(obj)
            dim = obj.dim;
        end
        
        function [mean, cov, covSqrt] = getMeanAndCov(obj)
            if isempty(obj.mean)
                [obj.mean, obj.cov] = Utils.getMeanAndCov(obj.samples, obj.weights);
            end
            
            mean = obj.mean;
            cov  = obj.cov;
            
            if nargout >= 3
                if isempty(obj.covSqrt)
                    [obj.covSqrt, isNonPos] = chol(obj.cov, 'Lower');
                    
                    if isNonPos
                        error('DiracMixture:InvalidCovariance', ...
                              'Covariance matrix is not positive definite.');
                    end
                end
                
                covSqrt = obj.covSqrt;
            end
        end
        
        function rndSamples = drawRndSamples(obj, numSamples)
            if ~Checks.isPosScalar(numSamples)
                error('DiracMixture:InvalidNumberOfSamples', ...
                      'numSamples must be a positive scalar.');
            end
            
            if isempty(obj.cumWeights)
                obj.cumWeights = cumsum(obj.weights);
            end
            
            rndSamples = Utils.resampling(obj.samples, obj.cumWeights, numSamples);
        end
        
        function logPdf(~, ~)
            error('DiracMixture:LogPdfNotSupported', ...
                  'A Dirac mixture has no proper pdf.');
        end
        
        function numComps = getNumComponents(obj)
            % Get the number of Dirac mixture components.
            %
            % Returns:
            %   << numComps (Scalar )
            %      The number of Dirac mixture components.
            
            numComps = obj.numComps;
        end
        
        function [samples, weights] = getComponents(obj)
            % Get the Dirac mixture components.
            %
            % Returns:
            %   << samples (Matrix)
            %      Column-wise arranged sample positions.
            %
            %   << weights (Row vector)
            %      Row-wise arranged weights of the samples.
            
            samples = obj.samples;
            weights = obj.weights;
        end
    end
    
    properties (Access = 'private')
        % Distribution's dimension.
        dim;
        
        % Number of Dirac mixture components.
        numComps;
        
        % Component positions (samples).
        samples;
        
        % Normalized component weights.
        weights;
        
        % Cumulative sum of component weights.
        cumWeights;
        
        % Mean vector of the Dirac mixture distribution.
        mean;
        
        % Covariance matrix of the Dirac mixture distribution.
        cov;
        
        % Lower Cholesky decomposition of covariance matrix of the Dirac mixture distribution.
        covSqrt;
    end
end


classdef GaussianMixture < Distribution
    % This class represents a multivariate Gaussian mixture distribution.
    %
    % GaussianMixture Methods:
    %   GaussianMixture  - Class constructor.
    %   copy             - Copy a distribution instance.
    %   set              - Set the parameters of the Gaussian mixture distribution.
    %   getDim           - Get the dimension of the distribution.
    %   getMeanAndCov    - Get mean and covariance matrix of the distribution.
    %   drawRndSamples   - Draw random samples from the distribution.
    %   logPdf           - Evaluate the logarithmic probability density function (PDF) of the distribution.
    %   getNumComponents - Get the number of Gaussian mixture components.
    %   getComponents    - Get the Gaussian mixture components.
    
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
        function obj = GaussianMixture(means, covariances, weights)
            % Class constructor.
            %
            % The default constructor results an uninitialized Gaussian mixture
            % distribution of zero dimension consisting of zero compoentns.
            %
            % Parameters
            %   >> means (Matrix)
            %      Column-wise arranged means of the Gaussian mixture components.
            %
            %   >> covariances (3D matrix containing positive definite matrices)
            %      Slice-wise arranged covariance matrices of the Gaussian
            %      mixture components.
            %
            %   >> weights (Row vector)
            %      Row-wise arranged weights of the Gaussian mixture components.
            %      If no weights are passed, all Gaussian mixture components
            %      are assumed to be equally weighted.
            
            if nargin == 3
                obj.set(means, covariances, weights);
            elseif nargin == 2
                obj.set(means, covariances);
            else
                % Default distribution information
                obj.dim          = 0;
                obj.numComps     = 0;
                obj.means        = [];
                obj.covs         = [];
                obj.covSqrts     = [];
                obj.invCovSqrts  = [];
                obj.weights      = [];
                obj.cumWeights   = [];
                obj.logPdfConsts = [];
                obj.mean         = [];
                obj.cov          = [];
                obj.covSqrt      = [];
            end
        end
        
        function set(obj, means, covariances, weights)
            % Set the parameters of the Gaussian mixture distribution.
            %
            % Parameters
            %   >> means (Matrix)
            %      Column-wise arranged means of the Gaussian mixture components.
            %
            %   >> covariances (3D matrix containing positive definite matrices)
            %      Slice-wise arranged covariance matrices of the Gaussian
            %      mixture components.
            %
            %   >> weights (Row vector)
            %      Row-wise arranged weights of the Gaussian mixture components.
            %      If no weights are passed, all Gaussian mixture components
            %      are assumed to be equally weighted.
            
            try
                % Check means
                if ~Checks.isMat(means)
                    error('GaussianMixture:InvalidMeans', ...
                          'means must be a matrix.');
                end
                
                [obj.dim, obj.numComps] = size(means);
                
                obj.means = means;
                
                % Check covariances
                if obj.numComps == 1
                    [isCov, obj.covSqrts] = Checks.isCov(covariances, obj.dim);
                    
                    if ~isCov
                        error('GaussianMixture:InvalidCovariances', ...
                              'covariances must be a positive definite matrix of dimension %dx%d.', ...
                              obj.dim, obj.dim);
                    end
                else
                    [isCov3D, obj.covSqrts] = Checks.isCov3D(covariances, obj.dim, obj.numComps);
                    
                    if ~isCov3D
                        error('GaussianMixture:InvalidCovariances', ...
                              ['covariances must be a matrix of dimension' ...
                               '%dx%dx%d containing positive definite matrices.'], ...
                              obj.dim, obj.dim, obj.numComps);
                    end
                end
                
                obj.covs = covariances;
                
                % Check weights
                if nargin == 4
                    if ~Checks.isNonNegativeRowVec(weights, obj.numComps)
                        error('GaussianMixture:InvalidWeights', ...
                              ['weights must be a row vector of dimension %d' ...
                               'containing only non-negative values.'], ...
                              obj.numComps);
                    end
                    
                    % Normalize component weights
                    sumWeights = sum(weights);
                    
                    if sumWeights <= 0
                        error('GaussianMixture:InvalidWeights', ...
                              'Sum of component weights is not positive.');
                    end
                    
                    obj.weights = weights / sumWeights;
                else
                    % Equally weighted components
                    obj.weights = repmat(1 / obj.numComps, 1, obj.numComps);
                end
            catch ex
                % Reset all distribution information
                obj.dim          = 0;
                obj.numComps     = 0;
                obj.means        = [];
                obj.covs         = [];
                obj.covSqrts     = [];
                obj.invCovSqrts  = [];
                obj.weights      = [];
                obj.cumWeights   = [];
                obj.logPdfConsts = [];
                obj.mean         = [];
                obj.cov          = [];
                obj.covSqrt      = [];
                
                ex.rethrow();
            end
        end
        
        function dim = getDim(obj)
            dim = obj.dim;
        end
        
        function [mean, cov, covSqrt] = getMeanAndCov(obj)
            if isempty(obj.mean)
                [obj.mean, obj.cov] = Utils.getGMMeanAndCov(obj.means, ...
                                                            obj.covs, ...
                                                            obj.weights);
            end
            
            mean = obj.mean;
            cov  = obj.cov;
            
            if nargout >= 3
                if isempty(obj.covSqrt)
                    obj.covSqrt = chol(obj.cov, 'Lower');
                end
                
                covSqrt = obj.covSqrt;
            end
        end
        
        function [rndSamples, compIds] = drawRndSamples(obj, numSamples)
            % Draw random samples from the distribution.
            %
            % Parameters:
            %   >> numSamples (Positive scalar)
            %      Number of samples to draw from the distribution.
            %
            % Returns:
            %   << rndSamples (Matrix)
            %      Column-wise arranged random samples.
            %
            %   << compIds (Row vector)
            %      Column-wise arranged corresponding ID of the Gaussian mixture
            %      component from which a sample was drawn.
            
            if ~Checks.isPosScalar(numSamples)
                error('GaussianMixture:InvalidNumberOfSamples', ...
                      'numSamples must be positive scalar.');
            end
            
            if isempty(obj.cumWeights)
                obj.cumWeights = cumsum(obj.weights);
            end
            
            % Select Gaussian components randomly according to the mixture weights
            u = rand(1, numSamples);
            
            u = sort(u);
            
            numCompsSamples = zeros(1, obj.numComps);
            
            i = 1;
            
            for j = 1:numSamples
                while u(j) > obj.cumWeights(i)
                    i = i + 1;
                end
                
                numCompsSamples(i) = numCompsSamples(i) + 1;
            end
            
            % Generate random samples
            rndSamples = nan(obj.dim, numSamples);
            compIds    = nan(1, numSamples);
            
            a = 1;
            b = 0;
            
            for i = 1:obj.numComps
                numCompSamples = numCompsSamples(i);
                
                if numCompSamples > 0
                    b = b + numCompSamples;
                    
                    % Compute random samples for ith component
                    rndSamples(:, a:b) = Utils.drawGaussianRndSamples(obj.means(:, i), ...
                                                                      obj.covSqrts(:, :, i), ...
                                                                      numCompSamples);
                    
                    % Save corresponding component id
                    compIds(a:b) = i;
                    
                    a = b + 1;
                end
            end
        end
        
        function logValues = logPdf(obj, values)
            obj.checkValues(values);
            
            if isempty(obj.logPdfConsts)
                logNormConst = obj.dim * 0.5 * log(2 * pi);
                
                obj.invCovSqrts  = nan(obj.dim, obj.dim, obj.numComps);
                obj.logPdfConsts = nan(1, obj.numComps);
                
                logWeights = log(obj.weights);
                
                for i = 1:obj.numComps
                    obj.invCovSqrts(:, :, i) = obj.covSqrts(:, :, i) \ eye(obj.dim);
                    
                    logSqrtDetCov = sum(log(diag(obj.covSqrts(:, :, i))));
                    
                    obj.logPdfConsts(i) = logWeights(i) - (logSqrtDetCov + logNormConst);
                end
            end
            
            compValues = nan(obj.numComps, size(values, 2));
            
            for i = 1:obj.numComps
                s = bsxfun(@minus, values, obj.means(:, i));
                
                v = obj.invCovSqrts(:, :, i) * s;
                
                compValues(i, :) = obj.logPdfConsts(i) - 0.5 * sum(v.^2, 1);
            end
            
            maxLogCompValues = max(compValues);
            
            compValues = bsxfun(@minus, compValues, maxLogCompValues);
            
            compValues = exp(compValues);
            
            logValues = maxLogCompValues + log(sum(compValues, 1));
        end
        
        function numComps = getNumComponents(obj)
            % Get the number of Gaussian mixture components.
            %
            % Returns:
            %   << numComps (Non-negative scalar)
            %      The number of Gaussian mixture components.
            
            numComps = obj.numComps;
        end
        
        function [means, covariances, weights] = getComponents(obj)
            % Get the Gaussian mixture components.
            %
            % Returns:
            %   << means (Matrix)
            %      Column-wise arranged means of the Gaussian mixture components.
            %
            %   << covariances (3D matrix containing positive definite matrices)
            %      Slice-wise arranged covariance matrices of the Gaussian mixture components.
            %
            %   << weights (Row vector)
            %      Row-wise arranged weights of the Gaussian mixture components.
            
            means       = obj.means;
            covariances = obj.covs;
            weights     = obj.weights;
        end
    end
    
    properties (Access = 'private')
        % Distribution's dimension.
        dim;
        
        % Number of Gaussian mixture components.
        numComps;
        
        % Component mean vectors.
        means;
        
        % Component covariance matrices.
        covs;
        
        % Lower Cholesky decompositions of component covariance matrices.
        covSqrts;
        
        % Inverse lower Cholesky decompositions of component covariance matrices.
        invCovSqrts;
        
        % Normalized component weights.
        weights;
        
        % Cumulative sum of component weights.
        cumWeights;
        
        % Logarithms of the Gaussian component PDF's normalization constants.
        logPdfConsts;
        
        % Mean vector of the Gaussian mixture distribution.
        mean;
        
        % Covariance matrix of the Gaussian mixture distribution.
        cov;
        
        % Lower Cholesky decomposition of covariance matrix of the Gaussian mixture distribution.
        covSqrt;
    end
end

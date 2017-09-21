
classdef Gaussian < Distribution
    % This class represents a multivariate Gaussian distribution.
    %
    % Gaussian Methods:
    %   Gaussian       - Class constructor.
    %   copy           - Copy a distribution instance.
    %   set            - Set the parameters of the Gaussian distribution.
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
    
    methods (Sealed)
        function obj = Gaussian(mean, covariance)
            % Class constructor.
            %
            % The default constructor results an uninitialized Gaussian
            % distribution of zero dimension.
            %
            % Parameters:
            %   >> mean (Column vector)
            %      Mean vector of the distribution.
            %
            %   >> covariance (Positive definite matrix, vector, or 3D matrix of positive definite matrices)
            %      Covariance matrix of the distribution. If a vector is
            %      passed, its values are interpreted as the variances of
            %      a diagonal covariance matrix. If a 3D matrix is passed,
            %      a block diagonal covariance matrix will be constructed.
            
            if nargin == 2
                obj.set(mean, covariance);
            else
                % Default distribution information
                obj.dim         = 0;
                obj.mean        = [];
                obj.cov         = [];
                obj.covSqrt     = [];
                obj.invCovSqrt  = [];
                obj.logPdfConst = [];
            end
        end
        
        function set(obj, mean, covariance)
            % Set the parameters of the Gaussian distribution.
            %
            % Parameters:
            %   >> mean (Column vector)
            %      Mean vector of the distribution.
            %
            %   >> covariance (Positive definite matrix, vector, or 3D matrix of positive definite matrices)
            %      Covariance matrix of the distribution. If a vector is
            %      passed, its values are interpreted as the variances of
            %      a diagonal covariance matrix. If a 3D matrix is passed,
            %      a block diagonal covariance matrix will be constructed.
            
            try
                if ~Checks.isColVec(mean)
                    error('Gaussian:InvalidMean', ...
                          'mean must be a column vector.');
                end
                
                obj.dim         = size(mean, 1);
                obj.mean        = mean;
                obj.logPdfConst = [];
                obj.invCovSqrt  = [];
                
                if Checks.isPosVec(covariance, obj.dim);
                    obj.cov     = diag(covariance);
                    obj.covSqrt = diag(sqrt(covariance));
                else
                    [isCov, obj.covSqrt] = Checks.isCov(covariance, obj.dim);
                    
                    if isCov
                        obj.cov = covariance;
                    else
                        [isCov3D, covSqrts] = Checks.isCov3D(covariance);
                        
                        if isCov3D
                            [dimCov, ~, numCovs] = size(covariance);
                            
                            if (dimCov * numCovs) ~= obj.dim
                                error('Gaussian:InvalidCovariance', ...
                                      'Constructed block covariance matrix has to be of dimension %dx%d.', ...
                                      obj.dim, obj.dim);
                            end
                            
                            covs     = mat2cell(covariance, dimCov, dimCov, ones(1, numCovs));
                            covSqrts = mat2cell(covSqrts, dimCov, dimCov, ones(1, numCovs));
                            
                            obj.cov     = blkdiag(covs{:});
                            obj.covSqrt = blkdiag(covSqrts{:});
                        else
                            error('Gaussian:InvalidCovariance', ...
                                  ['covariance must be\n' ...
                                   '  * a vector of dimension %d containing positive values only, or\n' ...
                                   '  * a positive definite matrix of dimension %dx%d, or\n' ...
                                   '  * a 3D matrix of positive definite matrices.'], ...
                                  obj.dim, obj.dim, obj.dim);
                        end
                    end
                end
            catch ex
                % Reset all distribution information
                obj.dim         = 0;
                obj.mean        = [];
                obj.cov         = [];
                obj.covSqrt     = [];
                obj.invCovSqrt  = [];
                obj.logPdfConst = [];
                
                ex.rethrow();
            end
        end
        
        function dim = getDim(obj)
            dim = obj.dim;
        end
        
        function [mean, cov, covSqrt] = getMeanAndCov(obj)
            mean    = obj.mean;
            cov     = obj.cov;
            covSqrt = obj.covSqrt;
        end
        
        function rndSamples = drawRndSamples(obj, numSamples)
            if ~Checks.isPosScalar(numSamples)
                error('Gaussian:InvalidNumberOfSamples', ...
                      'numSamples must be a positive scalar.');
            end
            
            % Generate random samples
            rndSamples = Utils.drawGaussianRndSamples(obj.mean, obj.covSqrt, numSamples);
        end
        
        function logValues = logPdf(obj, values)
            obj.checkValues(values);
            
            if isempty(obj.logPdfConst)
                obj.invCovSqrt = obj.covSqrt \ eye(obj.dim);
                
                logNormConst = obj.dim * 0.5 * log(2 * pi);
                
                logSqrtDetCov = sum(log(diag(obj.covSqrt)));
                
                obj.logPdfConst = logSqrtDetCov + logNormConst;
            end
            
            s = bsxfun(@minus, values, obj.mean);
            
            v = obj.invCovSqrt * s;
            
            logValues = -0.5 * sum(v.^2, 1) - obj.logPdfConst;
        end
    end
    
    properties (Access = 'private')
        % Distribution's dimension.
        dim;
        
        % Mean vector.
        mean;
        
        % Covariance matrix.
        cov;
        
        % Lower Cholesky decomposition of the covariance matrix.
        covSqrt;
        
        % Inverse of the lower Cholesky decomposition of the covariance matrix.
        invCovSqrt;
        
        % Logarithm of the Gaussian PDF's normalization constant.
        logPdfConst;
    end
end

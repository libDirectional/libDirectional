
classdef Gaussian < Distribution
    % This class represents a multivariate Gaussian distribution.
    % 
    % Gaussian Methods:
    %   Gaussian             - Class constructor.
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
    %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    methods
        function obj = Gaussian(mean, covariance)
            % Class constructor.
            %
            % Parameters:
            %   >> mean (Column vector)
            %      Mean vector of the distribution.
            %      Default: 0.
            %
            %   >> covariance (Positive definite matrix or vector)
            %      Covariance matrix of the distribution. If a vector is
            %      passed, its values are interpreted as the variances of
            %      a diagonal covariance matrix.
            %      Default: 1.
            
            if nargin == 2
                obj.set(mean, covariance);
            else
                obj.set(0, 1);
            end
            
            obj.logPdfConst = [];
            obj.invCovSqrt  = [];
        end
        
        function dimension = getDimension(obj)
            dimension = obj.dimension;
        end
        
        function [mean, covariance, covSqrt] = getMeanAndCovariance(obj)
            mean       = obj.mean;
            covariance = obj.covariance;
            
            if nargout >= 3
                covSqrt = obj.covSqrt;
            end
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
                obj.invCovSqrt = obj.covSqrt \ eye(obj.dimension);
                
                logNormConst = obj.dimension * 0.5 * log(2 * pi);
                
                logSqrtDetCov = sum(log(diag(obj.covSqrt)));
                
                obj.logPdfConst = logSqrtDetCov + logNormConst;
            end
            
            s = bsxfun(@minus, values, obj.mean);
            
            v = obj.invCovSqrt * s;
            
            logValues = -0.5 * sum(v.^2, 1) - obj.logPdfConst;
        end
    end
    
    methods (Access = 'private')
        function set(obj, mean, covariance)
            if ~Checks.isColVec(mean)
                error('Gaussian:InvalidMean', ...
                      'mean must be a column vector.');
            end
            
            dim = size(mean, 1);
            
            obj.dimension = dim;
            obj.mean      = mean;
            
            if Checks.isPosVec(covariance, dim);
                obj.covariance = diag(covariance);
                obj.covSqrt    = diag(sqrt(covariance));
            else
                [isCov, obj.covSqrt] = Checks.isCov(covariance, dim);
                
                if isCov
                    obj.covariance = covariance;
                else
                    [isCov3D, covSqrts] = Checks.isCov3D(covariance);
                    
                    if isCov3D
                        [dimCov, ~, numCovs] = size(covariance);
                        
                        if dim ~= (dimCov * numCovs)
                            error('Gaussian:InvalidCovariance', ...
                                  'Constructed block covariance matrix has to be of dimension %dx%d.', ...
                                  dim, dim);
                        end
                        
                        covs     = mat2cell(covariance, dimCov, dimCov, ones(1, numCovs));
                        covSqrts = mat2cell(covSqrts, dimCov, dimCov, ones(1, numCovs));
                        
                        obj.covariance = blkdiag(covs{:});
                        obj.covSqrt    = blkdiag(covSqrts{:});
                    else
                        error('Gaussian:InvalidCovariance', ...
                              ['covariance must be\n' ...
                               '  * a vector of dimension %d containing positive values only, or\n' ...
                               '  * a positive definite matrix of dimension %dx%d, or\n' ...
                               '  * a 3D matrix with positive definite submatrices.'], ...
                              dim, dim, dim);
                    end
                end
            end
        end
    end
    
    properties (Access = 'private')
        dimension;
        mean;
        covariance;
        covSqrt;
        
        logPdfConst;
        invCovSqrt;
    end
end

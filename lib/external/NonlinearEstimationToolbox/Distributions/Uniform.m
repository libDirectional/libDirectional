
classdef Uniform < Distribution
    % This class represents a multivariate axis-aligned uniform distribution.
    %
    % Uniform Methods:
    %   Uniform        - Class constructor.
    %   copy           - Copy a distribution instance.
    %   set            - Set the parameters of the uniform distribution.
    %   getDim         - Get the dimension of the distribution.
    %   getMeanAndCov  - Get mean and covariance matrix of the distribution.
    %   drawRndSamples - Draw random samples from the distribution.
    %   logPdf         - Evaluate the logarithmic probability density function (PDF) of the distribution.
    %   getInterval    - Get the support of the uniform distribution.
    
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
        function obj = Uniform(a, b)
            % Class constructor
            %
            % The default constructor results an uninitialized uniform
            % distribution of zero dimension.
            %
            % Parameters
            %   >> a (Vector)
            %      Lower bounds of the multivariate uniform distribution.
            %
            %   >> b (Vector)
            %      Upper bounds of the multivariate uniform distribution.
            %      Must have the same length as a.
            
            if nargin == 2
                obj.set(a, b);
            else
                % Default distribution information
                obj.dim           = 0;
                obj.a             = [];
                obj.b             = [];
                obj.lengths       = [];
                obj.validLogValue = [];
                obj.mean          = [];
                obj.cov           = [];
                obj.covSqrt       = [];
            end
        end
        
        function set(obj, a, b)
            % Set the parameters of the uniform distribution.
            %
            % Parameters
            %   >> a (Vector)
            %      Lower bounds of the multivariate uniform distribution.
            %
            %   >> b (Vector)
            %      Upper bounds of the multivariate uniform distribution.
            %      Must have the same length as a.
            
            try
                if ~Checks.isVec(a)
                    error('Uniform:InvalidMinimum', ...
                          'a must be vector.');
                end
                
                obj.a = a(:);
                
                obj.dim = size(obj.a, 1);
                
                if ~Checks.isVec(b, obj.dim)
                    error('Uniform:InvalidMaximum', ...
                          'b must be vector of length %d.', obj.dim);
                end
                
                if any(a >= b)
                    error('Uniform:InvalidMinimum', ...
                          'All entries of a must be smaller than their corresponding entries in b.');
                end
                
                obj.b = b(:);
                
                obj.lengths       = obj.b - obj.a;
                obj.validLogValue = -log(prod(obj.lengths));
            catch ex
                % Reset all distribution information
                obj.dim           = 0;
                obj.a             = [];
                obj.b             = [];
                obj.lengths       = [];
                obj.validLogValue = [];
                obj.mean          = [];
                obj.cov           = [];
                obj.covSqrt       = [];
                
                ex.rethrow();
            end
        end
        
        function dim = getDim(obj)
            dim = obj.dim;
        end
        
        function [mean, cov, covSqrt] = getMeanAndCov(obj)
            if isempty(obj.mean)
                obj.mean = 0.5 * (obj.a + obj.b);
                obj.cov  = diag((obj.b - obj.a).^2 / 12);
            end
            
            mean = obj.mean;
            cov  = obj.cov;
            
            if nargout >= 3
                if isempty(obj.covSqrt)
                    obj.covSqrt = sqrt(obj.cov);
                end
                
                covSqrt = obj.covSqrt;
            end
        end
        
        function rndSamples = drawRndSamples(obj, numSamples)
            if ~Checks.isPosScalar(numSamples)
                error('Uniform:InvalidNumberOfSamples', ...
                      'numSamples must be positive scalar.');
            end
            
            rndSamples = bsxfun(@times, obj.lengths, rand(obj.dim, numSamples));
            rndSamples = bsxfun(@plus, obj.a, rndSamples);
        end
        
        function logValues = logPdf(obj, values)
            obj.checkValues(values);
            
            numValues = size(values, 2);
            
            idxValid                   = zeros(2 * obj.dim, numValues);
            idxValid(1:obj.dim, :)     = bsxfun(@ge, values, obj.a);
            idxValid(obj.dim+1:end, :) = bsxfun(@le, values, obj.b);
            
            idx = all(idxValid, 1);
            
            logValues      = -inf(1, numValues);
            logValues(idx) = obj.validLogValue;
        end
        
        function [a, b] = getInterval(obj)
            % Get the support of the uniform distribution.
            %
            % Returns:
            %   << a (Vector)
            %      Lower bounds of the multivariate uniform distribution.
            %
            %   << b (Vector)
            %      Upper bounds of the multivariate uniform distribution.
            
            a = obj.a;
            b = obj.b;
        end
    end
    
    properties (Access = 'private')
        % Distribution's dimension.
        dim;
        
        % Lower bounds.
        a;
        
        % Upper bounds.
        b;
        
        % Lengths of the respective intervals.
        lengths;
        
        % Logarithm of the uniform's valid PDF values.
        validLogValue;
        
        % Mean vector of the uniform distribution.
        mean;
        
        % Covariance matrix of the uniform distribution.
        cov;
        
        % Lower Cholesky decomposition of covariance matrix of the uniform distribution.
        covSqrt;
    end
end

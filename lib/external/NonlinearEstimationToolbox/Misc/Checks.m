
classdef Checks
    % This class provides various runtime checks.
    %
    % Checks Methods:
    %   Scalar checks
    %     isScalar            - Checks for a scalar.
    %     isNonNegativeScalar - Checks for a non-negative scalar.
    %     isPosScalar         - Checks for a positive scalar.
    %     isScalarIn          - Checks for a scalar in a specified interval.
    %
    %   Vector checks
    %     isVec               - Checks for a column or row vector.
    %     isNonNegativeVec    - Checks for a column or row vector with non-negative entries.
    %     isPosVec            - Checks for a column or row vector with positive entries.
    %     isColVec            - Checks for a column vector.
    %     isNonNegativeColVec - Checks for a column vector with non-negative entries.
    %     isPosColVec         - Checks for a column vector with positive entries.
    %     isRowVec            - Checks for a row vector.
    %     isNonNegativeRowVec - Checks for a row vector with non-negative entries.
    %     isPosRowVec         - Checks for a row vector with positive entries.
    %
    %   Matrix checks
    %     isMat               - Checks for a numeric and non-empty matrix.
    %     isFixedRowMat       - Checks for a numeric matrix with a given number of rows.
    %     isFixedColMat       - Checks for a numeric matrix with a given number of columns.
    %     isSquareMat         - Checks for a numeric and square matrix.
    %
    %   Covariance checks
    %     isCov               - Checks for a numeric, square, and positive definite matrix.
    %     isCov3D             - Checks for numeric, square, and positive definite matrices stored along the third dimension.
    %
    %   Misc checks
    %     isFlag              - Checks for a boolean scalar.
    %     isClass             - Checks for a given class type.
    
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
    
    methods (Static)
        function ret = isScalar(scalar)
            ret = isnumeric(scalar) && ...
                  isscalar(scalar);
        end
        
        function ret = isNonNegativeScalar(scalar)
            ret = isnumeric(scalar) && ...
                  isscalar(scalar) && ...
                  scalar >= 0;
        end
        
        function ret = isPosScalar(scalar)
            ret = isnumeric(scalar) && ...
                  isscalar(scalar) && ...
                  scalar > 0;
        end
        
        function ret = isScalarIn(scalar, a, b)
            ret = isnumeric(scalar) && ...
                  isscalar(scalar) && ...
                  scalar >= a && ...
                  scalar <= b;
        end
        
        
        function ret = isVec(vec, dim)
            if nargin == 1
                ret = isnumeric(vec) && ...
                      isvector(vec) && ...
                      ~isempty(vec);
            else
                ret = isnumeric(vec) && ...
                      isvector(vec) && ...
                      length(vec) == dim;
            end
        end
        
        function ret = isNonNegativeVec(vec, dim)
            if nargin == 1
                ret = Checks.isVec(vec) && ...
                      all(vec >= 0);
            else
                ret = Checks.isVec(vec, dim) && ...
                      all(vec >= 0);
            end
        end
        
        function ret = isPosVec(vec, dim)
            if nargin == 1
                ret = Checks.isVec(vec) && ...
                      all(vec > 0);
            else
                ret = Checks.isVec(vec, dim) && ...
                      all(vec > 0);
            end
        end
        
        function ret = isColVec(vec, dim)
            if nargin == 1
                ret = isnumeric(vec) && ...
                      iscolumn(vec) && ...
                      ~isempty(vec);
            else
                ret = isnumeric(vec) && ...
                      iscolumn(vec) && ...
                      length(vec) == dim;
            end
        end
        
        function ret = isNonNegativeColVec(vec, dim)
            if nargin == 1
                ret = Checks.isColVec(vec) && ...
                      all(vec >= 0);
            else
                ret = Checks.isColVec(vec, dim) && ...
                      all(vec >= 0);
            end
        end
        
        function ret = isPosColVec(vec, dim)
            if nargin == 1
                ret = Checks.isColVec(vec) && ...
                      all(vec > 0);
            else
                ret = Checks.isColVec(vec, dim) && ...
                      all(vec > 0);
            end
        end
        
        function ret = isRowVec(vec, dim)
            if nargin == 1
                ret = isnumeric(vec) && ...
                      isrow(vec) && ...
                      ~isempty(vec);
            else
                ret = isnumeric(vec) && ...
                      isrow(vec) && ...
                      length(vec) == dim;
            end
        end
        
        function ret = isNonNegativeRowVec(vec, dim)
            if nargin == 1
                ret = Checks.isRowVec(vec) && ...
                      all(vec >= 0);
            else
                ret = Checks.isRowVec(vec, dim) && ...
                      all(vec >= 0);
            end
        end
        
        function ret = isPosRowVec(vec, dim)
            if nargin == 1
                ret = Checks.isRowVec(vec) && ...
                      all(vec > 0);
            else
                ret = Checks.isRowVec(vec, dim) && ...
                      all(vec > 0);
            end
        end
        
        
        function ret = isMat(mat, rows, cols)
            if nargin == 1
                ret = isnumeric(mat) && ...
                      ismatrix(mat) && ...
                      ~isempty(mat);
            else
                ret = isnumeric(mat) && ...
                      ismatrix(mat) && ...
                      size(mat, 1) == rows && ...
                      size(mat, 2) == cols;
            end
        end
        
        function ret = isFixedRowMat(mat, rows)
                ret = isnumeric(mat) && ...
                      ismatrix(mat) && ...
                      size(mat, 1) == rows;
        end
        
        function ret = isFixedColMat(mat, cols)
                ret = isnumeric(mat) && ...
                      ismatrix(mat) && ...
                      size(mat, 2) == cols;
        end
        
        function ret = isSquareMat(mat, dim)
            if nargin == 1
                ret = isnumeric(mat) && ...
                      ismatrix(mat) && ...
                      ~isempty(mat) && ...
                      size(mat, 1) == size(mat, 2);
            else
                ret = isnumeric(mat) && ...
                      ismatrix(mat) && ...
                      size(mat, 1) == dim && ...
                      size(mat, 2) == dim;
            end
        end
        
        function ret = isMat3D(mat, rows, cols, slices)
            nDims = numel(size(mat));
            
            if nargin == 1
                ret = isnumeric(mat) && ...
                      (nDims == 2 || nDims == 3) && ...
                      ~isempty(mat);
            else
                ret = isnumeric(mat) && ...
                      (nDims == 2 || nDims == 3) && ...
                      size(mat, 1)     == rows && ...
                      size(mat, 2)     == cols && ...
                      size(mat, 3)     == slices;
            end
        end
        
        
        function [ret, covSqrt] = isCov(cov, dim)
            if nargin == 1
                if Checks.isSquareMat(cov)
                    [covSqrt, isNonPos] = chol(cov, 'Lower');
                    
                    if isNonPos
                        ret = false;
                        covSqrt = [];
                        return;
                    end
                else
                    ret = false;
                    covSqrt = [];
                    return;
                end
            else
                if Checks.isSquareMat(cov, dim)
                    [covSqrt, isNonPos] = chol(cov, 'Lower');
                    
                    if isNonPos
                        ret = false;
                        covSqrt = [];
                        return;
                    end
                else
                    ret = false;
                    covSqrt = [];
                    return;
                end
            end
            
            ret = true;
        end
        
        function [ret, covSqrts] = isCov3D(cov, dim, numCovs)
            if ~isnumeric(cov)
                ret = false;
                covSqrts = [];
                return;
            end
            
            [dimA, dimB, N] = size(cov);
            
            if nargin < 2
                if dimA ~= dimB || dimA == 0
                    ret = false;
                    covSqrts = [];
                    return;
                end
            else
                if ~(dimA == dim && ...
                     dimB == dim)
                    ret = false;
                    covSqrts = [];
                    return;
                end
            end
            
            if nargin == 3
                if N ~= numCovs
                    ret = false;
                    covSqrts = [];
                    return;
                end
            end
            
            covSqrts = nan(dimA, dimB, N);
            
            for i = 1:N
                [covSqrt, isNonPos] = chol(cov(:, :, i), 'Lower');
                
                if isNonPos
                    ret = false;
                    covSqrts = [];
                    return;
                else
                    covSqrts(:, :, i) = covSqrt;
                end
            end
            
            ret = true;
        end
        
        
        function ret = isFlag(flag)
            ret = islogical(flag) && ...
                  isscalar(flag);
        end
        
        function ret = isClass(obj, className, numObjs)
            if nargin == 2
                ret = isa(obj, className) && ...
                      numel(obj) == 1;
            else
                ret = isa(obj, className) && ...
                      numel(obj) == numObjs;
            end
        end
    end
end

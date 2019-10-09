classdef FIGDistribution < AbstractCircularDistribution
    properties
        gridValues
        transformation
        enforcePdfNonnegative
    end
    methods
        function this = FIGDistribution(gridValues, transformation, enforcePdfNonnegative)
            if nargin<3, enforcePdfNonnegative=true; end
            if isa(gridValues, 'AbstractCircularDistribution')
                error('You gave a distribution as the first argument. To convert distributions to a distribution in grid representation, use .fromDistribution');
            end
            assert(size(gridValues, 2) == 1, 'Provided grid values need to be given as a column vector.');
            assert(isreal(gridValues) && all(gridValues >= 0), ...
                'Grid values must be positive real numbers.');
            this.gridValues = gridValues;
            this.transformation = transformation;
            this.enforcePdfNonnegative = enforcePdfNonnegative;
            % Check if normalized. If not: Normalize!
            this = this.normalize;
        end
        
        function p = pdf(this, xa)
            % Value when interpolated using sinc
            noCoeffs = numel(this.gridValues);
            if mod(numel(this.gridValues), 2) == 0
                warnStruct = warning('off', 'Truncate:TooFewCoefficients');
                noCoeffs = noCoeffs + 1;
            end
            if ~this.enforcePdfNonnegative && strcmp(this.transformation, 'identity')
                % Use square root to enforce negative values
                fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'identity');
            elseif this.enforcePdfNonnegative && strcmp(this.transformation, 'identity')
                fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'sqrt');
            elseif strcmp(this.transformation, 'sqrt')
                fd = FourierDistribution.fromFunctionValues((this.gridValues').^2, noCoeffs, this.transformation);
            else
                error('Transformation not supported');
            end
            if mod(numel(this.gridValues), 2) == 0
                warning(warnStruct);
            end
            p = fd.pdf(xa);
        end
        
        function [vals, xGrid] = pdfOnGrid(this, noOfGridpoints)
            xGrid = 0:2 * pi / numel(this.gridValues):2 * pi - 2 * pi / numel(this.gridValues);
            step = numel(this.gridValues) / noOfGridpoints;
            assert(mod(step, 1) == 0, ...
                'Number of function values has to be a multiple of noOfGridpoints');
            
            vals=this.gridValues(1:step:end);
            switch this.transformation
                case 'identity'
                    vals = vals';
                case 'sqrt'
                    vals = vals.^2';
                case 'square'
                    vals = sqrt(vals)';
                otherwise
                    error('PdfOnGrid:unsupportedTransformation', 'Transformation not recognized or unsupported');
            end
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            noCoeffs = numel(this.gridValues);
            switch this.transformation
                case 'identity'
                    fvals = this.gridValues;
                case 'sqrt'
                    fvals = this.gridValues.^2;
                otherwise
                    error('Moment:unsupportedTransformation', 'Transformation not recognized or unsupported');
            end
            if mod(numel(this.gridValues), 2) == 0
                warnStruct = warning('off', 'Truncate:TooFewCoefficients');
                noCoeffs = noCoeffs + 1;
            end
            fd = FourierDistribution.fromFunctionValues(fvals', noCoeffs, 'identity');
            if mod(numel(this.gridValues), 2) == 0
                warning(warnStruct);
            end
            m = fd.trigonometricMoment(n);
        end
        
        function p = plot(this, varargin)
            gridPoints = 0:2 * pi / numel(this.gridValues):2 * pi - 2 * pi / numel(this.gridValues);
            holdStatus = ishold;
            p(1) = plot@AbstractCircularDistribution(this, varargin{:});
            hold on
            switch this.transformation
                case 'identity'
                    p(2) = plot(gridPoints, this.gridValues, 'x');
                case 'sqrt'
                    p(2) = plot(gridPoints, this.gridValues.^2, 'x');
                otherwise
                    error('Transformation unsupported');
            end
            if ~holdStatus
                hold off
            end
        end
        
        function p = value(this, xa)
            % Value when interpreted as a pmf
            gridPoints = (0:2 * pi / numel(this.gridValues):2 * pi - 2 * pi / numel(this.gridValues))';
            p = (repmat(xa, 1, 5) == gridPoints)' * this.gridValues;
        end
        
        function f = multiply(this, f2)
            % Multiplies two transformed fourier pdfs (returns transformed result)
            if ~strcmp(this.transformation, f2.transformation)
                error('Multiply:differentTransformations', 'Transformations do not match, transform before using multiply');
            end
            if ~(strcmp(this.transformation, 'identity') || strcmp(this.transformation, 'sqrt'))
                error('Multiply:unsupportedTransformation', 'Transformation not recognized or unsupported');
            end
            assert(this.enforcePdfNonnegative == f2.enforcePdfNonnegative);
            assert(numel(this.gridValues) == numel(f2.gridValues));
            warnStruct = warning('off', 'Normalization:notNormalized');
            f = FIGDistribution(this.gridValues.*f2.gridValues, this.transformation, this.enforcePdfNonnegative);
            warning(warnStruct);
        end
        
        function f = convolve(this, f2)
            % Multiplies two transformed fourier pdfs (returns transformed result)
            if ~strcmp(this.transformation, f2.transformation)
                error('Convolve:differentTransformations', 'Transformations do not match, transform before using multiply');
            end
            if ~(strcmp(this.transformation, 'identity') || strcmp(this.transformation, 'sqrt'))
                error('Convolve:unsupportedTransformation', 'Transformation not recognized or unsupported');
            end
            assert(this.enforcePdfNonnegative == f2.enforcePdfNonnegative);
            assert(numel(this.gridValues) == numel(f2.gridValues))
            switch this.transformation
                case 'identity'
                    convResult = cconv(this.gridValues, f2.gridValues, numel(f2.gridValues)) * 2 * pi / numel(f2.gridValues);
                    % Due to numerical inaccuracies, negative values can
                    % occur
                    convResult(convResult < 0) = 0;
                    f = FIGDistribution(convResult, this.transformation, this.enforcePdfNonnegative);
                case 'sqrt'
                    convResult = cconv(this.gridValues.^2, f2.gridValues.^2, numel(f2.gridValues)) * 2 * pi / numel(f2.gridValues);
                    % Due to numerical inaccuracies, negative values can
                    % occur
                    convResult(convResult < 0) = 0;
                    f = FIGDistribution(sqrt(convResult), this.transformation, this.enforcePdfNonnegative);
                otherwise
                    error('Transformation not supported');
            end
        end
        
        function f = truncate(this, noOfCoefficients)
            % Truncates Fourier series. Fills up if there are less coefficients
            % Expects number of complex coefficients (or sum of number of real
            % coefficients)
            assert(noOfCoefficients - 1 > 0, 'Number of coefficients must be an integer greater zero');
            step = numel(this.gridValues) / noOfGridpoints;
            if mod(step, 1) == 0
                f = FIGDistribution(this.gridValues(1:step:end),this.transformation);
            elseif numel(this.gridPoints) < noOfCoefficients
                warning('Truncate:TooFewGridPoints', 'Less coefficients than desired, interpolating using Fourier sqrt');
                error('To be implemented.');
            else
                warning('Truncate:DownsampleViaFourier', 'Transforming to Fourier to interpolate');
                error('To be implemented.');
            end
            % Truncation can void normalization if transformation is not
            % identity
            warnStruct = warning('off', 'Normalization:notNormalized');
            f = f.normalize;
            warning(warnStruct);
        end
        
        function f = normalize(this)
            % This is just an approximate normalization so that all values
            % sum up to 1/(2*pi), it does not necessarily normalize the
            % continuous density interpolated using the sinc
            tol = 1e-4;
            switch this.transformation
                case 'identity'
                    integral = 2 * pi * mean(this.gridValues);
                    normalizationFactor = integral;
                case 'sqrt'
                    integral = 2 * pi * sum(this.gridValues.^2/numel(this.gridValues));
                    normalizationFactor = sqrt(integral);
                otherwise
                    warning('Normalization:cannotTest', 'Unable to test if normalized');
                    f = this;
                    return
            end
            
            f = this;
            if any(this.gridValues < 0)
                warning('Normalization:negative', 'There are negative values. This usualy points to a user error.');
            elseif abs(integral) < 1e-200 % Tolerance has to be that low to avoid unnecessary errors on multiply
                error('Normalization:almostZero', 'Sum of grid vallues is too close to zero, this usually points to a user error.');
            elseif abs(integral-1) > tol
                warning('Normalization:notNormalized', 'Grid values apparently do not belong to normalized density. Normalizing...');
            else
                return % Normalized, return original density
            end
            f = FIGDistribution(f.gridValues/normalizationFactor, this.transformation);
        end
        
        function f = transformViaValues(this, desiredTransformation)
            % Calculates transformations using Fourier coefficients
            
            switch desiredTransformation
                case 'identity'
                    f = this;
                case 'square'
                    switch this.transformation
                        case 'sqrt'
                            newTrans = 'identity';
                        case 'identity'
                            newTrans = 'square';
                        otherwise
                            newTrans = 'multiple';
                    end
                    f = FIGDistribution(this.gridValues.^2, newTrans);
                case 'sqrt'
                    switch this.transformation
                        case 'square'
                            newTrans = 'identity';
                        case 'identity'
                            newTrans = 'sqrt';
                        otherwise
                            newTrans = 'multiple';
                    end
                    f = FIGDistribution(sqrt(this.gridValues), newTrans);
                otherwise
                    error('Desired transformation not supported via coefficients');
            end
        end
        
        function gd = shift(this, angle)
            % Shift distribution by the given angle
            %
            % Parameters:
            %   shiftAngles (scalar)
            %       angle to shift by
            % Returns:
            %   gd (FIGDistribution)
            %       shifted distribution
            assert(isscalar(angle));
            fd = FourierDistribution.fromFunctionValues(this.gridValues', numel(this.gridValues), this.transformation);
            fd = fd.shift(angle);
            gd = FIGDistribution(fd.pdfOnGrid(numel(this.gridValues))', this.transformation);
        end
    end
    
    methods(Static)
        function f = fromDistribution(distribution, noOfGridpoints, desiredTransformation)
            f = FIGDistribution.fromFunction(@(x)distribution.pdf(x), noOfGridpoints, desiredTransformation);
        end
        
        function f = fromFunction(fun, noOfCoefficients, desiredTransformation)
            gridPoints = 0:2 * pi / noOfCoefficients:2 * pi - 2 * pi / noOfCoefficients;
            gridValues = fun(gridPoints)';
            if strcmp(desiredTransformation, 'sqrt')
                gridValues = sqrt(gridValues);
            end
            f = FIGDistribution(gridValues, desiredTransformation);
        end
        
        function gd = fromFunctionValues(fvals, noOfGridpoints, desiredTransformation)
            % Mostly useful for compatibility with FourierDistribution
            step = numel(fvals) / noOfGridpoints;
            assert(mod(step, 1) == 0, ...
                'Number of function values has to be a multiple of noOfGridpoints');
            fvals = fvals(1:step:end);
            if strcmp(desiredTransformation, 'sqrt')
                fvals = sqrt(fvals);
            end
            gd = FIGDistribution(fvals, desiredTransformation);
        end
        
        
    end
    
    
end
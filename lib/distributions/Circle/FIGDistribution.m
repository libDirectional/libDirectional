classdef FIGDistribution < AbstractCircularDistribution
    % This density representation is used by the Fourier-interpreted grid
    % filter. The density is represented using function values on a grid.
    % The interpolation is performed using Fourier series. This class is
    % thus related to FourierDistribution
    
    % see Florian Pfaff, Kailai Li, and Uwe D. Hanebeck,
    % Fourier Filters, Grid Filters, and the Fourier-Interpreted Grid Filter,
    % Proceedings of the 22nd International Conference on Information Fusion (Fusion 2019), Ottawa, Canada, July, 2019.
    properties
        gridValues (:,1) double {mustBeReal, mustBeNonnegative}
        enforcePdfNonnegative logical
    end
    methods
        function this = FIGDistribution(gridValues, enforcePdfNonnegative)
            arguments
                gridValues
                enforcePdfNonnegative logical = true
            end
            if isa(gridValues, 'AbstractCircularDistribution')
                error('You gave a distribution as the first argument. To convert distributions to a distribution in grid representation, use .fromDistribution');
            end
            this.gridValues = gridValues;
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
            if this.enforcePdfNonnegative
               % Use square root to enforce nonnegative values
               fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'sqrt');
            else
                fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'identity');
            end
            if mod(numel(this.gridValues), 2) == 0
                warning(warnStruct);
            end
            p = fd.pdf(xa);
        end
        
        function [vals, xGrid] = pdfOnGrid(this, noOfDesiredGridpoints)
            xGrid = 0:2 * pi / numel(this.gridValues):2 * pi - 2 * pi / numel(this.gridValues);
            step = numel(this.gridValues) / noOfDesiredGridpoints;
            assert(mod(step, 1) == 0, ...
                'Number of function values has to be a multiple of noOfGridpoints');
            
            vals=this.gridValues(1:step:end)';
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
            if mod(numel(this.gridValues), 2) == 0
                warnStruct = warning('off', 'Truncate:TooFewCoefficients');
                noCoeffs = noCoeffs + 1;
            end
            fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'identity');
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
            p(2) = plot(gridPoints, this.gridValues, 'x');
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
            assert(this.enforcePdfNonnegative == f2.enforcePdfNonnegative);
            assert(numel(this.gridValues) == numel(f2.gridValues));
            warnStruct = warning('off', 'Normalization:notNormalized');
            f = FIGDistribution(this.gridValues.*f2.gridValues, this.enforcePdfNonnegative);
            warning(warnStruct);
        end
        
        function f = convolve(this, f2)
            % Multiplies two transformed fourier pdfs (returns transformed result)
            assert(this.enforcePdfNonnegative == f2.enforcePdfNonnegative);
            assert(numel(this.gridValues) == numel(f2.gridValues))
            convResult = cconv(this.gridValues, f2.gridValues, numel(f2.gridValues)) * 2 * pi / numel(f2.gridValues);
            % Due to numerical inaccuracies, negative values can
            % occur
            convResult(convResult < 0) = 0;
            f = FIGDistribution(convResult, this.enforcePdfNonnegative);
        end
        
        function f = truncate(this, noOfGridpoints)
            % Truncates Fourier series. Fills up if there are less coefficients
            % Expects number of complex coefficients (or sum of number of real
            % coefficients)
            assert(noOfGridpoints - 1 > 0, 'Number of coefficients must be an integer greater zero');
            step = numel(this.gridValues) / noOfGridpoints;
            % Truncation can void normalization
            warnStruct = warning('off', 'Normalization:notNormalized');
            if mod(step, 1) == 0
                f = FIGDistribution(this.gridValues(1:step:end),this.enforcePdfNonnegative);
            elseif numel(this.gridValues) < noOfGridpoints
                warning('Truncate:TooFewGridPoints', 'Less coefficients than desired, interpolate using Fourier while ensuring nonnegativity.');
                this.enforcePdfNonnegative = true;
                f = FIGDistribution.fromDistribution(this, noOfGridpoints, this.enforcePdfNonnegative);
            else
                warning('Truncate:DownsampleViaFourier', 'Cannot downsample directly. Transforming to Fourier to interpolate.');
                f = FIGDistribution.fromDistribution(this, noOfGridpoints, this.enforcePdfNonnegative);
            end
            warning(warnStruct);
        end
        
        function f = normalize(this)
            tol = 1e-4;
            integral = 2 * pi * mean(this.gridValues);
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
            f = FIGDistribution(f.gridValues/integral, this.enforcePdfNonnegative);
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
            fd = FourierDistribution.fromFunctionValues(this.gridValues', numel(this.gridValues), 'identity');
            fd = fd.shift(angle);
            gd = this;
            gd.gridValues = fd.pdfOnGrid(numel(this.gridValues))';
        end
    end
    
    methods(Static)
        function f = fromDistribution(distribution, noOfGridpoints, enforcePdfNonnegative)
            f = FIGDistribution.fromFunction(@(x)distribution.pdf(x), noOfGridpoints, enforcePdfNonnegative);
        end
        
        function f = fromFunction(fun, noOfCoefficients, enforcePdfNonnegative)
            gridPoints = 0:2 * pi / noOfCoefficients:2 * pi - 2 * pi / noOfCoefficients;
            gridValues = fun(gridPoints)';
            f = FIGDistribution(gridValues, enforcePdfNonnegative);
        end
        
        function gd = fromFunctionValues(fvals, noOfGridpoints, enforcePdfNonnegative)
            % Mostly useful for compatibility with FourierDistribution
            step = numel(fvals) / noOfGridpoints;
            assert(mod(step, 1) == 0, ...
                'Number of function values has to be a multiple of noOfGridpoints');
            fvals = fvals(1:step:end);
            gd = FIGDistribution(fvals, enforcePdfNonnegative);
        end
        
        
    end
    
    
end
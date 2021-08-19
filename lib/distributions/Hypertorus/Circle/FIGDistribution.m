classdef FIGDistribution < AbstractCircularDistribution & HypertoroidalGridDistribution
    % This density representation is used by the Fourier-interpreted grid
    % filter. The density is represented using function values on a grid.
    % The interpolation is performed using Fourier series. This class is
    % thus related to FourierDistribution
    
    % see Florian Pfaff, Kailai Li, and Uwe D. Hanebeck,
    % Fourier Filters, Grid Filters, and the Fourier-Interpreted Grid Filter,
    % Proceedings of the 22nd International Conference on Information Fusion (Fusion 2019), Ottawa, Canada, July, 2019.
    methods
        function this = FIGDistribution(gridValues, enforcePdfNonnegative)
            arguments
                gridValues {mustBeNonnegative}
                enforcePdfNonnegative logical = true
            end
            if isa(gridValues, 'AbstractCircularDistribution')
                error('You gave a distribution as the first argument. To convert distributions to a distribution in grid representation, use .fromDistribution');
            end
            % The normalization routine called in the constructor will be
            % the normalization routine of FIGDistribution, not the one
            % from HypertoroidalGridDistribution, so all is fine!
            this@HypertoroidalGridDistribution([],gridValues,enforcePdfNonnegative);
        end
        
        function p = pdf(this, xs, useSinc, sincRepetitions)
            arguments
                this (1,1) FIGDistribution
                xs (1,:) double
                useSinc (1,1) logical = false
                sincRepetitions (1,1) double = 5
            end
            if useSinc
                assert(mod(sincRepetitions,2)==1);
                step_size = 2*pi/numel(this.gridValues);
                range = -floor(sincRepetitions/2)*numel(this.gridValues):ceil(sincRepetitions/2)*numel(this.gridValues)-1;
                sincVals = sinc((xs/step_size)'-range);
                if this.enforcePdfNonnegative
                    p = sum(repmat(sqrt(this.gridValues'),[1,sincRepetitions]).*sincVals,2).^2;
                else
                    p = sum(repmat(this.gridValues',[1,sincRepetitions]).*sincVals,2);
                end
            else
                noCoeffs = numel(this.gridValues);
                if mod(numel(this.gridValues), 2) == 0
                    noCoeffs = noCoeffs + 1; % Extra coefficient because we add another bk with 0 (see FourierDistribution.fromFunctionValues)
                end
                if this.enforcePdfNonnegative
                    % Use square root to enforce nonnegative values
                    fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'sqrt');
                else
                    fd = FourierDistribution.fromFunctionValues(this.gridValues', noCoeffs, 'identity');
                end
                p = fd.pdf(xs);
            end
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
            arguments
                this (1,1) FIGDistribution
                n (1,1) double {mustBeInteger}
            end
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
            arguments
                this (1,1) FIGDistribution
            end
            arguments (Repeating)
                varargin
            end
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
            arguments
                this (1,1) FIGDistribution
                xa (1,:) double
            end
            % Value when interpreted as a pmf
            % Cannot store this after determining it because it is not a handle class
            p = (xa == this.getGrid()')' * this.gridValues;
        end
        
        function gridPoints = getGrid(this)
            arguments
                this (1,1) FIGDistribution
            end
            % For compatibility with other GridDistributions
            gridPoints = 0:2 * pi / numel(this.gridValues):2 * pi - 2 * pi / numel(this.gridValues);
        end

        function gridPoints = getGridPoint(this, indices)
            % To avoid passing all points if only one or few are needed.
            % Overload if .grid should stay empty
            arguments
                this (1,1) AbstractGridDistribution
                indices double {mustBeInteger} = []
            end
            if isempty(indices)
                gridPoints = this.getGrid();
            else
                gridPoints = (indices-1)*2*pi/numel(this.gridValues);
            end
        end
        
        function f = convolve(this, f2)
            arguments
                this (1,1) FIGDistribution
                f2 (1,1) FIGDistribution
            end
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
            arguments
                this (1,1) FIGDistribution
                noOfGridpoints (1,1) double {mustBeInteger,mustBePositive}
            end
            % Reduce number of grid points. Fills up if there are less coefficients
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
        
        function f = normalize(this, opt)
            arguments
                this (1,1) HypertoroidalGridDistribution
                opt.tol (1,1) double = 1e-2
                opt.warnUnnorm (1,1) logical = true
            end
            f = normalize@AbstractGridDistribution(this,tol=opt.tol,warnUnnorm=opt.warnUnnorm);
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
            arguments
                this (1,1) FIGDistribution
                angle (1,1) double
            end
            assert(isscalar(angle));
            fd = FourierDistribution.fromFunctionValues(this.gridValues', numel(this.gridValues), 'identity');
            fd = fd.shift(angle);
            gd = this;
            gd.gridValues = fd.pdfOnGrid(numel(this.gridValues))';
        end
        
        function [points, indices] = getClosestPoint(this, xa)
            arguments
                this (1,1) FIGDistribution
                xa (1,:) double
            end
            % Overload if class does not have .grid
            % + 1 because indexed like that, the indices start at 0
            indices = mod(round(xa/(2*pi/this.noOfGridPoints)),this.noOfGridPoints)+1;
            points = this.getGridPoint(indices); % Because Matlab starts indexing at 1
        end
    end
    
    methods(Static)
        function f = fromDistribution(distribution, noOfGridpoints, enforcePdfNonnegative)
            arguments
                distribution AbstractCircularDistribution
                noOfGridpoints (1,1) {mustBeInteger,mustBePositive}
                enforcePdfNonnegative logical = true
            end
            if isa(distribution,'FourierDistribution')
                warnStruct = warning('off','Truncate:TooFewCoefficients');
                fdToConv = distribution.truncate(noOfGridpoints);
                warning(warnStruct);
                valsOnGrid = real(ifft(ifftshift(fdToConv.c), 'symmetric')) * (numel(fdToConv.a) + numel(fdToConv.b));
                switch fdToConv.transformation
                    % FourierIdentity -> FIG with identity interpolation is
                    % precise if no negative values. If negative values, it
                    % is not precise.
                    % FourierIdentity -> FIG with sqrt interpolation is never
                    % precise
                    % FourierSqrt -> FIG with sqrt can be precise if the
                    % ifft of the values are all nonnegative. Then, .^2 is
                    % precisely undone by sqrt(...). A transformation to
                    % identity can help! By convolving vector with itself,
                    % we get the (longer) coefficients vector precisely
                    % describing the same density
                    % FourierSqrt -> FIG with identity can be precise if
                    % sufficient coefficients are used so that the square
                    % of the trigonometric polynomial (i.e., the density)
                    % can be represented using the Fourier coefficients.
                    % This is related to the previous case
                    % See also the test cases
                    case 'identity'
                        if any(valsOnGrid<0)
                            warning('FourierToFIG:ImpreciseId','This is an inaccurate transformation because values smaller than 0 came up. These values were increased to 0.');
                            valsOnGrid = max(valsOnGrid,0);
                        end
                        if enforcePdfNonnegative
                            warning('FourierToFIG:OtherInterpolation','This is generally an inaccurate transformation because the interpolation is different. Set enforcePdfNonnegative to false false when transforming FourierDensities with transformation ''identity'' to obtain the same density again.');
                        end
                    case 'sqrt'
                        if any(valsOnGrid<0) && enforcePdfNonnegative
                            warning('FourierToFIG:ImpreciseSqrt','This is generally an inaccurate transformation because values smaller than 0 came up. Consider transforming the distribution first to be in the identity transformation.');
                        elseif ~enforcePdfNonnegative && noOfGridpoints < (2*(numel(distribution.a) + numel(distribution.b))-1)
                            % If too few coefficients for exact
                            % represetation, throw warning (must consider
                            % *original* Fourierdistribution, i.e.,
                            % the variable "distribution"
                            warning('FourierToFIG:OtherInterpolationWithInsufficientCoeffs','This may be an inaccurate transformation because the interpolation is different. For n coefficients, using 2n-1 grid points and setting enforcepdfNonnegative to false will result in a perfect conversion. Setting it to true will not always ensure an identical density as negative values in the interpolation of the square root of the density get lost when transforming to FIGDistribution');
                        end
                        valsOnGrid = valsOnGrid.^2;
                    otherwise
                        error('Transformation unsupported');
                end
                f = FIGDistribution(valsOnGrid, enforcePdfNonnegative);
            else
                f = FIGDistribution.fromFunction(@(x)distribution.pdf(x), noOfGridpoints, enforcePdfNonnegative);
            end
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
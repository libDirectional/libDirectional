classdef FourierDistribution < AbstractCircularDistribution
    % Used to represent circular densities with Fourier
    % series
    %#ok<*PROP>
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Multimodal Circular Filtering Using Fourier Series
    % Proceedings of the 18th International Conference on Information Fusion (Fusion 2015),
    % Washington D. C., USA, July 2015.
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Nonlinear Prediction for Circular Filtering Using Fourier Series
    % Proceedings of the 19th International Conference on Information Fusion (Fusion 2016),
    % Heidelberg, Germany, July 2016.
    
    properties
        a (1,:) double {mustBeReal}
        b (1,:) double {mustBeReal}
        transformation char = 'sqrt'
    end
    
    methods
        function this = FourierDistribution(a, b, transformation)
            if isa(a, 'AbstractCircularDistribution')
                error('You gave a distribution as the first argument. To convert distributions to a distribution in Fourier representation, use .fromDistribution');
            end
            assert(isreal(a) && isreal(b), 'Coefficient vecotrs must be real when using FourierDistribution. To use complex vectors, use .fromComplex or use HypertoroidalFourierDistribution.');
            assert(length(b) == (length(a) - 1), 'Coefficients have incompatible lengths');
            if nargin == 2 % Square root of density is standard case
                this.transformation = 'sqrt';
            else
                this.transformation = transformation;
            end
            this.a = a;
            this.b = b;
            % Check if normalized. If not: Normalize!
            this = this.normalize;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (1 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa, 1) == 1);
            % Evaluate actual pdf at xa (transformations need to be performed)
            val = value(this, xa);
            switch this.transformation
                case 'sqrt'
                    assert(all(imag(val) < 0.0001));
                    p = real(val).^2;
                case 'identity'
                    p = val;
                case 'log'
                    normConstLog = 1 / integral(@(x)exp(this.value(x)), 0, 2*pi); % Because this is a value class, we can only store values as attributes in the constructor
                    p = exp(val) * normConstLog;
                otherwise
                    error('Transformation not recognized or unsupported');
            end
        end
        
        function result = integral(this, l, r)
            % Calculates the integral of the pdf from l to r analytically
            % if possible, fall back to numerical calculation by default
            %
            % Parameters:
            %   l (scalar)
            %       left bound of integral, default 0
            %   r (scalar)
            %       right bound of integral, default 2*pi
            % Returns:
            %   result (scalar)
            %       value of the integral
            
            % This function uses the cdf. It is not implemented the other way
            % around because integral is only guaranteed to take scalar
            % arguments, whereas cdf always has to support vector-valued
            % inputs in libDirectional.
            if nargin < 2
                l = 0;
            end
            if nargin < 3
                r = 2 * pi;
            end
            if l <= r
                result = floor((r - l)/(2 * pi)); %each full 2pi interval contributes 1 to the result because of normalization
                r = mod(r, 2*pi);
                l = mod(l, 2*pi);
                if l <= r
                    result = result + this.cdf(r, l);
                else
                    result = result + this.cdf(r, 0) + this.cdf(2*pi, l); % Closer to integralNumerical than 1 - this.cdf(l,r) if not perfectly normalized
                end
            else
                result = -this.integral(r, l);
            end
        end
        
        function val = cdf(this, xa, startingPoint)
            % Evaluate cumulative distribution function
            %
            % Parameters:
            %   xa (1 x n)
            %       points where the cdf should be evaluated
            %   startingPoint (scalar)
            %       point where the cdf is zero (starting point can be
            %       [0,2pi) on the circle, default 0
            % Returns:
            %   val (1 x n)
            %       cdf evaluated at columns of xa
            if nargin <= 2, startingPoint = 0;end
            switch this.transformation
                case 'identity'
                    fd = this;
                case 'sqrt' % Transform to identity
                    fd = this.transformViaCoefficients('square', 4*length(this.a)-3);
                case 'log'
                    normConstLog = 1 / integral(@(x)exp(this.value(x)), 0, 2*pi); % Calculate factor once, do not use pdf to prevent calculating factor multiple times
                    startingPoint = mod(startingPoint, 2*pi);
                    xa = mod(xa, 2*pi);
                    if xa < startingPoint
                        val = 1 - normConstLog * arrayfun(@(xCurr)integral(@(x)exp(this.value(x)), xCurr, startingPoint), xa);
                    else
                        val = normConstLog * arrayfun(@(xCurr)integral(@(x)exp(this.value(x)), startingPoint, xCurr), xa);
                    end
                    return
                otherwise % Fall back to numerical integration for this transformation.
                    val = arrayfun(@(xCurr)this.integralNumerical(xCurr, startingPoint), xa);
                    return
            end
            xa = mod(xa-startingPoint, 2*pi) + startingPoint;
            c = fd.c;
            c0 = c((length(c) + 1)/2);
            cnew = fd.c ./ (1i * (-length(fd.b):length(fd.b))); %Calculate coefficients != 0 for antiderivative
            % To avoid unwanted normalization by the constructor, we set
            % c0=1/(2*pi). We do not have to address this further as +c to the
            % indefinite integral does not change the value of the definite integral.
            cnew((length(cnew) + 1)/2) = 1 / (2 * pi);
            fdInt = FourierDistribution.fromComplex(cnew, 'identity');
            val = fdInt.value(xa) - fdInt.value(startingPoint) + c0 * (xa - startingPoint);
        end
        
        function p = value(this, xa)
            arguments
                this (1,1) FourierDistribution
                xa (1,:) double
            end
            % Evalute current Fourier series without undoing transformations
            p = this.a(1)/2 + this.a(2:end)*cos(xa.*(1:length(this.a) - 1)') + this.b*sin(xa.*(1:length(this.b))');
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
            if n == 0
                m = 1;
            elseif n < 0
                m = conj(trigonometricMoment(this, -n));
            else
                switch this.transformation
                    case 'sqrt'
                        fdtmp = this.transformViaCoefficients('square', 4*length(this.b)+1);
                        atmp = fdtmp.a;
                        btmp = fdtmp.b;
                    case 'identity'
                        atmp = this.a;
                        btmp = this.b;
                    otherwise
                        m = trigonometricMomentNumerical(this, n);
                        return
                end
                if n > length(btmp)
                    m = 0;
                else
                    m = pi * conj(atmp(n+1)-1i*btmp(n));
                end
            end
        end
        
        function complexCoeffs = c(this)
            % Get complex coefficients
            complexCoeffs = NaN(1, length(this.a)+length(this.b), 'like', 1i);
            complexCoeffs((length(complexCoeffs) + 1)/2) = this.a(1) / 2;
            complexCoeffs((length(complexCoeffs) + 3)/2:end) = (this.a(2:end) - 1i * this.b) / 2;
            % We know pdf is real, so we can use complex conjugation
            complexCoeffs(1:(length(complexCoeffs) - 1)/2) = conj(fliplr(complexCoeffs((length(complexCoeffs) + 3)/2:end)));
        end
        
        function f = multiply(this, f2, noOfCoefficients)
            assert(isa(f2, 'FourierDistribution'));
            if nargin == 2
                noOfCoefficients = numel(this.a) + numel(this.b);
            end
            % Multiplies two transformed fourier pdfs (returns transformed result)
            if ~strcmp(this.transformation, f2.transformation)
                error('Multiply:differentTransformations', 'Transformations do not match, transform before using multiply');
            end
            if strcmp(this.transformation, 'log')
                % Normalization is not performed for log
                f = FourierDistribution(this.a+f2.a, this.b+f2.b, 'log');
                warning('Multiply:NotNormalizing', 'Not performing normalization when using log transformation.');
            elseif strcmp(this.transformation, 'identity') || strcmp(this.transformation, 'sqrt')
                % Calculate unnormalized result
                if noOfCoefficients <= numel(this.a) + numel(this.b)
                    c = conv(this.c, f2.c, 'same');
                else
                    c = conv(this.c, f2.c, 'full');
                end
                % Normalization is performed in constructor. Temporarily
                % disabling warning.
                warnStruct = warning('off', 'Normalization:notNormalized');
                f = FourierDistribution.fromComplex(c, this.transformation);
                warning(warnStruct);
                f.truncate(noOfCoefficients);
            else
                error('Multiply:unsupportedTransformation', 'Transformation not recognized or unsupported');
            end
            
        end
        
        function f = normalize(this, opt)
            arguments
                this (1,1) FourierDistribution
                opt.tol (1,1) double = 1e-4
                opt.warnUnnorm (1,1) logical = true
            end
            % Normalize Fourier density while taking its type into account
            % Avoid use of constructor to avoid normalizing twice
            f = this;
            switch this.transformation
                case 'sqrt'
                    % Calculate normalization factor and return normalized
                    % result
                    % Square root calculated later to use norm and not squared norm
                    sqrtOfa0 = norm([this.a(1)/sqrt(2),this.a(2:end),this.b]);
                    % Alternative with c: a0 = 2 * norm(this.c)^2; 
                    if sqrtOfa0 < 0
                        error('Normalization:negative', 'a0 is negative, this usually points to a user error');
                    elseif (sqrtOfa0^2 < 1e-200) % Tolerance has to be that low to avoid unnecessary errors on multiply
                        error('Normalization:almostZero', 'a0 is too close to zero, this usually points to a user error');
                    elseif (abs(sqrtOfa0^2-1/pi) > opt.tol)
                        if opt.warnUnnorm
                            warning('Normalization:notNormalized', 'Coefficients apparently do not belong to normalized density. Normalizing...');
                        end
                        f.a = this.a / (sqrtOfa0*sqrt(pi));
                        f.b = this.b / (sqrtOfa0*sqrt(pi));
                    end
                case 'identity'
                    % Calculate normalization factor and return normalized
                    % result
                    a0 = this.a(1);
                    if a0 < 0
                        warning('Normalization:negative', 'a0 is negative. This can either be caused by a user error or due to negativity caused by non-square rooted version');
                    elseif abs(a0) < 1e-200 % Tolerance has to be that low to avoid unnecessary errors on multiply
                        error('Normalization:almostZero', 'a0 is too close to zero, this usually points to a user error');
                    elseif (abs(a0-1/pi) > 1e-4)
                        warning('Normalization:notNormalized', 'Coefficients apparently do not belong to normalized density. Normalizing...');
                    else % Is already normalized, just return original density
                        return
                    end
                    f.a = this.a / (a0 * pi);
                    f.b = this.b / (a0 * pi);
                otherwise
                    warning('Normalization:cannotTest', 'Unable to test if normalized');
            end
        end
        
        function f = convolve(this, f2, noOfCoefficients)
            % Calculates convolution of two Fourier series
            % Expects number of complex coefficients (or sum of number of real
            % coefficients) to know how many points to sample for FFT
            if ~strcmp(this.transformation, f2.transformation)
                error('Convolve:differentTransformations', 'Transformations do not match, transform before using convolve');
            end
            if nargin == 2, noOfCoefficients = length(this.a) + length(this.b);end
            c1 = this.c;
            c2 = f2.c;
            switch this.transformation
                case 'sqrt'
                    % Calculate convolution in an exact fashion by first
                    % obtaining coefficients for the identity and then using
                    % the Hadamard product
                    cConv = 2 * pi * conv(c1, c1) .* conv(c2, c2);
                    ftmp = FourierDistribution.fromComplex(cConv, 'identity');
                    % Calculate coefficients for sqrt
                    f = ftmp.transformViaFFT('sqrt', noOfCoefficients);
                case 'log'
                    % Calculate function values and then calculate cyclic
                    % convolution via fft, this is already an approximation
                    fvals1 = exp(ifft(ifftshift(c1), 'symmetric')*length(c1));
                    fvals2 = exp(ifft(ifftshift(c2), 'symmetric')*length(c2));
                    ctmp = fftshift(fft(fvals1).*fft(fvals2));
                    ctmp = ctmp / length(ctmp);
                    ftmp = FourierDistribution.fromComplex(ctmp, 'identity');
                    % Calculate coefficients for log
                    f = ftmp.transformViaFFT('log', noOfCoefficients);
                case 'identity'
                    % This can be done in an exact fashion
                    cConv = 2 * pi * c1 .* c2;
                    ftmp = FourierDistribution.fromComplex(cConv, 'identity');
                    f = ftmp.truncate(noOfCoefficients);
                otherwise
                    error('Convolve:unsupportedTransformation', 'Transformation not recognized or unsupported');
            end
        end
        
        function f = truncate(this, noOfCoefficients)
            % Truncates Fourier series. Fills up if there are less coefficients
            % Expects number of complex coefficients (or sum of number of real
            % coefficients)
            arguments
                this (1,1) FourierDistribution
                noOfCoefficients (1,1) double {mustBePositive, mustBeInteger}
            end
            f = this;
            assert((noOfCoefficients - 1 > 0) && (mod(noOfCoefficients-1, 2) == 0), 'Invalid number of coefficients, number has to be odd');
            if ((noOfCoefficients + 1) / 2) <= length(this.a)
                f.a = this.a(1:((noOfCoefficients + 1) / 2));
                f.b = this.b(1:((noOfCoefficients - 1) / 2));
            else
                warning('Truncate:TooFewCoefficients', 'Less coefficients than desired, filling up with zeros')
                diff = (noOfCoefficients + 1) / 2 - length(this.a);
                f.a = [this.a, zeros(1, diff)];
                f.b = [this.b, zeros(1, diff)];
            end
            % Truncation can void normalization if transformation is not
            % identity
            if ~strcmp(f.transformation, 'identity')
                % Disable warning as we expect normalization to be
                % necessary
                warnStruct = warning('off', 'Normalization:notNormalized');
                f = f.normalize;
                warning(warnStruct);
            end
        end
        
        function [vals, xGrid] = pdfOnGrid(this, noOfGridpoints)
            arguments
                this (1,1) FourierDistribution
                noOfGridpoints (1,1) double = numel(this.a) + numel(this.b)
            end
            % Calculated function on grid from 0 to
            % 2*pi-2*pi/noOfGridpoints
            assert(mod(noOfGridpoints, 2) == 1, 'Evaluation on grid only supported for an odd number of grid points.')
            skippingFactor = (numel(this.a) + numel(this.b) - 1) / (noOfGridpoints - 1);
            assert(noOfGridpoints >= (numel(this.a) + numel(this.b)) || (floor(skippingFactor)) == (skippingFactor), ...
                'Number of grid points needs to be either higher than the number of coefficients or the grid points need to be obtainable via skipping.');
            if noOfGridpoints > (numel(this.a) + numel(this.b)) % Fill up with zeros to generate more function values via ifft
                warningSetting = [warning('off', 'Truncate:TooFewCoefficients'), warning('off', 'Normalization:cannotTest')];
                this = this.truncate(noOfGridpoints);
                warning(warningSetting);
            end
            switch this.transformation
                case 'identity'
                    vals = real(ifft(ifftshift(this.c), 'symmetric')) * (numel(this.a) + numel(this.b));
                case 'sqrt'
                    vals = (real(ifft(ifftshift(this.c), 'symmetric')) * (numel(this.a) + numel(this.b))).^2;
                case 'log'
                    vals = exp(real(ifft(ifftshift(this.c), 'symmetric'))*(numel(this.a) + numel(this.b))) / integral(@(x)exp(this.value(x)), 0, 2*pi);
                otherwise
                    error('Transformation not recognized or unsupported');
            end
            if noOfGridpoints < (numel(this.a) + numel(this.b)) % Use skipping if fewer gridpoints are desired
                vals = vals(1:skippingFactor:end);
            end
            if nargout == 2
                xGrid = 0:2 * pi / noOfGridpoints:2 * pi - 2 * pi / noOfGridpoints;
            end
        end
        
        function f = transformViaCoefficients(this, desiredTransformation, noOfCoefficients)
            % Calculates transformations using Fourier coefficients
            if nargin == 2, noOfCoefficients = length(this.a) + length(this.b);end
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
                    c = conv(this.c, this.c);
                    cTrunc = c(max((length(c) + 1)/2-(noOfCoefficients - 1)/2, 1):min((length(c) + 1)/2+(noOfCoefficients - 1)/2, end));
                    f = FourierDistribution.fromComplex(cTrunc, newTrans);
                otherwise
                    error('Desired transformation not supported via coefficients');
            end
            f = f.truncate(noOfCoefficients);
        end
        
        function f = transformViaFFT(this, desiredTransformation, noOfCoefficients, truncateAtEnd)
            % Calculates transformation of Fourier series via FFT.
            % Expects number of complex coefficients (or sum of number of real
            % coefficients). This is done by calculating the values of the
            % function on the grid (using inverse FFT in pdfOnGrid),
            % applying the transformation, and then using FFT again to
            % obtain the Fourier coefficients. If truncateAtEnd is set to
            % true, more coefficients are calculated if the original
            % density features more. If set to false, entries of the
            % inverse FFT result are discarded.
            if nargin == 2, noOfCoefficients = length(this.a) + length(this.b);end
            if nargin <= 3, truncateAtEnd = true;end
            if truncateAtEnd
                noOfGridpoints = max(noOfCoefficients, length(this.a)+length(this.b));
            else
                noOfGridpoints = noOfCoefficients;
            end
            switch this.transformation
                case 'identity'
                    f = FourierDistribution.fromFunctionValues(this.pdfOnGrid(noOfGridpoints), noOfGridpoints, desiredTransformation);
                case 'sqrt'
                    assert(strcmp(desiredTransformation, 'square'), 'Transformation:cannotCombine', ...
                        'For series of already transformed distributions, only transformation back to identity is supported');
                    f = FourierDistribution.fromFunctionValues(this.pdfOnGrid(noOfGridpoints), noOfGridpoints, 'identity');
                case 'log'
                    assert(strncmp(desiredTransformation, 'power', 3), 'Transformation:cannotCombine', ...
                        'For series of already transformed distributions, only transformation back to identity is supported');
                    f = FourierDistribution.fromFunctionValues(this.pdfOnGrid(noOfGridpoints), noOfGridpoints, 'identity');
                otherwise
                    error('Transformation not recognized or unsupported');
            end
            if truncateAtEnd
                f = f.truncate(noOfCoefficients);
            end
        end
        
        function f = transformViaVM(this, desiredTransformation, noOfCoefficients)
            % Calculates transformation of Fourier series via approximation with a von Mises distribution
            % Expects number of complex coefficients (or sum of number of real
            % coefficients)
            if ~strcmp(this.transformation, 'identity')
                error('Transformations of already transformed density via von Mises not supported');
            end
            if nargin == 2, noOfCoefficients = length(this.a) + length(this.b);end
            assert((noOfCoefficients - 1 > 0) && (mod(noOfCoefficients-1, 2) == 0), 'Invalid number of coefficients, number has to be odd');
            vmEquivalent = VMDistribution.fromMoment(this.a(2)*pi+1i*this.b(1)*pi);
            switch desiredTransformation
                case 'sqrt'
                    f = FourierDistribution.fromDistribution(vmEquivalent, noOfCoefficients, 'sqrt');
                case 'log'
                    f = FourierDistribution.fromDistribution(vmEquivalent, noOfCoefficients, 'log');
                otherwise
                    error('Transformation not recognized or unsupported');
            end
        end
        
        function dist = hellingerDistanceNumerical(this, other)
            % Numerically calculates the Hellinger distance to another
            % distribution.
            %
            % Parameters:
            %   other (AbstractHypertoroidalDistribution)
            %       distribution to compare with
            % Returns:
            %   dist (scalar)
            %       hellinger distance of this distribution to other distribution
            if isa(other, 'FourierDistribution') && strcmp(this.transformation, 'sqrt') && strcmp(other.transformation, 'sqrt')
                % For square root, we can use 1-int(sqrt(f(x)*g(x)),0,2*pi)
                cconv = conv(this.c, other.c);
                dist = 1 - 2 * pi * real(cconv((numel(cconv) + 1)/2));
            else % Fall back for other transformations
                dist = hellingerDistanceNumerical@AbstractHypertoroidalDistribution(this, other);
            end
        end
        
        function fd = shift(this, angle)
            % Shift distribution by the given angle
            %
            % Parameters:
            %   shiftAngles (scalar)
            %       angle to shift by
            % Returns:
            %   fd (FourierDistribution)
            %       shifted distribution
            assert(isscalar(angle));
            anew = [this.a(1), this.a(2:end) .* cos(-(1:length(this.b))*angle) + this.b(1:end) .* sin(-(1:length(this.b))*angle)];
            bnew = this.b(1:end) .* cos(-(1:length(this.b))*angle) - this.a(2:end) .* sin(-(1:length(this.b))*angle);
            fd = FourierDistribution(anew, bnew, this.transformation);
        end
    end
    
    methods(Static)
        function f = fromComplex(c, transformation)
            arguments
                c (1,:) double {mustBeNonempty}
                transformation char
            end
            % Create density from complex coefficients
            assert(size(c, 2) >= size(c, 1), 'c is expected to be a row vector');
            assert(abs(c((length(c) + 1)/2)) > 0, 'c0 is zero, cannot normalize to valid density.')
            % Need to double values. Adding with flipping variant to
            % neither favor negative nor positive indicies
            % (although c_k=conj(c_-k) should hold)
            tmp = c + fliplr(c);
            % Ensure all coefficients to be real due to numerical imprecision
            a = real(tmp(((length(c) + 1) / 2):end)); %a_0..a_n
            % Double the imaginary part using - since c_k=conj(c_-k) holds
            tmp = c - fliplr(c);
            b = -imag(tmp(((length(c) + 3) / 2):end)); %b_1..b_n
            f = FourierDistribution(a, b, transformation);
        end
        
        function f = fromFunction(fun, noOfCoefficients, desiredTransformation)
            % Creates Fourier distribution from function
            % Function must be able to take vector arguments
            assert(isa(fun, 'function_handle'));
            if nargin == 2, desiredTransformation = 'sqrt';end
            xvals = 0:2 * pi / noOfCoefficients:2 * pi - 2 * pi / noOfCoefficients;
            fvals = fun(xvals);
            assert(isequal(size(fvals), [1, noOfCoefficients]), 'Size of output of pdf is incorrect. Please ensure that the pdf returns only one scalar.');
            f = FourierDistribution.fromFunctionValues(fvals, noOfCoefficients, desiredTransformation);
        end
        
        function f = fromFunctionValues(fvals, noOfCoefficients, desiredTransformation)
            % Creates Fourier distribution from function values
            % Assumes fvals are not yet transformed, use custom if they already
            % are transformed
            assert((noOfCoefficients - 1 > 0) && (mod(noOfCoefficients-1, 2) == 0), 'Invalid number of coefficients');
            
            switch desiredTransformation
                case 'sqrt'
                    fvals = sqrt(fvals);
                case 'log'
                    fvals = log(fvals);
                case 'identity' %keep them unchanged
                case 'custom' %already transformed
                otherwise
                    error('Transformation not recognized or unsupported by transformation via FFT');
            end
            transformed = fftshift(fft(fvals)) / length(fvals);
            if mod(length(fvals), 2) == 0 % An additional a_k can be obtained
                aLast = transformed(1);
                transformed(1) = [];
            end
            coeffsTruncated = transformed(max((length(transformed) + 1)/2-(noOfCoefficients - 1)/2, 1):min((length(transformed) + 1)/2+(noOfCoefficients - 1)/2, end));
            ftmp = FourierDistribution.fromComplex(coeffsTruncated, desiredTransformation);
            if mod(length(fvals), 2) == 0
                ftmp.a = [ftmp.a, aLast];
                ftmp.b = [ftmp.b,0]; % No information is available. Adding it as zero to preserve the lengths
            end
            f = ftmp.truncate(noOfCoefficients);
        end
        
        function f = fromDistribution(distribution, noOfCoefficients, desiredTransformation)
            % Creates Fourier distribution from a different distribution
            assert(isa(distribution, 'AbstractCircularDistribution'), 'First argument has to be a circular distribution.');
            assert((noOfCoefficients > 0) && (mod(noOfCoefficients-1, 2) == 0), 'Invalid number of coefficients, number has to be odd');
            if nargin == 2
                desiredTransformation = 'sqrt';
            end
            lastk = (noOfCoefficients - 1) / 2;
            switch class(distribution)
                case 'VMDistribution'
                    switch desiredTransformation
                        case 'sqrt'
                            a = 2 / sqrt(2*pi*besseli(0, distribution.kappa)) * besseli(0:lastk, 0.5*distribution.kappa);
                        case 'identity'
                            a = [1 / pi, ...
                                1 / (pi * besseli(0, distribution.kappa)) * besseli(1:lastk, distribution.kappa)];
                        case 'log'
                            a = [-2 * log(2*pi*besseli(0, distribution.kappa)), distribution.kappa, zeros(1, lastk-1)];
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                    f = FourierDistribution(a, zeros(1, length(a)-1), desiredTransformation);
                    if ~(distribution.mu == 0)
                        f = f.shift(distribution.mu);
                    end
                case 'WNDistribution'
                    switch desiredTransformation
                        case 'sqrt'
                            warning('Conversion:NoFormulaSqrt', 'No explicit formula available, using FFT to get sqrt');
                            f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
                        case 'identity'
                            a = [1 / pi, ...
                                arrayfun(@(k)exp(-distribution.sigma^2*k^2/2)/pi, 1:lastk)];
                            f = FourierDistribution(a, zeros(1, length(a)-1), desiredTransformation);
                        case 'log'
                            % Using first 1000 components to approximate
                            % log of euler function/q-pochhammer
                            logeuler = sum(log(1-exp(-distribution.sigma^2*(1:1000))));
                            a0 = 2 * (-log(2*pi) + logeuler);
                            a1end = arrayfun(@(k)2*(-1)^(k)/k*(exp(0.5*k*distribution.sigma^2) / (1 - exp(k*distribution.sigma^2))), 1:lastk);
                            a = [a0, a1end];
                            f = FourierDistribution(a, zeros(1, length(a)-1), desiredTransformation);
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                    if ~(distribution.mu == 0) && ~strcmp(desiredTransformation, 'sqrt')
                        f = f.shift(distribution.mu);
                    end
                case 'WCDistribution'
                    switch desiredTransformation
                        case 'sqrt'
                            warning('Conversion:ApproximationHypergeometric', 'The implementation of the regularized hypergeometric function may not be accurate numerically. This can lead to unnormalized densities');
                            noOfSummands = 1000;
                            % Cannot use log of gamma (which would be
                            % better numerically) as MATLAB currently does
                            % not support negative values for gammaln
                            afun = @(n, k)gamma(n+1/2).^2 .* sech(distribution.gamma/2).^(2 .* n) ./ (gamma(1-k+n) .* gamma(1+k+n));
                            a = NaN(1, lastk+1);
                            for k = 0:lastk
                                vals = afun(0:noOfSummands, k);
                                a(k+1) = sum(vals(~isnan(vals))) * sqrt(2/(pi^3)*tanh(distribution.gamma/2));
                            end
                        case 'identity'
                            a = arrayfun(@(k)exp(-k*distribution.gamma)/pi, 0:lastk);
                        case 'log'
                            a = [2 * log((1 - exp(-2*distribution.gamma))/(2 * pi)), ...
                                arrayfun(@(k)2*exp(-k*distribution.gamma)/k, 1:lastk)];
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                    f = FourierDistribution(a, zeros(1, length(a)-1), desiredTransformation);
                    if ~(distribution.mu == 0)
                        f = f.shift(distribution.mu);
                    end
                case 'WEDistribution'
                    switch desiredTransformation
                        case 'sqrt'
                            afun = @(k)2 * distribution.lambda^(3 / 2) ./ (pi * sqrt(exp(2*pi*distribution.lambda)-1) * (4 * k.^2 + distribution.lambda^2)) ...
                                * (exp(pi*distribution.lambda) - 1);
                            bfun = @(k)4 * k * distribution.lambda^(1 / 2) ./ (pi * sqrt(exp(2*pi*distribution.lambda)-1) * (4 * k.^2 + distribution.lambda^2)) ...
                                * (exp(distribution.lambda*pi) - 1);
                            a = afun(0:lastk);
                            b = bfun(1:lastk);
                        case 'identity'
                            a = [1 / pi, ...
                                arrayfun(@(k)1/pi*distribution.lambda^2/(distribution.lambda^2 + k^2), 1:lastk)];
                            b = arrayfun(@(k)1/pi*distribution.lambda*k/(distribution.lambda^2 + k^2), 1:lastk);
                        case 'log'
                            a = [-2 * pi * distribution.lambda - 2 * log(1-exp(-2*pi*distribution.lambda)) + 2 * log(distribution.lambda), ...
                                zeros(1, lastk)];
                            b = arrayfun(@(k)2*distribution.lambda/k, 1:lastk);
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                    f = FourierDistribution(a, b, desiredTransformation);
                case 'WLDistribution'
                    switch desiredTransformation
                        case 'sqrt'
                            warning('Conversion:NoFormulaSqrt', 'No explicit formula available, using FFT to get sqrt');
                            f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
                        case 'identity'
                            a0 = 1 / pi;
                            a1end = arrayfun( ...
                                @(k)1/pi*distribution.kappa^2*distribution.lambda^2*(k^2 + distribution.lambda^2)/ ...
                                ((distribution.lambda^2 * distribution.kappa^2 + k^2) * (distribution.kappa^2 * k^2 + distribution.lambda^2)), ...
                                1:lastk);
                            b = arrayfun( ...
                                @(k)1/pi*distribution.kappa*distribution.lambda^3*k*(1 - distribution.kappa^2)/ ...
                                ((distribution.lambda^2 * distribution.kappa^2 + k^2) * (distribution.kappa^2 * k^2 + distribution.lambda^2)), ...
                                1:lastk);
                            f = FourierDistribution([a0, a1end], b, desiredTransformation);
                        case 'log'
                            warning('Conversion:NoFormulaLog', 'No explicit formula available, using FFT to get log');
                            f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                case 'CircularUniformDistribution'
                    switch desiredTransformation
                        case 'sqrt'
                            a = [sqrt(2/pi), zeros(1, lastk)];
                        case 'identity'
                            a = [1 / pi, zeros(1, lastk)];
                        case 'log'
                            a = [-2 * log(2*pi), zeros(1, lastk)];
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                    f = FourierDistribution(a, zeros(1, length(a)-1), desiredTransformation);
                case 'WDDistribution'
                    warning('Conversion:WDMatching', 'Approximating WD by matching moments');
                    switch desiredTransformation
                        case 'identity'
                            ctmp = arrayfun(@(i)distribution.trigonometricMoment(i), 1:lastk);
                            c = [fliplr(ctmp), 1, conj(ctmp)]/(2*pi);
                        case 'sqrt'
                            % If we normalize anyways and we consider
                            % equally-weighted samples, we do not need to
                            % take take the square roots since it a
                            % constant that can be pulled out. If they are
                            % not equally weighted, we *must* take the
                            % square root!
                            distribution.w = sqrt(distribution.w);
                            ctmp = arrayfun(@(i)distribution.trigonometricMoment(i), 1:lastk);
                            c = [fliplr(ctmp), distribution.trigonometricMoment(0), conj(ctmp)]/sqrt(noOfCoefficients*2*pi);
                        case 'log'
                            distribution.w = log(distribution.w);
                            ctmp = arrayfun(@(i)distribution.trigonometricMoment(i), 1:lastk);
                            % 2*pi: Factor between coefficients and
                            % trigonometric moments
                            c = [fliplr(ctmp), distribution.trigonometricMoment(0), conj(ctmp)]/(2*pi);
                            warning('The result is certainly unnormalized');
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                    f = FourierDistribution.fromComplex(c, desiredTransformation);
                    
                case 'CircularMixture'
                    switch desiredTransformation
                        case 'sqrt'
                            warning('Conversion:NoFormulaSqrt', 'No explicit formula available, using FFT to get sqrt');
                            f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
                        case 'identity'
                            a = zeros(1, lastk+1);
                            b = zeros(1, lastk);
                            for i = 1:length(distribution.dists)
                                fCurr = FourierDistribution.fromDistribution(distribution.dists{i}, noOfCoefficients, desiredTransformation);
                                a = a + fCurr.a * distribution.w(i);
                                b = b + fCurr.b * distribution.w(i);
                            end
                            f = FourierDistribution(a, b, 'identity');
                        case 'log'
                            warning('Conversion:NoFormulaLog', 'No explicit formula available, using FFT to get log');
                            f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
                        otherwise
                            error('Transformation not recognized or unsupported');
                    end
                case 'GvMDistribution'
                    if strcmp(desiredTransformation, 'log')
                        if numel(distribution.kappa) * 2 + 1 > noOfCoefficients
                            warning('Conversion:GvMTooFewCoefficients', 'Converting a GvM distribution to a Fourier distribution with log transformation and using too few coefficients severely impedes the performance of the approximation.');
                        end
                        a = zeros(1, lastk+1);
                        b = zeros(1, lastk);
                        a(1:min(end, numel(distribution.kappa)+1)) = [1;distribution.kappa(1:min(end, numel(b))) ...
                            .* cos(-(1:min(numel(b), numel(distribution.kappa)))'.*distribution.mu(1:min(end, numel(b))))];
                        b(1:min(end, numel(distribution.kappa))) = -distribution.kappa(1:min(end, numel(b))) ...
                            .* sin(-(1:min(numel(b), numel(distribution.kappa)))'.*distribution.mu(1:min(end, numel(b))));
                        f = FourierDistribution(a, b, 'log');
                    else
                        warning('Conversion:NoFormula', 'No explicit formula available, using FFT to get transformation');
                        f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
                    end
                case 'FIGDistribution'
                    % This also works if the number of coefficients is
                    % not identical to the number of grid points because
                    % padding and truncation is done automatically.
                    f = FourierDistribution.fromFunctionValues(distribution.gridValues', noOfCoefficients, desiredTransformation);
                otherwise
                    warning('Conversion:NoFormula', 'No explicit formula available, using FFT to get transformation');
                    f = FourierDistribution.fromFunction(@distribution.pdf, noOfCoefficients, desiredTransformation);
            end
        end
        
        function fd = fromSamples(samples, noOfCoefficients, desiredTransformation)
            if nargin == 2
                desiredTransformation = 'sqrt';
            end
            wd = WDDistribution(samples);
            fd = FourierDistribution.fromDistribution(wd, ...
                noOfCoefficients, desiredTransformation);
        end
    end
end

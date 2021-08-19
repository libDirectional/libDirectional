classdef ToroidalFourierDistribution < HypertoroidalFourierDistribution & AbstractToroidalDistribution
    % Distribution on the Torus using Fourier series
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Multivariate Angular Filtering Using Fourier Series
    % Journal of Advances in Information Fusion, 11(2):206-226, December 2016.
    
    methods
        function this = ToroidalFourierDistribution(C, transformation)
            if nargin == 1, transformation = 'sqrt';
            end
            if any(size(C) == 1) % Prevent C from being interpreted as being from a 1D series
                warning('ToroidalFourierDistribution:VectorGiven', 'ToroidalFourierDistributions with only one coefficient in a dimension are not allowed, filling up to 3 in each dimension.');
                C = blkdiag(0, C, 0);
            end
            this@HypertoroidalFourierDistribution(C, transformation);
        end
        
        function r = integral(this, l, r)
            % Calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (2 x 1 column vector)
            %       left bound of integral in each dimension, default 0
            %   r (2 x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi
            % Returns:
            %   result (scalar)
            %       value of the integral
            if nargin < 2
                l = [0; 0];
            end
            if nargin < 3
                r = [2 * pi; 2 * pi];
            end
            assert(all(size(l) == [this.dim, 1]));
            assert(all(size(r) == [this.dim, 1]));
            
            switch this.transformation
                case 'sqrt'
                    tfd = this.transformViaCoefficients('square', 2*size(this.C)-1);
                case 'identity'
                    tfd = this;
                otherwise
                    error('Transformation not recognized or unsupported');
            end
            maxInd = (size(tfd.C) - 1) / 2;
            jRanges = -maxInd(1):maxInd(1);
            kRanges = -maxInd(2):maxInd(2);
            % Calculate separate integrals
            jVector = -1i ./ jRanges .* (exp(1i*jRanges*r(1)) - exp(1i*jRanges*l(1)));
            jVector(maxInd(1)+1) = r(1) - l(1);
            kVector = -1i ./ kRanges .* (exp(1i*kRanges*r(2)) - exp(1i*kRanges*l(2)));
            kVector(maxInd(2)+1) = r(2) - l(2);
            % Create all combinations to and multiply factors for the
            % coefficients
            r = real(sum(tfd.C.*(jVector.' * kVector),[1,2]));
        end
        
        function Cov4D = covariance4D(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), cos(x2), sin(x2)]
            %
            % Returns:
            %   Cov4D (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), cos(x2), sin(x2)]
            Cov4D = this.covariance2dimD;
        end
        
        function rhoc = circularCorrelationJammalamadaka(this)
            % Calculates Jammalamadaka's correlation coefficient from
            % Fourier coefficients
            switch this.transformation
                case 'sqrt'
                    tfd = this.transformViaCoefficients('square', 2*size(this.C)-1);
                case 'identity'
                    tfd = this;
                otherwise
                    error('Transformation not recognized or unsupported');
            end
            m = tfd.circularMean();
            tfd = tfd.truncate([5, 5]); % Prevent index out of bounds
            % c_{0,0} is assumed to be the entry for which both exponents
            % are zero. Calculate indices for this entry
            j = (size(tfd.C, 1) + 1) / 2;
            k = (size(tfd.C, 2) + 1) / 2;
            % Only c_{1,1}, c_{-1,-1}, c{1,-1}, and
            % c{-1,1} yield non-zero terms for the integration
            EsinAsinB = -tfd.C(j+1, k+1) * pi^2 * exp(1i*(m(1) + m(2))) - tfd.C(j-1, k-1) * pi^2 * exp(-1i*(m(1) + m(2))) ...
                +tfd.C(j+1, k-1) * pi^2 * exp(1i*(m(1) - m(2))) + tfd.C(j-1, k+1) * pi^2 * exp(-1i*(m(1) - m(2)));
            % Only c_{0,0}, c{2,0}, and c{-2,0} yield non-zero terms
            EsinAsquared = 0.5 - tfd.C(j+2, k) * exp(2i*m(1)) * pi^2 - tfd.C(j-2, k) * exp(-2i*m(1)) * pi^2;
            % Only c_{0,0}, c{0,2}, and c{0,-2} yield non-zero terms
            EsinBsquared = 0.5 - tfd.C(j, k+2) * exp(2i*m(2)) * pi^2 - tfd.C(j, k-2) * exp(-2i*m(2)) * pi^2;
            rhoc = real(EsinAsinB) / sqrt(real(EsinAsquared)*real(EsinBsquared));
        end
        
        function twn = toTWN(this)
            if strcmp(this.transformation, 'identity')
                tfd = this;
            elseif strcmp(this.transformation, 'sqrt')
                tfd = this.transformViaCoefficients('square');
            else
                error('Conversion to TWN is currently not supported for this transformation');
            end
            % Convert to twn. This is works for TWNs but currently not for
            % other distributions
            tfdMeanZero = tfd.shift(-tfd.circularMean);
            tfdTrunc = tfdMeanZero.truncate([3, 3]);
            Cov = NaN(2);
            % The ones on the diagonal are easy because we can use the
            % characteristic function formula and solve for specific terms
            % via using C_1,0=(1/2*pi)^dim*exp(-0.5*[1,0]*Cov*[1;0]) and
            % for C_0,1 respectively.
            Cov(1, 1) = -2 * tfdTrunc.dim * log(2*pi) - 2 * log(tfdTrunc.C(3, 2));
            Cov(2, 2) = -2 * tfdTrunc.dim * log(2*pi) - 2 * log(tfdTrunc.C(2, 3));
            % By calculating C_1,1 / C_1,-1 we get
            % exp(-0.5*(2*Cov(1,2)+2*Cov(2,1)) and using that the off
            % diagonal entries are equal, derive the actual value.
            Cov(1, 2) = 0.5 * log(tfdTrunc.C(1, 3)) - 0.5 * log(tfdTrunc.C(3, 3));
            Cov(2, 1) = Cov(1, 2);
            % Remove imaginary parts due to imprecision
            assert(abs(sum(imag(Cov(:)))) < 0.01);
            Cov = real(Cov);
            twn = ToroidalWNDistribution(tfd.circularMean, Cov);
        end
    end
    
    methods(Static)
        function f = fromFunction(fun, noOfCoefficients, desiredTransformation)
            if nargin == 2, desiredTransformation = 'sqrt'; end
            if length(noOfCoefficients) == 1, noOfCoefficients = noOfCoefficients * ones(1, 2);
            end
            fhyper = fromFunction@HypertoroidalFourierDistribution(fun, noOfCoefficients, desiredTransformation);
            f = ToroidalFourierDistribution(fhyper.C, fhyper.transformation);
        end
        function f = fromFunctionValues(fvals, noOfCoefficients, desiredTransformation)
            if nargin == 2, desiredTransformation = 'sqrt';
            end
            if length(noOfCoefficients) == 1, noOfCoefficients = noOfCoefficients * ones(1, 2);
            end
            fhyper = fromFunctionValues@HypertoroidalFourierDistribution(fvals, noOfCoefficients, desiredTransformation);
            f = ToroidalFourierDistribution(fhyper.C, fhyper.transformation);
        end
        function f = fromDistribution(distribution, noOfCoefficients, desiredTransformation)
            if nargin == 2, desiredTransformation = 'sqrt';
            end
            if length(noOfCoefficients) == 1, noOfCoefficients = noOfCoefficients * ones(1, 2);
            end
            fhyper = fromDistribution@HypertoroidalFourierDistribution(distribution, noOfCoefficients, desiredTransformation);
            f = ToroidalFourierDistribution(fhyper.C, fhyper.transformation);
        end
    end
end
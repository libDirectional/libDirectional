classdef FourierDistributionComplex < HypertoroidalFourierDistribution & AbstractCircularDistribution
    % Circular variant
    %
    % Florian Pfaff, Gerhard Kurz, Uwe D. Hanebeck,
    % Multivariate Angular Filtering Using Fourier Series
    % Journal of Advances in Information Fusion, 11(2):206-226, December 2016.
    methods
        function this = FourierDistributionComplex(C, transformation)
            arguments
                C {mustBeNonempty} % Class is checked below to provide better error message
                transformation char = 'sqrt' % Square root of density is standard case
            end
            this@HypertoroidalFourierDistribution(C, transformation);
            assert(iscolumn(C));
        end
        
        
        function val = cdf(this, xa, startingPoint)
            fd = FourierDistribution.fromComplex(this.C.', this.transformation);
            val = fd.cdf(xa, startingPoint);
        end
    end
    
    methods(Static)
        function fd = fromFunction(fun, noOfCoefficients, desiredTransformation)
            hfd = HypertoroidalFourierDistribution.fromFunction(fun, noOfCoefficients, desiredTransformation);
            fd = hfd.toCircular;
        end
        
        function fd = fromFunctionValues(fvals, noOfCoefficients, desiredTransformation)
            hfd = HypertoroidalFourierDistribution.fromFunctionValues(fvals, noOfCoefficients, desiredTransformation);
            fd = hfd.toCircular;
        end
        
        function fd = fromDistribution(distribution, noOfCoefficients, desiredTransformation)
            hfd = HypertoroidalFourierDistribution.fromDistribution(distribution, noOfCoefficients, desiredTransformation);
            fd = hfd.toCircular;
        end
        
        function fd = fromSamples(samples, noOfCoefficients, desiredTransformation)
            hfd = HypertoroidalFourierDistribution.fromSamples(samples, noOfCoefficients, desiredTransformation);
            fd = hfd.toCircular;
        end
    end
end
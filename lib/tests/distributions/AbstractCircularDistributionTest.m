classdef AbstractCircularDistributionTest< matlab.unittest.TestCase
    properties
    end
    
    methods (Test)
        function testAbstractCircularDistribution(testCase)
            x = 0:6;
            
            for dist = {WNDistribution(2, 0.7), VMDistribution(6, 1.2)};
                current = dist{1};
                
                %% cdf numerical
                testCase.verifyEqual(current.cdfNumerical(x), current.cdf(x), 'RelTol', 1E-10);
                startingPoint = 2.1;
                testCase.verifyEqual(current.cdfNumerical(x, startingPoint), current.cdf(x, startingPoint), 'RelTol', 1E-10);
                
                %% angular moment numerical
                testCase.verifyEqual(current.trigonometricMoment(0), current.trigonometricMomentNumerical(0), 'RelTol', 1E-10);
                testCase.verifyEqual(current.trigonometricMoment(1), current.trigonometricMomentNumerical(1), 'RelTol', 1E-10);
                testCase.verifyEqual(current.trigonometricMoment(2), current.trigonometricMomentNumerical(2), 'RelTol', 1E-10);
                testCase.verifyEqual(current.trigonometricMoment(3), current.trigonometricMomentNumerical(3), 'RelTol', 1E-10);
                
                %% circular mean/variance
                m1 = current.trigonometricMoment(1);
                testCase.verifyEqual(current.circularMean, mod(atan2(imag(m1),real(m1)),2*pi), 'RelTol', 1E-10);
                testCase.verifyEqual(cos(current.circularMean), real(m1)/(1-current.circularVariance()), 'RelTol', 1E-10);
                testCase.verifyEqual(sin(current.circularMean), imag(m1)/(1-current.circularVariance()), 'RelTol', 1E-10);
                testCase.verifyEqual(1 - current.circularVariance, sqrt(real(m1)^2+imag(m1)^2), 'RelTol', 1E-10);
                
                %% integral numerical
                testCase.verifyEqual(current.integralNumerical(2), current.integral(2), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(2,3), current.integral(2,3), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(5,4), current.integral(5,4), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(0,4*pi), current.integral(0,4*pi), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(-pi,pi), current.integral(-pi,pi), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(0,4*pi), current.integral(0,4*pi), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(-3*pi,3*pi), current.integral(-3*pi,3*pi), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(-1, 20), current.integral(-1, 20), 'RelTol', 1E-10);
                testCase.verifyEqual(current.integralNumerical(12, -3), current.integral(12, -3), 'RelTol', 1E-10);
                
                %% entropy numerical
                testCase.verifyEqual(current.entropyNumerical(), current.entropy(), 'RelTol', 1E-10);
                
                %% squared distance
                testCase.verifyEqual(current.squaredDistanceNumerical(current), 0, 'RelTol', 1E-10);
                testCase.verifyGreaterThan(current.squaredDistanceNumerical(WNDistribution(3,2)), 0);
                
                %% kld numerical
                testCase.verifyEqual(current.kldNumerical(current), 0, 'RelTol', 1E-10);
                testCase.verifyGreaterThan(current.kldNumerical(WNDistribution(3,2)), 0);
                
                %% test sampling
                rng default
                n = 100;
                s = current.sampleMetropolisHastings(n);
                testCase.verifyEqual(size(s,1), 1);
                testCase.verifyEqual(size(s,2), n);
                testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
                testCase.verifyLessThan(s,2*pi*ones(size(s)));
                
                rng default
                n = 10;
                s = current.sampleCdf(n);
                testCase.verifyEqual(size(s,1), 1);
                testCase.verifyEqual(size(s,2), n);
                testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
                testCase.verifyLessThan(s,2*pi*ones(size(s)));
                
                %% test periodicity
                testCase.verifyEqual(current.pdf(linspace(-2*pi,0,100)),current.pdf(linspace(0,2*pi,100)),'RelTol',1E-10);
            end
        end
    end
end
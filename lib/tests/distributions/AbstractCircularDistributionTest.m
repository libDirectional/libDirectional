classdef AbstractCircularDistributionTest< matlab.unittest.TestCase
   
    methods (Test)
        function testAbstractCircularDistribution(testCase)
            x = 0:6;
            
            for dist = {WNDistribution(2, 0.7), VMDistribution(6, 1.2)}
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
                
                %% conversions
                testCase.verifyEqual(current.trigonometricMoment(1), current.toWN().trigonometricMoment(1), 'RelTol', 1E-10);
                testCase.verifyEqual(current.trigonometricMoment(1), current.toVM().trigonometricMoment(1), 'RelTol', 1E-10);
                testCase.verifyEqual(current.trigonometricMoment(1), current.toWC().trigonometricMoment(1), 'RelTol', 1E-10);
                
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
        
        function testMetropolisHastings(testCase)
            dist = VMDistribution(1.3,1.8);
            rng default
            n = 1000;
            s = dist.sampleMetropolisHastings(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
            wd = WDDistribution(s);
            dist2 = wd.toVM();
            testCase.verifyEqual(dist.mu, dist2.mu, 'AbsTol', 0.1);
            testCase.verifyEqual(dist.kappa, dist2.kappa, 'RelTol', 0.1);
        end
        
        function testToDiracBT(testCase)
            dist = PWCDistribution([0, 1]);
            
            for n = [1:10, 20, 100]
                for fixCircularMean = [true false]
                    for fixFirstTrigonometrictMoment = [true false]
                        wd = dist.toDiracBT(n, fixCircularMean, fixFirstTrigonometrictMoment);
                        testCase.verifyClass(wd, 'WDDistribution');
                        testCase.verifyEqual(length(wd.d), n);
                        testCase.verifyEqual(length(wd.w), n);
                        if fixCircularMean
                            testCase.verifyEqual(dist.circularMean(), wd.circularMean(), 'RelTol', 1E-10);
                        end
                        if ~fixCircularMean && ~fixFirstTrigonometrictMoment && n>=2
                            testCase.verifyEqual(wd.integral(0,pi), dist.integral(0,pi), 'RelTol', 1E-10);
                            testCase.verifyEqual(wd.integral(pi,2*pi), dist.integral(pi,2*pi), 'RelTol', 1E-10);
                        end
                        if fixFirstTrigonometrictMoment && n>=2 %it is not possible to fix the trigonometric moment for n=1
                            testCase.verifyEqual(dist.trigonometricMoment(1), wd.trigonometricMoment(1), 'RelTol', 1E-3);
                        else
                            testCase.verifyGreaterThanOrEqual(wd.d, pi*ones(size(wd.d)));
                        end
                    end
                end
            end
        end
        
        function testToDirac5Superposition(testCase)
            dist = VMDistribution(1,2);
            
            % test scalar lambda
            for lambda = 1:5
                wd = dist.toDirac5SuperPosition(lambda);
                testCase.verifyEqual(dist.trigonometricMoment(1), wd.trigonometricMoment(1), 'RelTol', 1E-10);
                testCase.verifyEqual(dist.trigonometricMoment(2), wd.trigonometricMoment(2), 'RelTol', 1E-10);
                testCase.verifySize(wd.d, [1, 1 + 4*lambda]);
            end
            
            % test vector lambda
            lambda = [0.1 0.5 0.9];
            wd = dist.toDirac5SuperPosition(lambda);
            testCase.verifyEqual(dist.trigonometricMoment(1), wd.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(dist.trigonometricMoment(2), wd.trigonometricMoment(2), 'RelTol', 1E-10);
            testCase.verifySize(wd.d, [1, 1 + 4*length(lambda)]);
        end
        
        function testl2distanceCdfNumerical(testCase)
            vm1=VMDistribution(0,10);
            vm2=VMDistribution(0,9.8);
            vm3=VMDistribution(pi,10);
            
            testCase.verifyLessThan(vm1.l2distanceCdfNumerical(vm2), 0.001);
            testCase.verifyGreaterThan(vm1.l2distanceCdfNumerical(vm3), 0.5);
        end
    end
end
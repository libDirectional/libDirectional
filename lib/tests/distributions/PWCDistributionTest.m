classdef PWCDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testPWCDistribution(testCase)
            w = [1:5 10:-1:1];
            normal = 1 / (2 * pi * mean(w));
            p = PWCDistribution(w);
            
            % test pdf
            testCase.verifyEqual(p.pdf(0), 1 * normal, 'RelTol', 1E-10);
            testCase.verifyEqual(p.pdf(4.2), 5 * normal, 'RelTol', 1E-10);
            testCase.verifyEqual(p.pdf(10.9), 4 * normal, 'RelTol', 1E-10);
            
            % test integral
            testCase.verifyEqual(p.integral, 1, 'RelTol', 1E-6);
            testCase.verifyEqual(p.integral(0,pi)+p.integral(pi,2*pi), 1, 'RelTol', 1E-6);
            
            % test trigonometric moments
            testCase.verifyEqual(p.trigonometricMoment(1), p.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(p.trigonometricMoment(2), p.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(p.trigonometricMoment(3), p.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
            
            % test interval borders
            testCase.verifyEqual(PWCDistribution.leftBorder(1,2), 0*2*pi);
            testCase.verifyEqual(PWCDistribution.intervalCenter(1,2), 1/4*2*pi);
            testCase.verifyEqual(PWCDistribution.rightBorder(1,2), 1/2*2*pi);
            testCase.verifyEqual(PWCDistribution.leftBorder(2,2), 1/2*2*pi);
            testCase.verifyEqual(PWCDistribution.intervalCenter(2,2), 3/4*2*pi);
            testCase.verifyEqual(PWCDistribution.rightBorder(2,2), 1*2*pi);
            
            % test numerical calculation
            % we check if more samples lead to a better angular moment matching
            wn = WNDistribution(2, 1.3);
            w1 = PWCDistribution.calculateParametersNumerically(@(x) wn.pdf(x), 40);
            w2 = PWCDistribution.calculateParametersNumerically(@(x) wn.pdf(x), 45);
            w3 = PWCDistribution.calculateParametersNumerically(@(x) wn.pdf(x), 50);
            p1 = PWCDistribution(w1);
            p2 = PWCDistribution(w2);
            p3 = PWCDistribution(w3);
            delta1 = abs(wn.trigonometricMoment(1) - p1.trigonometricMoment(1));
            delta2 = abs(wn.trigonometricMoment(1) - p2.trigonometricMoment(1));
            delta3 = abs(wn.trigonometricMoment(1) - p3.trigonometricMoment(1));
            testCase.verifyLessThanOrEqual(delta2, delta1);
            testCase.verifyLessThanOrEqual(delta3, delta2);
            
            % test entropy
            testCase.verifyEqual(p.entropy(), p.entropyNumerical(), 'RelTol', 1E-6);
            
            % test sampling
            p.sample(1);
        end
    end
end
classdef CircularUniformDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testCircularUniformDistribution(testCase)
            cu = CircularUniformDistribution();
            x = [1 2 3 4 5 6];
            
            %% test pdf
            testCase.verifyEqual(cu.pdf(x), 1/(2*pi)*ones(size(x)));
            
            %% test shift
            cu2 = cu.shift(3);
            testCase.verifyEqual(cu2.pdf(x), 1/(2*pi)*ones(size(x)));
            
            %% test cdf
            testCase.verifyEqual(cu.cdf(x), cu.cdfNumerical(x), 'RelTol', 1E-10);
            testCase.verifyEqual(cu.cdf(x, 3), cu.cdfNumerical(x,3), 'RelTol', 1E-10);
            
            %% test trigonometric moments
            testCase.verifyEqual(cu.trigonometricMoment(0),cu.trigonometricMomentNumerical(0),'AbsTol', 1E-10);
            testCase.verifyEqual(cu.trigonometricMoment(0),1,'RelTol', 1E-10);
            
            testCase.verifyEqual(cu.trigonometricMoment(1),cu.trigonometricMomentNumerical(1),'AbsTol', 1E-10);
            testCase.verifyEqual(cu.trigonometricMoment(1),0,'RelTol', 1E-10);
            
            testCase.verifyEqual(cu.trigonometricMoment(2),cu.trigonometricMomentNumerical(2),'AbsTol', 1E-10);
            testCase.verifyEqual(cu.trigonometricMoment(2),0,'RelTol', 1E-10);
            
            testCase.verifyEqual(cu.trigonometricMoment(3),cu.trigonometricMomentNumerical(3),'AbsTol', 1E-10);
            testCase.verifyEqual(cu.trigonometricMoment(3),0,'RelTol', 1E-10);
            
            %% test integral
            testCase.verifyEqual(cu.integral(), cu.integralNumerical(), 'RelTol', 1E-10);
            testCase.verifyEqual(cu.integral(), 1, 'RelTol', 1E-10);
            testCase.verifyEqual(cu.integral(1,4), cu.integralNumerical(1,4), 'RelTol', 1E-10);
            testCase.verifyEqual(cu.integral(-4,11), cu.integralNumerical(-4,11), 'RelTol', 1E-10);
            testCase.verifyEqual(cu.integral(2*pi,-1), cu.integralNumerical(2*pi,-1), 'RelTol', 1E-10);
            
            %% test mean
            testCase.verifyWarning(@cu.circularMean,'MEAN:UNDEFINED');
            
            %% test entropy
            testCase.verifyEqual(cu.entropy(), cu.entropyNumerical(), 'RelTol', 1E-10);
            
            %% test sampling
            n = 10;
            s = cu.sample(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
    end
end
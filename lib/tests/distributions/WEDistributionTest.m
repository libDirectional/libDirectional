classdef WEDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testWEDistribution(testCase)
            lambda = 2;
            we = WEDistribution(lambda);
            
            pdftemp = @(x) lambda*exp(-lambda*x)/(1-exp(-2*pi*lambda));
            x = [0 1 2 3 4];
            testCase.verifyEqual(we.pdf(x), pdftemp(x), 'RelTol', 1E-10);
            
            %% test integral
            testCase.verifyEqual(we.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(we.integralNumerical, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(we.integral(0,pi)+we.integral(pi,2*pi), 1, 'RelTol', 1E-10);

            %% test angular moments
            testCase.verifyEqual(we.trigonometricMoment(1), we.trigonometricMomentNumerical(1), 'RelTol', 1E-10)
            testCase.verifyEqual(we.trigonometricMoment(2), we.trigonometricMomentNumerical(2), 'RelTol', 1E-10)
            testCase.verifyEqual(we.trigonometricMoment(3), we.trigonometricMomentNumerical(3), 'RelTol', 1E-10)
            testCase.verifyEqual(we.circularMean, atan(1/lambda));
            
            %% test entropy
            testCase.verifyEqual(we.entropy, we.entropyNumerical, 'RelTol', 1E-10)

            %% test sampling
            n = 10;
            s = we.sample(10);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        
            %% test periodicity
            testCase.verifyEqual(we.pdf(linspace(-2*pi,0,100)),we.pdf(linspace(0,2*pi,100)),'RelTol',1E-10);
            
        end
    end
end
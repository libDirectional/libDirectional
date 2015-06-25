classdef WLDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testWLDistribution(testCase)
            lambda = 2;
            kappa = 1.3;
            wl = WLDistribution(lambda, kappa);
            
            function y = laplace (x)
                % see Sreenivasa Rao Jammalamadaka, Tomasz J. Kozubowski,
                % New Families of Wrapped Distributions for Modeling Skew Circular Data, 
                % COMMUNICATIONS IN STATISTICS, Theory and Methods
                % Vol. 33, No. 9, pp. 2059�2074, 2004
                % eq (2.2)
                if x>=0
                    y = lambda/(1/kappa+kappa)*exp(-(abs(x)*lambda*kappa));
                else
                    y = lambda/(1/kappa+kappa)*exp(-(abs(x)*lambda/kappa));
                end
            end
            pdftemp = @(x) sum(arrayfun(@(z) laplace(z), x + 2*pi * (-20:20)));
            for x = [0 1 2 3 4]
                testCase.verifyEqual(wl.pdf(x), pdftemp(x), 'RelTol', 1E-10);
            end
            
            %% test integral
            testCase.verifyEqual(wl.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wl.integralNumerical, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(wl.integral(0,pi)+wl.integral(pi,2*pi), 1, 'RelTol', 1E-10);
                                   
            %% test angular moments
            testCase.verifyEqual(wl.trigonometricMoment(1), wl.trigonometricMomentNumerical(1), 'RelTol', 1E-10)
            testCase.verifyEqual(wl.trigonometricMoment(2), wl.trigonometricMomentNumerical(2), 'RelTol', 1E-10)
            testCase.verifyEqual(wl.trigonometricMoment(3), wl.trigonometricMomentNumerical(3), 'RelTol', 1E-10)
        
            %% test periodicity
            testCase.verifyEqual(wl.pdf(linspace(-2*pi,0,100)),wl.pdf(linspace(0,2*pi,100)),'RelTol',1E-10);
            
        end
    end
end
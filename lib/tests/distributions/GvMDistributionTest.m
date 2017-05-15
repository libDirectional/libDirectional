classdef GvMDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function compareWithVm(testCase)
            vm=VMDistribution(1,10);
            gvm=GvMDistribution(1,10);
            gvm2=GvMDistribution([1;2],[10;1E-15]);
            xvals=linspace(0,2*pi,100);
            testCase.verifyEqual(gvm.pdf(xvals),vm.pdf(xvals),'AbsTol',1E-10);
            testCase.verifyEqual(gvm2.pdf(xvals),vm.pdf(xvals),'AbsTol',1E-10);
        end
        
        function testPdfUnnormalized(testCase)
            gvm=GvMDistribution([3;2],[10;5]);
            xvals=linspace(0,2*pi,100);
            pdfUnnormalized = @(x) exp(cos(x-gvm.mu(1))*gvm.kappa(1) + cos(2*(x-gvm.mu(2)))*gvm.kappa(2));
            testCase.verifyEqual(gvm.pdfUnnormalized(xvals),pdfUnnormalized(xvals),'AbsTol',1E-10);
        end        
        
        function testNormalized(testCase)
            gvm=GvMDistribution([1;2],[10;20]);
            testCase.verifyEqual(gvm.integral(0,2*pi),1,'AbsTol',1E-10);
        end
    end
end

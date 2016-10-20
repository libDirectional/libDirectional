classdef GvMDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function compareWithVm(testCase)
            vm=VMDistribution(1,10);
            gvm=GvMDistribution(1,10);
            xvals=linspace(0,2*pi,100);
            testCase.verifyEqual(gvm.pdf(xvals),vm.pdf(xvals),'AbsTol',1E-10);
        end
        function testNormalized(testCase)
            gvm=GvMDistribution([1;2],[10;20]);
            testCase.verifyEqual(gvm.integral(0,2*pi),1,'AbsTol',1E-10);
        end
    end
end

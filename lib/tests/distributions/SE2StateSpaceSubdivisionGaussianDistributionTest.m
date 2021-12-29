classdef SE2StateSpaceSubdivisionGaussianDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testConstructor(testCase)
            n = 100;
            vm = VMDistribution(0,1);
            testCase.verifyWarningFree(@()SE2StateSpaceSubdivisionGaussianDistribution(...
                FIGDistribution.fromDistribution(vm,n),...
                repmat(GaussianDistribution([0;0],eye(2)),[1,n])));
        end
        
        function testPlotting(testCase)
            n = 100;
            vm = VMDistribution(0,1);
            apd = SE2StateSpaceSubdivisionGaussianDistribution(...
                FIGDistribution.fromDistribution(vm,n),...
                repmat(GaussianDistribution([0;0],eye(2)),[1,n]));
            
            testCase.verifyWarningFree(@()apd.plotState());
        end
    end
        
end
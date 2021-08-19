classdef GSSVMDistributionTest < matlab.unittest.TestCase
    
    methods(Test)
        function testConstructor(testCase)
            testCase.verifyWarningFree(@()GSSVMDistribution(2,0.3,2));
        end
    end
end

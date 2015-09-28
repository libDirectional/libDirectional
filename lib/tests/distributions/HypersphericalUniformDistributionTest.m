classdef HypersphericalUniformDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testIntegral(testCase)
            for dim=1:4
                hud=HypersphericalUniformDistribution(dim);
                testCase.verifyEqual(hud.integral,1,'AbsTol',1E-6);
            end
        end
    end
end

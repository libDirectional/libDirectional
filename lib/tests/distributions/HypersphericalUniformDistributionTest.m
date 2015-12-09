classdef HypersphericalUniformDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testIntegral(testCase)
            for dim=1:4
                hud=HypersphericalUniformDistribution(dim);
                testCase.verifyEqual(hud.integral,1,'AbsTol',1E-6);
            end
        end
        
        function testPdf(testCase)
            rng default
            for dim=1:4
                hud=HypersphericalUniformDistribution(dim);
                x = rand(dim,1);
                x = x/norm(x);
                testCase.verifyEqual(hud.pdf(x),1/AbstractHypersphericalDistribution.computeUnitSphereSurface(dim),'AbsTol',1E-10);
            end
        end
    end
end

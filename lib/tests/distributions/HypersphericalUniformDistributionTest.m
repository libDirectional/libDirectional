classdef HypersphericalUniformDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testIntegral(testCase)
            for dim = 2:4
                hud = HypersphericalUniformDistribution(dim);
                testCase.verifyEqual(hud.integral,1,'AbsTol',1E-6);
            end
        end
        
        function testPdf(testCase)
            rng default
            for dim = 2:4
                hud = HypersphericalUniformDistribution(dim);
                x = rand(dim,1);
                x = x/norm(x);
                testCase.verifyEqual(hud.pdf(x),1/AbstractHypersphericalDistribution.computeUnitSphereSurface(dim),'AbsTol',1E-10);
            end
        end
        
        function testSample(testCase)
            for dim = 2:4
                hud = HypersphericalUniformDistribution(dim);
                n = 10;
                samples = hud.sample(n);
                testCase.verifySize(samples, [hud.dim, n]);
                testCase.verifyEqual(sum(samples.*samples), ones(1,n), 'RelTol', 1E-10);
            end
        end
    end
end

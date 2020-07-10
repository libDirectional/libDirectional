classdef AbstractHyperhemisphericalDistributionTest< matlab.unittest.TestCase    
    methods (Test)
        function testGetManifoldArea(testCase)
            % Circle
            hud = HyperhemisphericalUniformDistribution(2);
            testCase.verifyEqual(hud.getManifoldSize, pi, 'AbsTol', 1e-16);
            % Sphere
            hud = HyperhemisphericalUniformDistribution(3);
            testCase.verifyEqual(hud.getManifoldSize, 2*pi, 'AbsTol', 1e-16);
        end
    end
end
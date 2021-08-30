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
        
        function testSampleMetropolisHastingsBasicsOnly(testCase)
            vmf = VMFDistribution([1 0]', 2);
            chd = CustomHyperhemisphericalDistribution(@(x)vmf.pdf(x)+vmf.pdf(-x),2);
            n = 10;
            s = chd.sampleMetropolisHastings(n);
            testCase.verifySize(s, [chd.dim, n]);
            testCase.verifyEqual(sum(s.^2,1), ones(1,n), 'RelTol', 1E-10);
            
            s2 = chd.sample(n);
            testCase.verifySize(s2, [chd.dim, n]);
            testCase.verifyEqual(sum(s2.^2,1), ones(1,n), 'RelTol', 1E-10);
        end
    end
end
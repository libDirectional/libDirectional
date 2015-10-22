classdef AbstractHypersphericalDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testAbstractHypersphericalDistribution(testCase)           
            %% integral 2d
            alpha = 0.3;
            mu = [cos(alpha), sin(alpha)]';
            kappa = 1.2;
            vmf = VMFDistribution(mu,kappa);
            testCase.verifyEqual(vmf.integral, 1, 'RelTol', 1E-10);

            %% integral 3d
            mu = [cos(alpha), sin(alpha), 0]';
            vmf = VMFDistribution(mu,kappa);
            testCase.verifyEqual(vmf.integral, 1, 'RelTol', 1E-8);

            %% integral 4d
            mu = [cos(alpha), sin(alpha), 0, 0]';
            vmf = VMFDistribution(mu,kappa);
            testCase.verifyEqual(vmf.integral, 1, 'RelTol', 1E-7);
           
            %% integral 5d (monte carlo, not very precise)
            rng default
            mu = [cos(alpha), sin(alpha), 0, 0, 0]';
            vmf = VMFDistribution(mu,kappa);
            testCase.verifyEqual(vmf.integral, 1, 'RelTol', 1E-2);

            %% unit sphere surface area
            testCase.verifyEqual(AbstractHypersphericalDistribution.computeUnitSphereSurface(2), 2*pi, 'RelTol', 1E-10);
            testCase.verifyEqual(AbstractHypersphericalDistribution.computeUnitSphereSurface(3), 4*pi, 'RelTol', 1E-10);
            testCase.verifyEqual(AbstractHypersphericalDistribution.computeUnitSphereSurface(4), 2*pi^2, 'RelTol', 1E-10);
        end
        
        function testEntropyNumerical(testCase)
            %% 2D
            vmf = VMFDistribution([1 0]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.toVM().entropy(), 'RelTol', 1E-10);
            vmf = VMFDistribution([0 1]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.toVM().entropy(), 'RelTol', 1E-10);
            
            %% 3D
            vmf = VMFDistribution([0 0 1]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-7);
            vmf = VMFDistribution([0 1 0]', 0.4);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-7);
            
            %% 4D
            vmf = VMFDistribution([0 0 0 1]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-7);
        end
    end
end
classdef AbstractHypersphericalDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testIntegral(testCase)           
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
        end

        function testUnitSphereSurface(testCase)
            % unit sphere surface area
            testCase.verifyEqual(AbstractHypersphericalDistribution.computeUnitSphereSurface(2), 2*pi, 'RelTol', 1E-10);
            testCase.verifyEqual(AbstractHypersphericalDistribution.computeUnitSphereSurface(3), 4*pi, 'RelTol', 1E-10);
            testCase.verifyEqual(AbstractHypersphericalDistribution.computeUnitSphereSurface(4), 2*pi^2, 'RelTol', 1E-10);
        end
        
        function testEntropyNumerical(testCase)
            % 2D
            vmf = VMFDistribution([1 0]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.toVM().entropy(), 'RelTol', 1E-10);
            vmf = VMFDistribution([0 1]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.toVM().entropy(), 'RelTol', 1E-10);
            
            % 3D
            vmf = VMFDistribution([0 0 1]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-7);
            vmf = VMFDistribution([0 1 0]', 0.4);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-7);
            
            % 4D
            vmf = VMFDistribution([0 0 0 1]', 2);
            testCase.verifyEqual(vmf.entropyNumerical(), vmf.entropy(), 'RelTol', 1E-7);
        end
        
        function testSampleMetropolisHastingsBasicsOnly(testCase)
            vmf = VMFDistribution([1 0]', 2);
            n = 10;
            s = vmf.sampleMetropolisHastings(n);
            testCase.verifySize(s, [vmf.dim, n]);
            testCase.verifyEqual(sum(s.^2,1), ones(1,n), 'RelTol', 1E-10);
            
            s2 = vmf.sample(n);
            testCase.verifySize(s2, [vmf.dim, n]);
            testCase.verifyEqual(sum(s2.^2,1), ones(1,n), 'RelTol', 1E-10);
        end
        
        function testSampleDeterministic(testCase)
            % define VMFDistribution
            mu = rand(3,1);
            mu = mu / vecnorm(mu);
            kappa = 3*randn(1)^2 + 1;        
            distr = VMFDistribution(mu, kappa);
            % obtain deterministic samples
            n = randi(100);
            s = distr.sampleDeterministicLCD(n);
            % check size & norm of vectors
            testCase.verifySize(s, [3, n]);
            testCase.verifyEqual(vecnorm(s,2,1), ones(1,n), 'RelTol',1E-10)
            % check mean
            dd = HypersphericalDiracDistribution(s);
            testCase.verifyEqual(dd.meanDirection(), mu, 'AbsTol', 1E-1)
        end
        
        function testMeanDirectionNumerical(testCase)
            mu = 1/sqrt(2)*[1;1;0];
            vmf = VMFDistribution(mu,1);
            testCase.verifyEqual(vmf.meanDirectionNumerical, mu, 'AbsTol', 1e-6);
        end
        
        function testDistances(testCase)
            for dim=2:4
                dist1=VMFDistribution([1;0;zeros(dim-2,1)],2);
                testCase.verifyEqual(dist1.hellingerDistanceNumerical(dist1), 0, 'AbsTol', 1E-5);
                testCase.verifyEqual(dist1.totalVariationDistanceNumerical(dist1), 0, 'AbsTol', 1E-5);
                
                dist2=VMFDistribution([0;1;zeros(dim-2,1)],2);
                testCase.verifyGreaterThan(dist1.hellingerDistanceNumerical(dist2), 0.2);
                testCase.verifyGreaterThan(dist1.totalVariationDistanceNumerical(dist2), 0.2);
            end
        end

        function testModeNumerical2D(testCase)
            M = eye(2,2);
            Z = [-3 0]';
            bd = BinghamDistribution(Z,M);
            bd.F=bd.F*bd.integralNumerical();
            testCase.verifyEqual(bd.momentNumerical(), bd.moment(),'AbsTol',0.001);

            phi = 0.7;
            M = [cos(phi),-sin(phi);sin(phi),cos(phi)];
            Z = [-5 0]';
            bd = BinghamDistribution(Z,M);
            bd.F=bd.F*bd.integralNumerical();
            testCase.verifyEqual(bd.momentNumerical(), bd.moment(),'AbsTol',0.001);
        end

        function testModeNumerical3D(testCase)
            M = eye(4,4);
            Z = [-10 -2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F=bd.F*bd.integralNumerical();
            testCase.verifyEqual(bd.momentNumerical(), bd.moment(),'AbsTol',0.003);
        end

        function testModeNumerical4D(testCase)
            q = [1,2,3,4]';
            q = q/norm(q);
            M = [quaternionMultiplication(q, [1 0 0 0]'), quaternionMultiplication(q, [0 1 0 0]'), quaternionMultiplication(q, [0 0 1 0]'), quaternionMultiplication(q, [0 0 0 1]')];
            Z = [-10 -2 -1 0]';
            bd = BinghamDistribution(Z,M);
            % Improve normalization constant for Bingham distribution
            bd.F=bd.F*bd.integralNumerical();

            testCase.verifyEqual(bd.momentNumerical(), bd.moment(),'AbsTol',0.001);
        end
    end
end
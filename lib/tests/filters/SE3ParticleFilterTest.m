classdef SE3ParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function testInitialization(testCase)
            testCase.verifyWarningFree(@()SE3ParticleFilter(10));
        end
        function testSetState(testCase)
            rng default
            n = 5000;
            cpsd = SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;3],diag([3,2,1]))});
            ddist = SE3DiracDistribution.fromDistribution(cpsd,n);
            
            lpf = SE3ParticleFilter(10);
            testCase.verifyWarningFree(@()lpf.setState(ddist));
        end

        function testgetPointEstimate(testCase)
            rng default
            n = 10000;
            cpsd = SE3CartProdStackedDistribution({HyperhemisphericalWatsonDistribution([0;0;0;1],1),...
                GaussianDistribution([1;2;3],diag([3,2,1]))});
            ddist = SE3DiracDistribution.fromDistribution(cpsd,n);
            
            lpf = SE3ParticleFilter(10);
            lpf.setState(ddist);

            testCase.verifyEqual(lpf.getPointEstimate(),cpsd.mode(),'AbsTol',0.05);
        end
        
        function testPredictUpdateCycle3D(testCase)
            % With prediction update cycles, we should converge to a forced
            % (noise-free measurement)
            rng default
            n = 3000;
            cpsdPrior = SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;3],diag([3,2,1]))});
            ddistPrior = SE3DiracDistribution.fromDistribution(cpsdPrior,n);
            
            hpf=SE3ParticleFilter(n);
            hpf.setState(ddistPrior);        
            sysNoise = SE3CartProdStackedDistribution({HyperhemisphericalWatsonDistribution([0;0;0;1],1),...
                    GaussianDistribution([0;0;0],diag([3,2,1]))});
            measNoise = SE3CartProdStackedDistribution({HyperhemisphericalWatsonDistribution([0;0;0;1],1),...
                    GaussianDistribution([0;0;0],diag([3,2,1]))});
            
            forcedBound=[1;2;3;4]/norm([1;2;3;4]);
            forcedLin=[-10;0;10];
            forcedMean = [forcedBound;forcedLin];
            for i=1:10
                hpf.predictIdentity(sysNoise);
                testCase.verifySize(hpf.getPointEstimate, [7,1]);
                hpf.updateIdentity(measNoise,forcedMean);
                hpf.updateIdentity(measNoise,forcedMean);
                hpf.updateIdentity(measNoise,forcedMean);
            end
            testCase.verifySize(hpf.getPointEstimate, [7,1]);
            testCase.verifyEqual(hpf.getPointEstimate,forcedMean,'AbsTol',0.1);
        end
        
        function testUpdate3DForcedParticlePosNoPred(testCase)
            % Without predictions, converge to closest particle
            rng default
            n = 3000;
            cpsdPrior = SE3CartProdStackedDistribution({HyperhemisphericalUniformDistribution(4),GaussianDistribution([1;2;3],diag([3,2,1]))});
            measNoise = SE3CartProdStackedDistribution({HyperhemisphericalWatsonDistribution([0;0;0;1],1),...
                    GaussianDistribution([0;0;0],diag([3,2,1]))});
            
            hpf=SE3ParticleFilter(n);
            hpf.setState(cpsdPrior)
            forcedBound=[1;2;3;4]/norm([1;2;3;4]);
            forcedLin=[-10;0;10];
            forcedMean = [forcedBound;forcedLin];
            forcedFirstParticlePos = forcedMean + [zeros(6,1);0.1];
            hpf.dist.d(:,1)=forcedFirstParticlePos;
            for i=1:5
                testCase.verifySize(hpf.getPointEstimate, [7,1]);
                hpf.updateIdentity(measNoise,forcedMean);
                hpf.updateIdentity(measNoise,forcedMean);
                hpf.updateIdentity(measNoise,forcedMean);
            end
            testCase.verifySize(hpf.getPointEstimate(), [7,1]);
            testCase.verifyEqual(hpf.getPointEstimate(),forcedFirstParticlePos,'AbsTol',1e-12);
        end
    end
end


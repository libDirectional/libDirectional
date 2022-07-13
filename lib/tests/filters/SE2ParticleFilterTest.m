classdef SE2ParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function testInitialization(testCase)
            testCase.verifyWarningFree(@()SE2ParticleFilter(10));
        end
        function testSetState(testCase)
            rng default
            n = 5000;
            hwn = HypercylindricalWNDistribution([1;2;3],diag([1,2,3]),1);
            ddist = SE2PWDDistribution.fromDistribution(hwn,n);
            
            lpf = SE2ParticleFilter(10);
            testCase.verifyWarningFree(@()lpf.setState(ddist));
            testCase.verifyEqual(lpf.getPointEstimate,hwn.mu,'AbsTol',0.05);
        end
        
        function testPredictUpdateCycle3D(testCase)
            % With prediction update cycles, we should converge to a forced
            % (noise-free measurement)
            rng default
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;1;1]+pi/2;
            boundD = 1;
            hwnd=HypercylindricalWNDistribution(mu,C, boundD);
            hpf=SE2ParticleFilter(500);
            hpf.setState(hwnd)
            forcedMean=[1;10;20];
            for i=1:50
                hpf.predictIdentity(HypercylindricalWNDistribution([0;0;0],C, boundD));
                testCase.verifySize(hpf.getPointEstimate, [3,1]);
                hpf.updateIdentity(HypercylindricalWNDistribution([0;0;0],0.5*C, boundD),forcedMean);
                hpf.updateIdentity(HypercylindricalWNDistribution([0;0;0],0.5*C, boundD),forcedMean);
                hpf.updateIdentity(HypercylindricalWNDistribution([0;0;0],0.5*C, boundD),forcedMean);
            end
            testCase.verifySize(hpf.getPointEstimate, [3,1]);
            testCase.verifyEqual(hpf.getPointEstimate,forcedMean,'AbsTol',0.1);
            
            n=5;
            samples = rand(3,n);
            weights = ones(1,n)/n;
            f = @(x,w) [mod(x(1:boundD)+w(1:boundD),2*pi);mod(x(boundD+1:end)+w(boundD+1:end),2*pi)];
            hpf.setState(hwnd)
            hpf.predictNonlinearNonAdditive(f, samples, weights)
            est = hpf.getPointEstimate();
            testCase.verifySize(hpf.getPointEstimate, [3,1]);
            testCase.verifyEqual(est, mod(hwnd.mu + mean(samples,2), 2*pi), 'RelTol', 0.1);
        end
        
        function testUpdate3DForcedPartilePosNoPred(testCase)
            % Without predictions, converge to closest particle
            rng default
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;1;1]+pi/2;
            boundD=1;
            hwnd=HypercylindricalWNDistribution(mu,C, boundD);
            hpf=SE2ParticleFilter(500);
            hpf.setState(hwnd)
            forcedMean=[1;10;20];
            forceFirstParticlePos = [1.1;10;20];
            hpf.dist.d(:,1)=forceFirstParticlePos;
            for i=1:50
                testCase.verifySize(hpf.getPointEstimate, [3,1]);
                hpf.updateIdentity(HypercylindricalWNDistribution([0;0;0],0.5*C, boundD),forcedMean);
                hpf.updateIdentity(HypercylindricalWNDistribution([0;0;0],0.5*C, boundD),forcedMean);
                hpf.updateIdentity(HypercylindricalWNDistribution([0;0;0],0.5*C, boundD),forcedMean);
            end
            testCase.verifySize(hpf.getPointEstimate, [3,1]);
            testCase.verifyEqual(hpf.getPointEstimate,forceFirstParticlePos,'AbsTol',1e-12);
        end
    end
end


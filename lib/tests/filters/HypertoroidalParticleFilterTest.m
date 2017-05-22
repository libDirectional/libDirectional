classdef HypertoroidalParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function test3D(testCase)
            rng default
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;1;1]+pi/2;
            hwnd=HypertoroidalWNDistribution(mu,C);
            hpf=HypertoroidalParticleFilter(500,3);
            hpf.setState(hwnd)
            forcedMean=[1;2;3];
            for i=1:50
                hpf.predictIdentity(HypertoroidalWNDistribution([0;0;0],C));
                testCase.verifySize(hpf.getEstimateMean, [3,1]);
                hpf.updateIdentity(HypertoroidalWNDistribution([0;0;0],0.5*C),forcedMean);
                hpf.updateIdentity(HypertoroidalWNDistribution([0;0;0],0.5*C),forcedMean);
                hpf.updateIdentity(HypertoroidalWNDistribution([0;0;0],0.5*C),forcedMean);
            end
            testCase.verifySize(hpf.getEstimateMean, [3,1]);
            testCase.verifyEqual(hpf.getEstimateMean,forcedMean,'AbsTol',0.1);
            
            n=5;
            samples = rand(3,n);
            weights = ones(1,n)/n;
            f = @(x,w) mod(x+w,2*pi);
            hpf.setState(hwnd)
            hpf.predictNonlinearNonAdditive(f, samples, weights)
            est = hpf.getEstimateMean();
            testCase.verifySize(hpf.getEstimateMean, [3,1]);
            testCase.verifyEqual(est, mod(hwnd.mu + mean(samples,2), 2*pi), 'RelTol', 0.1);
        end
    end
end


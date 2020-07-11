classdef ToroidalParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function testToroidalParticleFilter(testCase)
            rng default
            C=[0.7,0.4;0.4,0.6];
            mu=[1;1]+pi/2;
            hwnd=ToroidalWNDistribution(mu,C);
            tpf=ToroidalParticleFilter(200);
            tpf.setState(hwnd)
            forcedMean=[1;1];
            for i=1:50
                tpf.predictIdentity(HypertoroidalWNDistribution([0;0],C));
                tpf.updateIdentity(HypertoroidalWNDistribution([0;0],0.5*C),forcedMean);
                tpf.updateIdentity(HypertoroidalWNDistribution([0;0],0.5*C),forcedMean);
                tpf.updateIdentity(HypertoroidalWNDistribution([0;0],0.5*C),forcedMean);
            end
            testCase.verifyEqual(tpf.getEstimateMean,forcedMean,'AbsTol',0.1);
        end
    end
end


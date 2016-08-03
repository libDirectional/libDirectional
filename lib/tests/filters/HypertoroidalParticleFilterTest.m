classdef HypertoroidalParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function test3D(testCase)
            C=[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            mu=[1;1;1]+pi/2;
            hwnd=HypertoroidalWNDistribution(mu,C);
            hpf=HypertoroidalParticleFilter(500,3);
            hpf.setState(hwnd)
            forcedMean=[1;1;1];
            for i=1:50
                hpf.predictIdentity(HypertoroidalWNDistribution([0;0;0],C));
                hpf.updateIdentity(HypertoroidalWNDistribution([0;0;0],0.5*C),forcedMean);
                hpf.updateIdentity(HypertoroidalWNDistribution([0;0;0],0.5*C),forcedMean);
                hpf.updateIdentity(HypertoroidalWNDistribution([0;0;0],0.5*C),forcedMean);
            end
            testCase.verifyEqual(hpf.getEstimateMean,[1;1;1],'AbsTol',1);
        end
    end
end


classdef HypersphericalParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function test3D(testCase)
            nSamples=5000;
            vmfInit=VMFDistribution([1;0;0],10);
            vmfSys=VMFDistribution([0;0;1],10);
            vmfMeas=VMFDistribution([1;0;0],3);

            vmfFilter=VMFFilter;
            hpf=HypersphericalParticleFilter(nSamples,3);
            vmfFilter.setState(vmfInit);
            hpf.setState(vmfInit);

            hpf.predictIdentity(vmfSys);
            vmfFilter.predictIdentity(vmfSys);

            hpf.updateIdentity(vmfMeas);
            hpf.updateIdentity(vmfMeas);
            hpf.updateIdentity(vmfMeas);
            vmfFilter.updateIdentity([1;0;0],vmfMeas);
            vmfFilter.updateIdentity([1;0;0],vmfMeas);
            vmfFilter.updateIdentity([1;0;0],vmfMeas);

            testCase.verifyEqual(norm(hpf.getEstimateMean-vmfFilter.getEstimateMean),0,'AbsTol',0.05);
        end
    end
end


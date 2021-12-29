classdef HypersphericalParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function test3D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            rng(1);
            nSamples=20000;
            vmfInit=VMFDistribution([1;0;0],10);
            vmfSys=VMFDistribution([0;0;1],10);
            vmfMeas=VMFDistribution([0;0;1],3);

            vmfFilter=VMFFilter;
            hpf=HypersphericalParticleFilter(nSamples,3);
            vmfFilter.setState(vmfInit);
            hpf.setState(vmfInit);
            est = hpf.getEstimate();
            testCase.verifyClass(est, 'HypersphericalDiracDistribution');
            testCase.verifyEqual(vmfInit.moment(),sum(diag(est.w)*est.d')', 'AbsTol', 1E-2);

            hpf.predictIdentity(vmfSys);
            vmfFilter.predictIdentity(vmfSys);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('updateIdentity:muReplaced'));
            hpf.updateIdentity(vmfMeas,[0;0;1]);
            hpf.updateIdentity(vmfMeas,[0;0;1]);
            hpf.updateIdentity(vmfMeas,[0;0;1]);
            fixture.teardown();
            vmfFilter.updateIdentity(vmfMeas,[0;0;1]);
            vmfFilter.updateIdentity(vmfMeas,[0;0;1]);
            vmfFilter.updateIdentity(vmfMeas,[0;0;1]);

            testCase.verifySize(hpf.getEstimateMean, [3,1]);
            testCase.verifyEqual(norm(hpf.getEstimateMean-vmfFilter.getEstimateMean),0,'AbsTol',0.05);
            
            hpf.setState(vmfInit);
            f = @(z,x) ones(1,size(x,2));
            z = 3;
            est = hpf.getEstimateMean();
            testCase.verifyEqual(norm(est), 1, 'RelTol', 1E-10);
            hpf.updateNonlinear(f, z);
            testCase.verifyEqual(est, hpf.getEstimateMean(), 'AbsTol', 1E-2);
        end
    end
end


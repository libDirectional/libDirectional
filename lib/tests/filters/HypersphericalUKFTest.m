classdef HypersphericalUKFTest< matlab.unittest.TestCase    
    methods (Test)
        function test3D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmfInit=VMFDistribution([1;0;0],10);
            vmfSys=VMFDistribution([0;0;1],10);
            vmfMeas=VMFDistribution([0;0;1],20);

            vmfFilter=VMFFilter;
            ukf=HypersphericalUKF();
            vmfFilter.setState(vmfInit);
            ukf.setState(vmfInit);

            ukf.predictIdentity(vmfSys);
            vmfFilter.predictIdentity(vmfSys);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('updateIdentity:muReplaced'));
            ukf.updateIdentity(vmfMeas,[0;0;1]);
            ukf.updateIdentity(vmfMeas,[0;0;1]);
            ukf.updateIdentity(vmfMeas,[0;0;1]);
            fixture.teardown;
            vmfFilter.updateIdentity(vmfMeas,[0;0;1]);
            vmfFilter.updateIdentity(vmfMeas,[0;0;1]);
            vmfFilter.updateIdentity(vmfMeas,[0;0;1]);

            est = ukf.getEstimate();
            testCase.verifyClass(est, 'GaussianDistribution');
            testCase.verifyEqual(norm(est.mu), 1, 'RelTol', 1E-10);
            testCase.verifyEqual(norm(ukf.getEstimateMean-vmfFilter.getEstimateMean),0,'AbsTol',0.05);
            
            %test nonadd. noise
            ukf.setState(vmfInit);
            vmfFilter.setState(vmfInit);
            f = @(x,w) (x + w)/norm(x+w);
            noiseSamples = rand(3,5);
            noiseWeights = ones(1,5)/5;
            ukf.predictNonlinearArbitraryNoise(f, noiseSamples, noiseWeights)
            vmfFilter.predictNonlinearArbitraryNoise(f, noiseSamples, noiseWeights);
            testCase.verifyEqual(norm(ukf.getEstimateMean-vmfFilter.getEstimateMean),0,'AbsTol',0.05);
        end
    end
end


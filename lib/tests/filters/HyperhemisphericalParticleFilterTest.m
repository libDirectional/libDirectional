classdef HyperhemisphericalParticleFilterTest< matlab.unittest.TestCase    
    methods (Test)
        function test3D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            rng(1);
            nSamples=20000;

            watsonHemiInit=HyperhemisphericalWatsonDistribution([1;0;0],10);
            watsonHemiSys=HyperhemisphericalWatsonDistribution([0;0;1],10);
            watsonHemiMeas=HyperhemisphericalWatsonDistribution([0;0;1],3);

            hgf=HyperhemisphericalGridFilter(151,3);
            hpf=HyperhemisphericalParticleFilter(nSamples,3);
            hgf.setState(watsonHemiInit);
            hpf.setState(watsonHemiInit);
            est = hpf.getEstimate();
            testCase.verifyClass(est, 'HyperhemisphericalDiracDistribution');
            % Minimum of distance of point estimate to mu or -mu should be
            % close to 0
            testCase.verifyEqual(min(vecnorm(hpf.getPointEstimate()-[watsonHemiInit.mode(),-watsonHemiInit.mode()])),0, 'AbsTol', 1E-2);
            
            testCase.verifyEqual(min(vecnorm(hpf.getPointEstimate()-[hgf.getPointEstimate(),-hgf.getPointEstimate()])),0,'AbsTol',0.05);
            hgf.predictIdentity(watsonHemiSys);
            hpf.predictIdentity(watsonHemiSys);
            testCase.verifyEqual(min(vecnorm(hpf.getPointEstimate()-[hgf.getPointEstimate(),-hgf.getPointEstimate()])),0,'AbsTol',0.05);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('updateIdentity:muReplaced'));
            hgf.updateIdentity(watsonHemiMeas,[0;0;1]);
            hgf.updateIdentity(watsonHemiMeas,[0;0;1]);
            hgf.updateIdentity(watsonHemiMeas,[0;0;1]);
            hpf.updateIdentity(watsonHemiMeas,[0;0;1]);
            hpf.updateIdentity(watsonHemiMeas,[0;0;1]);
            hpf.updateIdentity(watsonHemiMeas,[0;0;1]);
            fixture.teardown();

            testCase.verifySize(hpf.getPointEstimate(), [3,1]);
            testCase.verifyEqual(min(vecnorm(hpf.getPointEstimate()-[hgf.getPointEstimate(),-hgf.getPointEstimate()])),0,'AbsTol',0.1);
            
            hpf.setState(watsonHemiInit);
            f = @(z,x) ones(1,size(x,2));
            z = 3;
            est = hpf.getEstimateMean();
            testCase.verifyEqual(norm(est), 1, 'RelTol', 1E-10);
            hpf.updateNonlinear(f, z);
            testCase.verifyEqual(est, hpf.getEstimateMean(), 'AbsTol', 1E-2);
        end
    end
end


classdef CircularParticleFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testCircularParticleFilter(testCase)
            nParticles = 30;
            filter = CircularParticleFilter(nParticles);
            wd = filter.getEstimate();
            testCase.verifyEqual(wd.trigonometricMoment(1), 0, 'AbsTol', 1E-10);
            
            %% sanity check
            filter.setState(wd);
            wd1 = filter.getEstimate();
            testCase.verifyClass(wd1, 'WDDistribution');
            testCase.verifyEqual(wd.d, wd.d);
            testCase.verifyEqual(wd.w, wd.w);
            
            filter.setState(VMDistribution(0,1));
            testCase.verifyNumElements(filter.wd.d,nParticles);
            
            %% test sampling
            % check wether only valid dirac positions are sampled
            positions = (0:.1:1);
            wd3 = WDDistribution(positions);
            RandStream.setGlobalStream(RandStream.create('mt19937ar'));
            numSamples = 20;
            samples = wd3.sample(numSamples);
            testCase.verifyEqual(size(samples), [1 numSamples]);
            for i=1:numSamples
                testCase.verifyTrue(ismember(samples(i), positions));
            end
            
            %% test prediciton
            filter.setState(wd);
            f = @(x) x;
            wn = WNDistribution(1.3, 0.8);
            filter.predictNonlinear(f, wn);
            wd2 = filter.getEstimate();
            testCase.verifyClass(wd2, 'WDDistribution');
            
            filter.setState(wd);
            filter.predictIdentity(wn);
            wd2identity = filter.getEstimate();
            testCase.verifyClass(wd2identity, 'WDDistribution');
            testCase.verifyEqual(wd2.w, wd2identity.w);
            
            %% nonlinear test without noise
            filter.setState(wd);
            f = @(x) x^2;
            noNoise = WDDistribution((0));
            filter.predictNonlinear(f, noNoise);
            predicted = filter.getEstimate();
            testCase.verifyClass(predicted, 'WDDistribution');
            wdF = wd.applyFunction(f);
            testCase.verifyEqual(predicted, wdF, 'RelTol', 1E-10);
            
            %% prediction with non-additive noise
            filter.setState(wd);
            f = @(x,w) x + norm(w);
            wdNoise = wn.toDirac3();
            filter.predictNonlinearNonAdditive(f, wdNoise.d, wdNoise.w)
            predicted = filter.getEstimate();
            testCase.verifyClass(predicted, 'WDDistribution');
           
            %% test update
            rng default
            filter.setState(wd);
            h = @(x) x;
            z = 0;
            filter.updateNonlinear(LikelihoodFactory.additiveNoiseLikelihood(h, wn),z);
            wd3a = filter.getEstimate();
            testCase.verifyClass(wd3a, 'WDDistribution');
            
            filter.setState(wd);
            filter.updateIdentity(wn, z);
            wd3b = filter.getEstimate();
            testCase.verifyClass(wd3b, 'WDDistribution');
            
            filter.setState(WDDistribution(0:0.1:1));
            likelihood = @(z, x) x == 0.5;
            filter.updateNonlinear(likelihood, 42);
            estimation = filter.getEstimate();
            testCase.verifyClass(estimation, 'WDDistribution');
            for i=1:length(estimation.d)
                testCase.verifyEqual(estimation.d(i), 0.5);
            end
            
            %% test update with single parameter likelihood
            rng default
            filter = CircularParticleFilter(nParticles);
            filter.setState(wd);
            filter.updateNonlinear(@(x) wn.pdf(-x));
            wd3c = filter.getEstimate();
            testCase.verifyClass(wd3c, 'WDDistribution');
            testCase.verifyEqual(wd3a.d, wd3c.d);
            testCase.verifyEqual(wd3a.w, wd3c.w);
        end
        function testAssociationLikelihood(testCase)
            wd=WDDistribution([1,2,3],[1/3,1/3,1/3]);
            pf=CircularParticleFilter(3);
            pf.setState(wd);

            testCase.verifyEqual(pf.associationLikelihood(CircularUniformDistribution),...
                1/(2*pi), 'AbsTol', 1E-10);
            testCase.verifyGreaterThan(pf.associationLikelihood(VMDistribution(2,1)),1/(2*pi));
        end
    end
end

classdef CircularParticleFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testCircularParticleFilter(testCase)
            nParticles = 30;
            filter = CircularParticleFilter(nParticles);
            dist = filter.getEstimate();
            testCase.verifyEqual(dist.trigonometricMoment(1), 0, 'AbsTol', 1E-10);
            
            %% sanity check
            filter.setState(dist);
            dist1 = filter.getEstimate();
            testCase.verifyClass(dist1, 'WDDistribution');
            testCase.verifyEqual(dist.d, dist.d);
            testCase.verifyEqual(dist.w, dist.w);
            
            filter.setState(VMDistribution(0,1));
            testCase.verifyNumElements(filter.dist.d,nParticles);
            
            %% test sampling
            % check wether only valid dirac positions are sampled
            positions = (0:.1:1);
            dist3 = WDDistribution(positions);
            RandStream.setGlobalStream(RandStream.create('mt19937ar'));
            numSamples = 20;
            samples = dist3.sample(numSamples);
            testCase.verifyEqual(size(samples), [1 numSamples]);
            for i=1:numSamples
                testCase.verifyTrue(ismember(samples(i), positions));
            end
            
            %% test prediciton
            filter.setState(dist);
            f = @(x) x;
            wn = WNDistribution(1.3, 0.8);
            filter.predictNonlinear(f, wn);
            dist2 = filter.getEstimate();
            testCase.verifyClass(dist2, 'WDDistribution');
            
            filter.setState(dist);
            filter.predictIdentity(wn);
            dist2identity = filter.getEstimate();
            testCase.verifyClass(dist2identity, 'WDDistribution');
            testCase.verifyEqual(dist2.w, dist2identity.w);
            
            %% nonlinear test without noise
            filter.setState(dist);
            f = @(x) x.^2;
            noNoise = WDDistribution((0));
            filter.predictNonlinear(f, noNoise);
            predicted = filter.getEstimate();
            testCase.verifyClass(predicted, 'WDDistribution');
            distF = dist.applyFunction(f);
            testCase.verifyEqual(predicted, distF, 'RelTol', 1E-10);
            
            %% prediction with non-additive noise
            filter.setState(dist);
            f = @(x,w) x + norm(w);
            distNoise = wn.toDirac3();
            filter.predictNonlinearNonAdditive(f, distNoise.d, distNoise.w)
            predicted = filter.getEstimate();
            testCase.verifyClass(predicted, 'WDDistribution');
           
            %% test update
            rng default
            filter.setState(dist);
            h = @(x) x;
            z = 0;
            filter.updateNonlinear(LikelihoodFactory.additiveNoiseLikelihood(h, wn),z);
            dist3a = filter.getEstimate();
            testCase.verifyClass(dist3a, 'WDDistribution');
            
            filter.setState(dist);
            filter.updateIdentity(wn, z);
            dist3b = filter.getEstimate();
            testCase.verifyClass(dist3b, 'WDDistribution');
            
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
            filter.setState(dist);
            filter.updateNonlinear(@(x) wn.pdf(-x));
            dist3c = filter.getEstimate();
            testCase.verifyClass(dist3c, 'WDDistribution');
            testCase.verifyEqual(dist3a.d, dist3c.d);
            testCase.verifyEqual(dist3a.w, dist3c.w);
        end
        function testAssociationLikelihood(testCase)
            dist=WDDistribution([1,2,3],[1/3,1/3,1/3]);
            pf=CircularParticleFilter(3);
            pf.setState(dist);

            testCase.verifyEqual(pf.associationLikelihood(CircularUniformDistribution),...
                1/(2*pi), 'AbsTol', 1E-10);
            testCase.verifyGreaterThan(pf.associationLikelihood(VMDistribution(2,1)),1/(2*pi));
        end
    end
end

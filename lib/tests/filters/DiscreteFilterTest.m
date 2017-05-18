classdef DiscreteFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testDiscreteFilter(testCase)
            nParticles = 30;
            filter = DiscreteFilter(nParticles);
            wd = filter.getEstimate();
            testCase.verifyEqual(wd.trigonometricMoment(1), 0, 'AbsTol', 1E-10);
            
            %% sanity check
            filter.setState(wd);
            wd1 = filter.getEstimate();
            testCase.verifyClass(wd1, 'WDDistribution');
            testCase.verifyEqual(wd.d, wd.d);
            testCase.verifyEqual(wd.w, wd.w);
            
            %% test sampling
            % check wether only valid dirac positions are sampled
            positions = (0:.1:1);
            wd3a = WDDistribution(positions);
            RandStream.setGlobalStream(RandStream.create('mt19937ar'));
            numSamples = 20;
            samples = wd3a.sample(numSamples);
            testCase.verifyEqual(size(samples), [1 numSamples]);
            for i=1:numSamples
                testCase.verifyTrue(ismember(samples(i), positions));
            end
            
            %% test conversion from continous distribution
            vm = VMDistribution(2,3);
            filter.setState(vm);
            testCase.verifyClass(filter.getEstimate(), 'WDDistribution');
            testCase.verifyEqual(filter.getEstimateMean(), vm.circularMean(), 'RelTol', 1E-8);
            testCase.verifyEqual(filter.getEstimate().trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-8);
            
            %% test prediciton
            filter.setState(wd);
            f = @(x) x;
            wn = WNDistribution(1.3, 0.8);
            filter.predictNonlinear(f, wn);
            wd2 = filter.getEstimate();
            testCase.verifyClass(wd2, 'WDDistribution');
            
            %% prediction with additive noise
            filter.setState(wd);
            f = @(x) x^2;
            wnNoise = WNDistribution(0,0.5);
            filter.predictNonlinear(f, wnNoise);
            predicted = filter.getEstimate();
            testCase.verifyClass(predicted, 'WDDistribution');
            
            % WD noise
            filter.setState(wd);
            filter.predictNonlinear(f, wnNoise.toDirac5());
            predictedWD = filter.getEstimate();
            testCase.verifyClass(predictedWD, 'WDDistribution');
            testCase.verifyEqual(predictedWD.d, predicted.d);
            testCase.verifyEqual(predictedWD.w, predicted.w, 'AbsTol', 0.02);

            %% prediction with non-additive noise
            filter.setState(wd);
            f = @(x,w) x + norm(w);
            wdNoise = wn.toDirac3();
            filter.predictNonlinearNonAdditive(f, wdNoise.d, wdNoise.w)
            predictedNonadditive = filter.getEstimate();
            testCase.verifyClass(predictedNonadditive, 'WDDistribution');
           
            %% test update
            filter.setState(wd);
            h = @(x) x;
            filter.updateNonlinear(LikelihoodFactory.additiveNoiseLikelihood(h, wn),0);
            wd3a = filter.getEstimate();
            testCase.verifyClass(wd3a, 'WDDistribution');
            
            filter.setState(WDDistribution((0:20)/21*(2*pi)));
            likelihood = @(z, x) abs(x - 2*pi/21)<1E-5;
            filter.updateNonlinear(likelihood, 42);
            estimation = filter.getEstimate();
            testCase.verifyClass(estimation, 'WDDistribution');
            testCase.verifyEqual(estimation.w(abs(estimation.d - 2*pi/21)<1E-5), 1);
            testCase.verifyEqual(estimation.w(abs(estimation.d - 2*pi/21)>1E-5), zeros(1,sum(abs(estimation.d - 2*pi/21)>1E-5)));
            
            %% test update with single parameter likelihood
            filter = DiscreteFilter(nParticles);
            filter.setState(wd);
            filter.updateNonlinear(@(x) wn.pdf(-x));
            wd3c = filter.getEstimate();
            testCase.verifyClass(wd3c, 'WDDistribution');
            testCase.verifyEqual(wd3a.d, wd3c.d);
            testCase.verifyEqual(wd3a.w, wd3c.w);            
        end
    end
end

classdef WNFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testWNFilter(testCase)
            filter = WNFilter();
            wn = WNDistribution(1.3, 0.8);
            
            %% sanity check
            filter.setState(wn);
            wn1 = filter.getEstimate();
            testCase.verifyClass(wn1, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wn1.mu);
            testCase.verifyEqual(wn.sigma, wn1.sigma);            
        end
        
        function testPrediction(testCase)
            filter = WNFilter();
            wn = WNDistribution(1.3, 0.8);
            
            %% predict identity
            filter.setState(wn);
            filter.predictIdentity(WNDistribution(0, wn.sigma));
            wnIdentity = filter.getEstimate();
            testCase.verifyClass(wnIdentity, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnIdentity.mu);
            testCase.verifyLessThan(wn.sigma, wnIdentity.sigma);
            
            %% predict nonlinear with identity as function
            filter.setState(wn);
            filter.predictNonlinear(@(x) x, WNDistribution(0, wn.sigma));
            wnNonlinIdentity = filter.getEstimate();
            testCase.verifyClass(wnNonlinIdentity, 'WNDistribution');
            testCase.verifyEqual(wnIdentity.mu, wnNonlinIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wnIdentity.sigma, wnNonlinIdentity.sigma, 'RelTol', 1E-10);
            
            %% use linear system function
            wnPrior = WNDistribution(0, 0.8);
            filter.setState(wnPrior);
            filter.predictNonlinear(@(x) 2*x + 1, wnPrior);
            wnLinear = filter.getEstimate();
            testCase.verifyEqual(wnLinear.mu, 2*wnPrior.mu + 1, 'RelTol', 1E-10);
            testCase.verifyGreaterThan(wnLinear.sigma, wnPrior.sigma);
                                   
            %% prediction with non-additive noise
            filter.setState(wn);
            f = @(x,w) x + norm(w);
            wd = WNDistribution(0, wn.sigma).toDirac3();
            filter.predictNonlinearNonAdditive(f, wd.d, wd.w)
            wnNonadditive = filter.getEstimate();
            testCase.verifyClass(wnNonadditive, 'WNDistribution');
            testCase.verifyEqual(wnIdentity.mu, wnNonadditive.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(wnIdentity.sigma, wnNonadditive.sigma, 'RelTol', 1E-10);              
        end
        
        function testUpdate(testCase)
            filter = WNFilter();
            wn = WNDistribution(1.3, 0.8);
            measNoise = WNDistribution(0, 0.9);
            
            %% update identity
            filter.setState(wn);
            filter.updateIdentity(measNoise, wn.mu);
            wnIdentity = filter.getEstimate();
            testCase.verifyClass(wnIdentity, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnIdentity.mu);
            testCase.verifyGreaterThan(wn.sigma, wnIdentity.sigma);
            
            %% update identity with different measurement
            filter.setState(wn);
            filter.updateIdentity(measNoise, wn.mu + 0.1);
            wnIdentity2 = filter.getEstimate();
            testCase.verifyClass(wnIdentity2, 'WNDistribution');
            testCase.verifyLessThan(wn.mu, wnIdentity2.mu);
            testCase.verifyGreaterThan(wn.sigma, wnIdentity2.sigma);            
            
            %% nonlinear simple
            likelihood = LikelihoodFactory.additiveNoiseLikelihood(@(x) x, measNoise);
            filter.setState(wn);
            filter.updateNonlinear(likelihood,wn.mu);
            wnNonlinIdentity = filter.getEstimate();
            testCase.verifyClass(wnNonlinIdentity, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnNonlinIdentity.mu);
            testCase.verifyEqual(wnNonlinIdentity.sigma, wnIdentity.sigma, 'RelTol', 0.05);
            
            %% nonlinear particle
            rng default
            filter.setState(wn);
            filter.updateNonlinearParticle(likelihood,wn.mu);
            wnNonlinParticleIdentity = filter.getEstimate();
            testCase.verifyClass(wnNonlinParticleIdentity, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnNonlinParticleIdentity.mu, 'RelTol', 1E-1);
            testCase.verifyEqual(wnNonlinParticleIdentity.sigma, wnIdentity.sigma, 'RelTol', 0.2);
            
            %% nonlinear progressive
            filter.setState(wn);
            filter.updateNonlinearProgressive(likelihood,wn.mu);
            wnNonlinProgressiveIdentity = filter.getEstimate();
            testCase.verifyClass(wnNonlinProgressiveIdentity, 'WNDistribution');
            testCase.verifyEqual(wn.mu, wnNonlinProgressiveIdentity.mu, 'RelTol', 1E-5);
            testCase.verifyEqual(wnNonlinProgressiveIdentity.sigma, wnIdentity.sigma, 'RelTol', 0.05);
        end
    end
end

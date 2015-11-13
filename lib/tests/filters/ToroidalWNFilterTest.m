classdef ToroidalWNFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testToroidalWNFilter(testCase)
            filter = ToroidalWNFilter();
            mu = [5 2.5]';
            C = [1.3 1.4;
                1.4 2];
            twn = ToroidalWNDistribution(mu,C);
            
            %% sanity check
            filter.setState(twn);
            twn1 = filter.getEstimate();
            testCase.verifyClass(twn1, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn.mu, twn1.mu);
            testCase.verifyEqual(twn.C, twn1.C);
                       
            %% predict identity
            filter.setState(twn);
            filter.predictIdentity(twn);
            twn2 = filter.getEstimate();
            testCase.verifyClass(twn2, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn2.mu, mod(twn.mu+twn.mu,2*pi));
            testCase.verifyEqual(twn2.C, twn.C + twn.C);
            
            %% update identity
            filter.setState(twn);
            filter.updateIdentity(ToroidalWNDistribution([0,0]',twn.C),twn.mu);
            twn6 = filter.getEstimate();
            testCase.verifyClass(twn6, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn6.mu, twn.mu, 'RelTol', 1E-9);
            testCase.verifyLessThan(twn6.C, twn.C);
            
            %% update nonlinear
            z = 0.4;
            likelihood = LikelihoodFactory.additiveNoiseLikelihood(@(x) x(1,:)+x(2,:), WNDistribution(0,1));
                        
            rng default
            filter.setState(twn);
            filter.updateNonlinearParticle(likelihood, z);
            twn9 = filter.getEstimate();
            testCase.verifyGreaterThan(twn.C, twn9.C);
        end
    end
end

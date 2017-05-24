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
            
            %% predict nonlinear
            filter.setState(twn);
            filter.predictNonlinear(@(x) x, twn);
            twn3 = filter.getEstimate();
            testCase.verifyClass(twn3, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn3.mu, mod(twn.mu+twn.mu,2*pi), 'RelTol', 1E-10);
            testCase.verifyEqual(diag(twn3.C), diag(twn.C + twn.C), 'RelTol', 1E-10); %todo check off-diagonal entries with lower tolerance
            
            %% true nonlinear test
            filter.setState(ToroidalWNDistribution([0,0]',C));
            filter.predictNonlinear(@(x) [x(1)^3,x(2)^2], ToroidalWNDistribution([0,0]',C));
            twn5 = filter.getEstimate();
            testCase.verifyClass(twn5, 'ToroidalWNDistribution');
            testCase.verifyGreaterThan(eig(twn5.C-C), [0,0]');
            
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
            
            filter.setState(twn);
            filter.updateNonlinear(likelihood, z);
            twn8 = filter.getEstimate();
            testCase.verifyGreaterThan(twn.C, twn8.C);
            
            rng default
            filter.setState(twn);
            filter.updateNonlinearParticle(likelihood, z);
            twn9 = filter.getEstimate();
            testCase.verifyGreaterThan(twn.C, twn9.C);
            
            filter.setState(twn);
            filter.updateNonlinearProgressive(likelihood, z);
            twn10 = filter.getEstimate();
            testCase.verifyGreaterThan(twn.C, twn10.C);
            testCase.verifyEqual(twn9.mu, twn10.mu, 'RelTol', 1E-1);
            testCase.verifyEqual(twn8.mu, twn10.mu, 'RelTol', 1E-1);
        end
    end
end

classdef ConstrainedUKFTest < matlab.unittest.TestCase
   
    properties
    end
    
	methods (Test)
        function testConstrainedUKF(testCase)
            filter = ConstrainedUKF();
            g = GaussianDistribution(0.5,0.7);
            rng default
            
            %% sanity check
            filter.setState(g);
            g1 = filter.getEstimate();
            testCase.verifyClass(g1, 'GaussianDistribution');
            testCase.verifyEqual(g.mu, g1.mu);
            %testCase.verifyEqual(g.C, g1.C);
            
            %% predict identity
%             filter.setState(g);
%             filter.predictIdentity(g);
%             g2 = filter.getEstimate();
%             testCase.verifyClass(g2, 'GaussianDistribution');
%             testCase.verifyEqual(g2.mu, g.mu+g.mu);
%             testCase.verifyEqual(g2.C, g.C + g.C);
            
            %% predict nonlinear
            filter.setState(g);
            filter.predictNonlinear(@(x) x, g);
            g3 = filter.getEstimate();
            testCase.verifyClass(g3, 'GaussianDistribution');
            %testCase.verifyEqual(g3.mu, g.mu+g.mu, 'RelTol', 1E-10);
            %testCase.verifyGreaterThan(g3.C, g.C);
            
            %% true nonlinear test
            g4 = GaussianDistribution(0, 0.7);
            filter.setState(g4);
            filter.predictNonlinear(@(x) x^3, g4);
            g5 = filter.getEstimate();
            %testCase.verifyEqual(g5.mu, g4.mu, 'RelTol', 1E-10);
            %testCase.verifyGreaterThan(g5.C, g4.C);
            
            %% update identity
%             filter.setState(g);
%             filter.updateIdentity(GaussianDistribution(0,g.C),g.mu);
%             g6 = filter.getEstimate();
%             testCase.verifyClass(g6, 'GaussianDistribution');
%             testCase.verifyEqual(g6.mu, g.mu);
%             testCase.verifyEqual(g6.C, g.C/2);
            
            %% update nonlinear
            g7 = GaussianDistribution(0, 0.7);
            filter.setState(g7);
            z = 0.4;
            filter.updateNonlinear(@(x) x, g7, z);
            g8 = filter.getEstimate();
            testCase.verifyGreaterThan(g8.mu, g7.mu);
            testCase.verifyLessThan(g8.mu, z);
            testCase.verifyGreaterThan(g7.C, g8.C);
        end
    end
end

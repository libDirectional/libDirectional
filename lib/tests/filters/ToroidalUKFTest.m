classdef ToroidalUKFTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testToroidalUKF(testCase)
            filter = ToroidalUKF();
            mu = [0.5,1.5]';
            C = [1.2 0.9;
                0.9 1.4];
            g = GaussianDistribution(mu,C);
            
            %% sanity check
            filter.setState(g);
            g1 = filter.getEstimate();
            testCase.verifyClass(g1, 'GaussianDistribution');
            testCase.verifyEqual(g.mu, g1.mu);
            testCase.verifyEqual(g.C, g1.C);
            testCase.verifyEqual(g.mu, filter.getEstimateMean());
        end
        
        function testPrediction(testCase)
            filter = ToroidalUKF();
            mu = [0.5,1.5]';
            C = [1.2 0.9;
                0.9 1.4];
            g = GaussianDistribution(mu,C);
            
            %% predict identity
            filter.setState(g);
            filter.predictIdentity(g);
            g2 = filter.getEstimate();
            testCase.verifyClass(g2, 'GaussianDistribution');
            testCase.verifyEqual(g2.mu, mod(g.mu+g.mu, 2*pi));
            testCase.verifyEqual(g2.C, g.C + g.C);
            
            %% predict nonlinear
            filter.setState(g);
            filter.predictNonlinear(@(x) x, g);
            g3 = filter.getEstimate();
            testCase.verifyClass(g3, 'GaussianDistribution');
            testCase.verifyEqual(g3.mu, mod(g.mu+g.mu, 2*pi));
            testCase.verifyEqual(g3.C, g.C + g.C, 'RelTol', 1E-10);
            
            filter.setState(g);
            filter.predictNonlinear(@(x) [sin(x(1)); cos(x(2))], g);
            g3 = filter.getEstimate();
            testCase.verifyClass(g3, 'GaussianDistribution');
            testCase.verifyGreaterThanOrEqual(eig(g3.C-g.C), [0,0]');
        end
        
        function testUpdate(testCase)
            filter = ToroidalUKF();
            mu = [0.5,1.5]';
            C = [1.2 0.9;
                0.9 1.4];
            g = GaussianDistribution(mu,C);
            
            %% update identity
            filter.setState(g);
            measNoise = GaussianDistribution([0,0]',g.C);
            filter.updateIdentity(measNoise,g.mu);
            gIdentitiy = filter.getEstimate();
            testCase.verifyClass(gIdentitiy, 'GaussianDistribution');
            testCase.verifyEqual(gIdentitiy.mu, g.mu);
            testCase.verifyEqual(gIdentitiy.C, g.C/2, 'RelTol', 1E-10);
            
            %noise with non-zero mean
            filter.setState(g);
            filter.updateIdentity(GaussianDistribution([2;3],g.C),g.mu+[2;3]);
            gIdentity2 = filter.getEstimate();
            testCase.verifyClass(gIdentity2, 'GaussianDistribution');
            testCase.verifyEqual(gIdentity2.mu, g.mu);
            testCase.verifyEqual(gIdentity2.C, g.C/2, 'RelTol', 1E-10);
            
            %different measurement
            filter.setState(g);
            measNoise = GaussianDistribution([0,0]',g.C);
            filter.updateIdentity(measNoise,[1,2]');
            gIdentitiy3 = filter.getEstimate();
            testCase.verifyClass(gIdentitiy3, 'GaussianDistribution');
            testCase.verifyGreaterThan(g.C, gIdentitiy3.C);
            
            %% update nonlinear (with identity as function)
            filter.setState(g);
            filter.updateNonlinear(@(x) x, measNoise, g.mu);
            gNonlinIdentity = filter.getEstimate();
            testCase.verifyEqual(gNonlinIdentity.mu, gIdentitiy.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(gNonlinIdentity.C, gIdentitiy.C, 'RelTol', 1E-10);
            
            %with periodic measurement
            filter.setState(GaussianDistribution([0.1; 2*pi-0.1], C));
            filter.updateNonlinear(@(x) x, measNoise, [2*pi-0.5; 0.5] ,true);
            gNonlinIdentityPeriodic = filter.getEstimate();
            testCase.verifyGreaterThan(gNonlinIdentityPeriodic.mu(1), 2*pi-0.5);
            testCase.verifyLessThan(gNonlinIdentityPeriodic.mu(1),  2*pi);
            testCase.verifyGreaterThan(gNonlinIdentityPeriodic.mu(2), 0);
            testCase.verifyLessThan(gNonlinIdentityPeriodic.mu(2),  0.5);
            testCase.verifyGreaterThan(C, gNonlinIdentityPeriodic.C);
        end
    end
end

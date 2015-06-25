classdef CircularUKFTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testCircularUKF(testCase)
            filter = CircularUKF();
            g = GaussianDistribution(0.5,0.7);
            
            %% sanity check
            filter.setState(g);
            g1 = filter.getEstimate();
            testCase.verifyClass(g1, 'GaussianDistribution');
            testCase.verifyEqual(g.mu, g1.mu);
            testCase.verifyEqual(g.C, g1.C);
        end
        
        function testPrediction(testCase)
            filter = CircularUKF();
            g = GaussianDistribution(0.5,0.7);
                    
            %% predict identity
            filter.setState(g);
            filter.predictIdentity(g);
            gIdentity = filter.getEstimate();
            testCase.verifyClass(gIdentity, 'GaussianDistribution');
            testCase.verifyEqual(gIdentity.mu, g.mu+g.mu);
            testCase.verifyEqual(gIdentity.C, g.C + g.C);
            
            %% predict nonlinear
            filter.setState(g);
            filter.predictNonlinear(@(x) x, g);
            gNonlin = filter.getEstimate();
            testCase.verifyClass(gNonlin, 'GaussianDistribution');
            testCase.verifyEqual(gNonlin.mu, g.mu+g.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(gNonlin.C, g.C + g.C, 'RelTol', 1E-10);
            
            %% true nonlinear prediction
            g4 = GaussianDistribution(0, 0.7);
            filter.setState(g4);
            filter.predictNonlinear(@(x) x^3, g4);
            gTrueNonln = filter.getEstimate();
            testCase.verifyEqual(gTrueNonln.mu, g4.mu, 'RelTol', 1E-10);
            testCase.verifyGreaterThan(gTrueNonln.C, g4.C);
        end
        
        function testUpdate(testCase)
            filter = CircularUKF();
            g = GaussianDistribution(0.5,0.7);
                        
            %% update identity
            filter.setState(g);
            filter.updateIdentity(GaussianDistribution(0,g.C),g.mu);
            gIdentity = filter.getEstimate();
            testCase.verifyClass(gIdentity, 'GaussianDistribution');
            testCase.verifyEqual(gIdentity.mu, g.mu);
            testCase.verifyEqual(gIdentity.C, g.C/2);
            
            %noise with non-zero mean
            filter.setState(g);
            filter.updateIdentity(GaussianDistribution(2,g.C),g.mu+2);
            gIdentity2 = filter.getEstimate();
            testCase.verifyClass(gIdentity2, 'GaussianDistribution');
            testCase.verifyEqual(gIdentity2.mu, g.mu);
            testCase.verifyEqual(gIdentity2.C, g.C/2);
            
            %different measurement
            filter.setState(g);
            filter.updateIdentity(GaussianDistribution(0,g.C),2*pi-1);
            gIdentity3 = filter.getEstimate();
            testCase.verifyClass(gIdentity3, 'GaussianDistribution');
            testCase.verifyGreaterThan(gIdentity3.mu, 2*pi-1);
            testCase.verifyLessThan(gIdentity3.mu, 2*pi);
            testCase.verifyGreaterThan(g.C, gIdentity3.C);
            
            %% update nonlinear (with identity as function)
            g7 = GaussianDistribution(0, 0.7);
            filter.setState(g7);
            z = 0.4;
            filter.updateNonlinear(@(x) x, g7, z);
            g8 = filter.getEstimate();
            testCase.verifyGreaterThan(g8.mu, g7.mu);
            testCase.verifyLessThan(g8.mu, z);
            testCase.verifyGreaterThan(g7.C, g8.C);
            
            %with periodic measurement
            g9 = GaussianDistribution(0.1, 0.7);
            filter.setState(g9);
            z = 2*pi-0.4;
            filter.updateNonlinear(@(x) x, g7, z, true);
            g10 = filter.getEstimate();
            testCase.verifyGreaterThan(g10.mu, z);
            testCase.verifyLessThan(g10.mu, 2*pi);
            testCase.verifyGreaterThan(g7.C, g8.C);
        end
    end
end

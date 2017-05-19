classdef SE2UKFTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testSE2UKF(testCase)
            muInit = [1 0 0.3 0.4]';
            muSys = [1 0 0 0]';
            muMeas = [1 0 0 0]';
            C = [1, -0.1, 3 4;
                -0.1 1 5 6; 
                3 5  1 0; 
                4 6  0 1];
            C = C*C;
            filter = SE2UKF();
            filter.setState(GaussianDistribution(muInit, C));
            
            %% sanity check         
            est = filter.getEstimate();
            testCase.verifyClass(est, 'GaussianDistribution');
            testCase.verifyEqual(est.mu, muInit);
            testCase.verifyEqual(est.C, C)
            
            %% predict identity
            filter.predictIdentity(GaussianDistribution(muSys, C));
            gaussPredictIdentity = filter.getEstimate();
            testCase.verifyClass(gaussPredictIdentity, 'GaussianDistribution');
            testCase.verifyEqual(norm(gaussPredictIdentity.mu(1:2)), 1);
            testCase.verifyEqual(gaussPredictIdentity.mu(1:2), muInit(1:2), 'AbsTol', 1E-2);
            testCase.verifyEqual(gaussPredictIdentity.mu(3:4), muInit(3:4), 'AbsTol', 1E-10);            
            
            %% update identity
            z = muInit;
            filter2 = SE2UKF();
            filter2.setState(filter.getEstimate());
            filter.updateIdentity(GaussianDistribution(muMeas, C), z);
            filter2.updateIdentity(GaussianDistribution(muMeas, C), -z); %updating with -z should give the same result as updating with z
            
            gaussUpdateIdentity = filter.getEstimate();
            testCase.verifyClass(gaussUpdateIdentity, 'GaussianDistribution');
            testCase.verifyEqual(norm(gaussUpdateIdentity.mu(1:2)), 1, 'RelTol', 1E-10);
            testCase.verifyEqual(gaussUpdateIdentity.mu(1:2), muInit(1:2), 'AbsTol', 1E-2);
            testCase.verifyEqual(gaussUpdateIdentity.mu(3:4), muInit(3:4), 'AbsTol', 1E-1);            
            
            gaussUpdateIdentity2 = filter2.getEstimate();
            testCase.verifyClass(gaussUpdateIdentity2, 'GaussianDistribution');
            testCase.verifyEqual(gaussUpdateIdentity.mu, gaussUpdateIdentity2.mu);
            testCase.verifyEqual(gaussUpdateIdentity.C, gaussUpdateIdentity2.C);
        end
    end
end

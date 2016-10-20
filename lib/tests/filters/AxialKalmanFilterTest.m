classdef AxialKalmanFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testAxialKalmanFilterFilter4D(testCase)            
            filter = AxialKalmanFilter();
            mu = [1,2,3,4]';
            mu = mu/norm(mu);
            C = 0.3*eye(4,4);
                        
            %% sanity check
            filter.setState(GaussianDistribution(mu, C));
            gauss1 = filter.getEstimate();
            testCase.verifyEqual(mu, gauss1.mu);
            testCase.verifyEqual(C, gauss1.C);
            
            %% predict identity
            %"zero-mean"
            filter.setState(GaussianDistribution(mu, C));
            filter.predictIdentity(GaussianDistribution([1,0,0,0]', 0.1*eye(4,4)));
            gauss2 = filter.getEstimate();
            testCase.verifyEqual(mu, gauss2.mu);
            testCase.verifyGreaterThanOrEqual(gauss2.C, C);
            
            %non "zero-mean"
            filter.setState(GaussianDistribution(mu, C));
            filter.predictIdentity(GaussianDistribution(mu, 0.1*eye(4,4)));
            gauss3 = filter.getEstimate();
            testCase.verifyEqual(gauss3.mu, quaternionMultiplication(mu,mu));
            testCase.verifyGreaterThanOrEqual(gauss3.C, C);
            
            %% update identity
            filter.setState(GaussianDistribution(mu, C));
            z = mu;
            filter.updateIdentity(GaussianDistribution([1;0;0;0], C), z);
            gauss4 = filter.getEstimate();
            testCase.verifyEqual(mu, gauss4.mu, 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(gauss4.C, C); 
            
            filter.setState(GaussianDistribution(mu, C));
            z = -mu; %antipodally symmetric to prevous example
            filter.updateIdentity(GaussianDistribution([1;0;0;0], C), z);
            gauss5 = filter.getEstimate();
            testCase.verifyEqual(gauss4.mu, gauss5.mu);
            testCase.verifyEqual(gauss4.C, gauss5.C); 
            
            filter.setState(GaussianDistribution(mu, C));
            z = [0,0,0,1]';
            filter.updateIdentity(GaussianDistribution(mu, C), z);
            gauss6 = filter.getEstimate();
            testCase.verifyEqual(norm(gauss6.mu), 1, 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(gauss6.C, C); 
        end
        
        function testAxialKalmanFilterFilter2D(testCase)            
            filter = AxialKalmanFilter();
            mu = [1,2]';
            mu = mu/norm(mu);
            C = 0.3*eye(2,2);
                        
            %% sanity check
            filter.setState(GaussianDistribution(mu, C));
            gauss1 = filter.getEstimate();
            testCase.verifyEqual(mu, gauss1.mu);
            testCase.verifyEqual(C, gauss1.C);
            
            %% predict identity
            %"zero-mean"
            filter.setState(GaussianDistribution(mu, C));
            filter.predictIdentity(GaussianDistribution([1,0]', 0.1*eye(2,2)));
            gauss2 = filter.getEstimate();
            testCase.verifyEqual(mu, gauss2.mu);
            testCase.verifyGreaterThanOrEqual(gauss2.C, C);
            
            %non "zero-mean"
            filter.setState(GaussianDistribution(mu, C));
            filter.predictIdentity(GaussianDistribution(mu, 0.1*eye(2,2)));
            gauss3 = filter.getEstimate();
            testCase.verifyEqual(gauss3.mu, complexMultiplication(mu,mu));
            testCase.verifyGreaterThanOrEqual(gauss3.C, C);
            
            %% update identity
            filter.setState(GaussianDistribution(mu, C));
            z = mu;
            filter.updateIdentity(GaussianDistribution([1;0;], C), z);
            gauss4 = filter.getEstimate();
            testCase.verifyEqual(mu, gauss4.mu, 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(gauss4.C, C); 
            
            filter.setState(GaussianDistribution(mu, C));
            z = -mu; %antipodally symmetric to prevous example
            filter.updateIdentity(GaussianDistribution([1;0;], C), z);
            gauss5 = filter.getEstimate();
            testCase.verifyEqual(gauss4.mu, gauss5.mu);
            testCase.verifyEqual(gauss4.C, gauss5.C); 
            
            filter.setState(GaussianDistribution(mu, C));
            z = [0,1]';
            filter.updateIdentity(GaussianDistribution(mu, C), z);
            gauss6 = filter.getEstimate();
            testCase.verifyEqual(norm(gauss6.mu), 1);
            testCase.verifyLessThanOrEqual(gauss6.C, C); 
        end
    end
end

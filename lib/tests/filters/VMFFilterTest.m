classdef VMFFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testVMFFilter2d(testCase)
            filter = VMFFilter();
            vmFilter = VMFilter();
            phi = 0.3;
            mu = [cos(phi); sin(phi)];
            kappa = 0.7;
            vmf = VMFDistribution(mu,kappa);
            
            %% sanity check
            filter.setState(vmf);
            vmf1 = filter.getEstimate();
            testCase.verifyClass(vmf1, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmf1.mu);
            testCase.verifyEqual(vmf.kappa, vmf1.kappa);
                        
            %% predict identity
            filter.setState(vmf);
            vmFilter.setState(vmf.toVM());
            vmfNoise = VMFDistribution([1;0], 0.9);
            filter.predictIdentity(vmfNoise);
            vmFilter.predictIdentity(vmfNoise.toVM());
            vmfPredictIdentity = filter.getEstimate();
            vmPredictIdentity = vmFilter.getEstimate();
            testCase.verifyClass(vmfPredictIdentity, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(vmfPredictIdentity.kappa, vmf.kappa);
            testCase.verifyEqual([cos(vmPredictIdentity.mu); sin(vmPredictIdentity.mu)], vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmPredictIdentity.kappa, vmfPredictIdentity.kappa, 'RelTol', 1E-10);
            
            %todo: non-zero mean noise
            
            %% update identity
            filter.setState(vmf);
            vmFilter.setState(vmf.toVM());
            filter.updateIdentity(vmf.mu, vmfNoise);
            vmFilter.updateIdentity(vmf.toVM().mu, vmfNoise.toVM());
            vmfUpdateIdentity = filter.getEstimate();
            vmUpdateIdentity = vmFilter.getEstimate();
            testCase.verifyClass(vmfUpdateIdentity, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfUpdateIdentity.mu, 'RelTol', 1E-10); 
            testCase.verifyGreaterThanOrEqual(vmfUpdateIdentity.kappa, vmf.kappa);
            testCase.verifyEqual([cos(vmUpdateIdentity.mu); sin(vmUpdateIdentity.mu)], vmfUpdateIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmUpdateIdentity.kappa, vmfUpdateIdentity.kappa, 'RelTol', 1E-10);
            
            return
            
            %% update identity with different measurement
            filter.setState(vmf);
            z = vmf.mode()+[0.1,0]';
            filter.updateIdentity(vmfNoise, z/norm(z));
            B4 = filter.getEstimate();
            testCase.verifyClass(B4, 'BinghamDistribution');
            testCase.verifyLessThanOrEqual(B4.Z, vmf.Z);
            
            %% predict nonlinear
            filter.setState(vmf);
            vmfNoise = BinghamDistribution([-3 0]', [0 1; 1 0]);
            filter.predictNonlinear(@(x) x, vmfNoise);
            B5 = filter.getEstimate();
            testCase.verifyClass(B5, 'BinghamDistribution');
            testCase.verifyEqual(abs(vmf.M), abs(B5.M), 'RelTol', 1E-10); %each column of M is only determined up to sign
            testCase.verifyGreaterThanOrEqual(B5.Z, vmf.Z);
        end
        

    end
end

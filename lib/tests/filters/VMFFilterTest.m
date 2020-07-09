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
            vmfNoise = VMFDistribution([0;1], 0.9);
            filter.predictIdentity(vmfNoise);
            vmDisregardMu=vmfNoise.toVM();
            vmDisregardMu.mu=0;
            vmFilter.predictIdentity(vmDisregardMu);
            vmfPredictIdentity = filter.getEstimate();
            vmPredictIdentity = vmFilter.getEstimate();
            testCase.verifyClass(vmfPredictIdentity, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(vmfPredictIdentity.kappa, vmf.kappa);
            testCase.verifyEqual([cos(vmPredictIdentity.mu); sin(vmPredictIdentity.mu)], vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmPredictIdentity.kappa, vmfPredictIdentity.kappa, 'RelTol', 1E-10);
            
            %% predict nonlinear with identity as function
            filter.setState(vmf);
            filter.predictNonlinear(@(x) x, vmfNoise);
            vmfPredictNonLinear = filter.getEstimate();
            testCase.verifyClass(vmfPredictNonLinear, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmfPredictIdentity.mu, vmfPredictNonLinear.mu, 'RelTol', 1E-10);
            
            %todo: non-zero mean noise            
            
            %% update identity
            filter.setState(vmf);
            vmFilter.setState(vmf.toVM());
            filter.updateIdentity(vmfNoise, vmf.mu);
            vmFilter.updateIdentity(vmDisregardMu, vmf.toVM().mu);
            vmfUpdateIdentity = filter.getEstimate();
            vmUpdateIdentity = vmFilter.getEstimate();
            testCase.verifyClass(vmfUpdateIdentity, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfUpdateIdentity.mu, 'RelTol', 1E-10); 
            testCase.verifyGreaterThanOrEqual(vmfUpdateIdentity.kappa, vmf.kappa);
            testCase.verifyEqual([cos(vmUpdateIdentity.mu); sin(vmUpdateIdentity.mu)], vmfUpdateIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmUpdateIdentity.kappa, vmfUpdateIdentity.kappa, 'RelTol', 1E-10);
        end
        
        function testVMFFilter3d(testCase)
            filter = VMFFilter();
            phi = 0.3;
            mu = [cos(phi); sin(phi); 0];
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
            vmfNoise = VMFDistribution([0;0;1], 0.9);
            filter.predictIdentity(vmfNoise);
            vmfPredictIdentity = filter.getEstimate();
            testCase.verifyClass(vmfPredictIdentity, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(vmfPredictIdentity.kappa, vmf.kappa);
            
            %% predict nonlinear with identity as function
            filter.setState(vmf);
            filter.predictNonlinear(@(x) x, vmfNoise);
            vmfPredictNonLinear = filter.getEstimate();
            testCase.verifyClass(vmfPredictNonLinear, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfPredictIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmfPredictIdentity.mu, vmfPredictNonLinear.mu, 'RelTol', 1E-10);
            
            %todo: non-zero mean noise            
            
            %% update identity
            filter.setState(vmf);
            filter.updateIdentity(vmfNoise, vmf.mu);
            vmfUpdateIdentity = filter.getEstimate();
            testCase.verifyClass(vmfUpdateIdentity, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, vmfUpdateIdentity.mu, 'RelTol', 1E-10); 
            testCase.verifyGreaterThanOrEqual(vmfUpdateIdentity.kappa, vmf.kappa);
        end
    end
end

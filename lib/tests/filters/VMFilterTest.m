classdef VMFilterTest < matlab.unittest.TestCase
   
    properties
    end
    
    methods (Test)
        function testVMFilter(testCase)
            filter = VMFilter();
            vm = VMDistribution(2.1,1.3);
            
            %% sanity check
            filter.setState(vm);
            vm1 = filter.getEstimate();
            testCase.verifyClass(vm1, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vm1.mu);
            testCase.verifyEqual(vm.kappa, vm1.kappa);
        end
        
        function testPrediction(testCase)
            filter = VMFilter();
            vm = VMDistribution(2.1,1.3);
            vmSysnoise = VMDistribution(0,0.3);
                        
            %% predict identity
            filter.setState(vm);
            filter.predictIdentity(vmSysnoise);
            vmIdentity = filter.getEstimate();
            testCase.verifyClass(vmIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmIdentity.mu);
            testCase.verifyGreaterThanOrEqual(vm.kappa, vmIdentity.kappa);
            
            %% predict nonlinear with identity as function
            filter.setState(vm);
            filter.predictNonlinear(@(x) x, vmSysnoise);
            vm2NonlinIdentity = filter.getEstimate();
            testCase.verifyClass(vm2NonlinIdentity, 'VMDistribution');
            testCase.verifyEqual(vmIdentity.mu, vm2NonlinIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmIdentity.kappa, vm2NonlinIdentity.kappa, 'RelTol', 1E-10);
            
            %% use linear system function
            vmPrior = VMDistribution(0, 0.8);
            filter.setState(vmPrior);
            filter.predictNonlinear(@(x) 2*x + 1, vmPrior);
            vmLinear = filter.getEstimate();
            testCase.verifyClass(vmLinear, 'VMDistribution');
            testCase.verifyEqual(vmLinear.mu, 2*vmPrior.mu + 1, 'RelTol', 1E-10);
            testCase.verifyLessThan(vmLinear.kappa, vmPrior.kappa);
                        
            %% prediction with non-additive noise
            filter.setState(vm);
            f = @(x,w) x + norm(w);
            wd = vmSysnoise.toDirac3();
            filter.predictNonlinearNonAdditive(f, wd.d, wd.w)
            vmNonadditive = filter.getEstimate();
            testCase.verifyClass(vmNonadditive, 'VMDistribution');
            testCase.verifyEqual(vmIdentity.mu, vmNonadditive.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmIdentity.kappa, vmNonadditive.kappa, 'RelTol', 1E-10);
        end
        
        function testUpdate(testCase)
            filter = VMFilter();
            vm = VMDistribution(2.1,1.3);
            measNoise = VMDistribution(0, 0.5);
            
            %% update identity
            filter.setState(vm);
            filter.updateIdentity(measNoise, vm.mu);
            vmIdentity = filter.getEstimate();
            testCase.verifyClass(vmIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmIdentity.mu);
            testCase.verifyLessThan(vm.kappa, vmIdentity.kappa);
            
            %% update identity with different measurement
            filter.setState(vm);
            filter.updateIdentity(measNoise, vm.mu+0.1);
            vmIdentity2 = filter.getEstimate();
            testCase.verifyClass(vmIdentity2, 'VMDistribution');
            testCase.verifyLessThan(vm.mu, vmIdentity2.mu);
            testCase.verifyLessThan(vm.kappa, vmIdentity2.kappa);
            
            %% nonlinear simple
            likelihood = LikelihoodFactory.additiveNoiseLikelihood(@(x) x, measNoise);
            filter.setState(vm);
            filter.updateNonlinear(likelihood,vm.mu);
            vmNonlinIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinIdentity.kappa, vmIdentity.kappa, 'RelTol', 0.2);
            
            %% nonlinear particle
            rng default
            filter.setState(vm);
            filter.updateNonlinearParticle(likelihood,vm.mu);
            vmNonlinParticleIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinParticleIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinParticleIdentity.mu, 'RelTol', 1E-1);
            testCase.verifyEqual(vmNonlinParticleIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-2);
            
            %% nonlinear progressive
            filter.setState(vm);
            filter.updateNonlinearProgressive(likelihood,vm.mu);
            vmNonlinProgressiveIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinProgressiveIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinProgressiveIdentity.mu, 'RelTol', 1E-5);
            testCase.verifyEqual(vmNonlinProgressiveIdentity.kappa, vmIdentity.kappa, 'RelTol', 0.2);
            
            %% updateNonlinearShift, 
            z = vm.mu;
            h = @(x) x;
            filter.setState(vm);
            filter.updateNonlinearShift(z,h,measNoise)
            vmNonlinShiftIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinShiftIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinShiftIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinShiftIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);
            
            %% updateNonlinearStatisticalShift 
            filter.setState(vm);
            filter.updateNonlinearStatisticalShift(z,h,measNoise)
            vmNonlinStatShiftIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinStatShiftIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinStatShiftIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinStatShiftIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);
            
            filter.setState(vm);
            filter.updateNonlinearStatisticalShift(z,h,measNoise,'dirac5')
            vmNonlinStatShiftIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinStatShiftIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinStatShiftIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinStatShiftIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);
            
            filter.setState(vm);
            filter.updateNonlinearStatisticalShift(z,h,measNoise,'dirac5',0.7)
            vmNonlinStatShiftIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinStatShiftIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinStatShiftIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinStatShiftIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);

            %% updateNonlinearCorrectedNoise
            filter.setState(vm);
            filter.updateNonlinearCorrectedNoise(z,h,measNoise)
            vmNonlinCorrectedNoiseIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinCorrectedNoiseIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinCorrectedNoiseIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinCorrectedNoiseIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);            
            
            filter.setState(vm);
            filter.updateNonlinearCorrectedNoise(z,h,measNoise,'dirac5')
            vmNonlinCorrectedNoiseIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinCorrectedNoiseIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinCorrectedNoiseIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinCorrectedNoiseIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);            
            
            filter.setState(vm);
            filter.updateNonlinearCorrectedNoise(z,h,measNoise,'dirac5',0.7)
            vmNonlinCorrectedNoiseIdentity = filter.getEstimate();
            testCase.verifyClass(vmNonlinCorrectedNoiseIdentity, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmNonlinCorrectedNoiseIdentity.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmNonlinCorrectedNoiseIdentity.kappa, vmIdentity.kappa, 'RelTol', 1E-10);            
        end
        
        function testAssociationLikelihood(testCase)
            vm1=VMDistribution(3,3);
            vm2=VMDistribution(1,10);
            vmfilter=VMFilter;
            vmfilter.setState(vm1);
            testCase.verifyEqual(vmfilter.associationLikelihood(vm2), vmfilter.associationLikelihoodNumerical(vm2), 'AbsTol', 1E-10);
        end
    end
end

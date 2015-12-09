classdef HypertoroidalVMSineDistributionTest< matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        function test2D(testCase)
            mu = [1, 2]';
            kappa = [0.7, 1.4]';
            Lambda = [0 0.4; 0.4 0];
            htvm = HypertoroidalVMSineDistribution(mu, kappa, Lambda);
            tvm = ToroidalVMSineDistribution(mu, kappa, Lambda(1,2));
            tvm2 = htvm.toToroidalVMSine();
            testCase.verifyClass(tvm2, 'ToroidalVMSineDistribution');
            
            % sanity check
            testCase.verifyClass(htvm, 'HypertoroidalVMSineDistribution');
            testCase.verifyEqual(htvm.mu, mu);
            testCase.verifyEqual(htvm.kappa, kappa);
            testCase.verifyEqual(htvm.Lambda, Lambda);

            % test pdf
            unnormalizedPdf = @(x) exp(kappa(1)*cos(x(1)-mu(1))+kappa(2)*cos(x(2)-mu(2)) + Lambda(1,2) * sin(x(1)-mu(1)).* sin(x(2)-mu(2)));
            T = htvm.T;
            pdf = @(x) unnormalizedPdf(x)/T;
            testCase.verifyEqual(htvm.pdf([3;2]), pdf([3;2]), 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf([1;4]), pdf([1;4]), 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf([5;6]), pdf([5;6]), 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf([-3;11]), pdf([-3;11]), 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf([5 1;6 3]), [pdf([5;6]) pdf([1;3])] , 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf([10 1 -5;11 -33 2]), [pdf([10;11]) pdf([1;-33]) pdf([-5; 2])] , 'RelTol', 1E-10);
            
            rng default
            testpoints = rand(2,100)*2*pi;
            testCase.verifyEqual(htvm.pdf(testpoints), tvm.pdf(testpoints), 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf(testpoints), tvm2.pdf(testpoints), 'RelTol', 1E-10);
            
            % test integral
            testCase.verifyEqual(htvm.integral(), 1, 'RelTol', 1E-5);
            testCase.verifyEqual(htvm.trigonometricMomentNumerical(0), [1;1], 'RelTol', 1E-5)
            
            % test toGaussian
            g = htvm.toGaussian();
            testCase.verifyClass(g, 'GaussianDistribution');
            testCase.verifyEqual(htvm.mu, g.mu);
        end
        
        function test1D(testCase)
            mu = 2;
            kappa = 0.7;
            Lambda = 0;
            htvm = HypertoroidalVMSineDistribution(mu, kappa, Lambda);
            vm1 = VMDistribution(mu, kappa);
            vm2 = htvm.toVM();
            testCase.verifyClass(vm2, 'VMDistribution');
            
            % sanity check
            testCase.verifyClass(htvm, 'HypertoroidalVMSineDistribution');
            testCase.verifyEqual(htvm.mu, mu);
            testCase.verifyEqual(htvm.kappa, kappa);
            testCase.verifyEqual(htvm.Lambda, Lambda);

            % test pdf
            rng default
            testpoints = rand(1,100)*2*pi;
            testCase.verifyEqual(htvm.pdf(testpoints), vm1.pdf(testpoints), 'RelTol', 1E-10);
            testCase.verifyEqual(htvm.pdf(testpoints), vm2.pdf(testpoints), 'RelTol', 1E-10);
            
            % test integral
            testCase.verifyEqual(htvm.integral(), 1, 'RelTol', 1E-5);
            testCase.verifyEqual(htvm.trigonometricMomentNumerical(0), 1, 'RelTol', 1E-5)
            
            % test toGaussian
            g = htvm.toGaussian();
            testCase.verifyClass(g, 'GaussianDistribution');
            testCase.verifyEqual(htvm.mu, g.mu);
        end        
    end
end
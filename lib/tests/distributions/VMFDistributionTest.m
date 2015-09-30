classdef VMFDistributionTest< matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        function testVMFDistribution2d(testCase)
            phi = 0.3;
            mu = [cos(phi); sin(phi)];            
            kappa = 2;
            vmf = VMFDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(vmf, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, mu);
            testCase.verifyEqual(vmf.kappa, kappa);
            testCase.verifyEqual(vmf.d, length(mu));
            
            %% test pdf
            x = 0:6;
            vmSame = VMDistribution(phi, kappa);
            testCase.verifyEqual(vmSame.pdf(x), vmf.pdf([cos(x);sin(x)]),  'RelTol', 1E-10);

            %% test integral
            testCase.verifyEqual(vmf.integral(), 1, 'RelTol', 1E-5);
            
            %% test fitting to samples
            % compare to results von von Mises distribution
            vm = VMDistribution(0,1);
            n = 2000;
            s = vm.sample(n);
            testCase.verifyEqual(size(s), [1, n]);
            wd = WDDistribution(s);
            vmFitted = wd.toVM();
            vmfFitted =  VMFDistribution.fit([cos(s); sin(s)]);
            testCase.verifyEqual([cos(vmFitted.mu);sin(vmFitted.mu)], vmfFitted.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmFitted.kappa, vmfFitted.kappa, 'RelTol', 1E-10);
            
            %% test multiplication
            other = VMFDistribution([cos(2*phi); sin(2*phi)], kappa/3);
            vmfMul = vmf.multiply(other);
            vmfMul2 = other.multiply(vmf);
            c = vmfMul.pdf([1;0])/(vmf.pdf([1;0])*other.pdf([1;0]));
            testCase.verifyEqual(vmf.pdf([cos(x);sin(x)]).*other.pdf([cos(x);sin(x)])*c, vmfMul.pdf([cos(x);sin(x)]), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.pdf([cos(x);sin(x)]).*other.pdf([cos(x);sin(x)])*c, vmfMul2.pdf([cos(x);sin(x)]), 'RelTol', 1E-10);
            
            %% test colvolve
            vmfConv = vmf.convolve(other);
            vmfConv2 = other.convolve(vmf);            
            testCase.verifyEqual(vmfConv.mu, vmf.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmfConv2.mu, other.mu, 'RelTol', 1E-10);
            d = 2;
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv2.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            
            %% test stochastic sampling
            n = 10;
            s = vmf.sample(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
            
            %% test entropy
            e = vmf.entropy();
            testCase.verifyEqual(e, vmf.entropyNumerical(), 'RelTol', 1E-10);
            testCase.verifyEqual(e, vmf.toVM().entropy(), 'RelTol', 1E-10);
        end
        
        function testVMFDistribution3d(testCase)
            mu = [1, 2, 3]';
            mu = mu/norm(mu);
            kappa = 2;
            vmf = VMFDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(vmf, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, mu);
            testCase.verifyEqual(vmf.kappa, kappa);
            testCase.verifyEqual(vmf.d, length(mu));

            %% test integral
            testCase.verifyEqual(vmf.integral(), 1, 'RelTol', 1E-5);
            
            %% test multiplication
            other = VMFDistribution([0; 0; 1], kappa/3);
            vmfMul = vmf.multiply(other);
            vmfMul2 = other.multiply(vmf);
            c = vmfMul.pdf([1;0;0])/(vmf.pdf([1;0;0])*other.pdf([1;0;0]));
            x = [0,1,0]';
            testCase.verifyEqual(vmf.pdf(x).*other.pdf(x)*c, vmfMul.pdf(x), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.pdf(x).*other.pdf(x)*c, vmfMul2.pdf(x), 'RelTol', 1E-10);
            
            %% test colvolve
            vmfConv = vmf.convolve(other);
            vmfConv2 = other.convolve(vmf);
            testCase.verifyEqual(vmfConv.mu, vmf.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmfConv2.mu, other.mu, 'RelTol', 1E-10);
            d = 3;
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv2.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            
            %% test stochastic sampling
            n = 10;
            s = vmf.sample(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
            
            s = vmf.sampleCdf(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
            
            s = vmf.sampleCdfVectorized(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
            
            s = vmf.sampleWood(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
            
            %% test entropy
            e = vmf.entropy();
            testCase.verifyEqual(e, vmf.entropyNumerical(), 'RelTol', 1E-7);
        end
        
        function testVMFDistribution4d(testCase)
            mu = [1, 2, 3, -2]';
            mu = mu/norm(mu);
            kappa = 1.6;
            vmf = VMFDistribution(mu, kappa);
            
            %% sanity check
            testCase.verifyClass(vmf, 'VMFDistribution');
            testCase.verifyEqual(vmf.mu, mu);
            testCase.verifyEqual(vmf.kappa, kappa);
            testCase.verifyEqual(vmf.d, length(mu));

            %% test integral
            testCase.verifyEqual(vmf.integral(), 1, 'RelTol', 1E-5);
            
            %% test multiplication
            other = VMFDistribution([0; 0; 0; 1], kappa/3);
            vmfMul = vmf.multiply(other);
            vmfMul2 = other.multiply(vmf);
            c = vmfMul.pdf([1;0;0;0])/(vmf.pdf([1;0;0;0])*other.pdf([1;0;0;0]));
            x = [0,1,0,0]';
            testCase.verifyEqual(vmf.pdf(x).*other.pdf(x)*c, vmfMul.pdf(x), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.pdf(x).*other.pdf(x)*c, vmfMul2.pdf(x), 'RelTol', 1E-10);
            
            %% test colvolve
            vmfConv = vmf.convolve(other);
            vmfConv2 = other.convolve(vmf);
            testCase.verifyEqual(vmfConv.mu, vmf.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmfConv2.mu, other.mu, 'RelTol', 1E-10);
            d = 4;
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv2.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);            
            
            %% test stochastic sampling
            n = 10;
            s = vmf.sample(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
        end
    end
end
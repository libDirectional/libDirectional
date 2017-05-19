classdef VMDistributionTest< matlab.unittest.TestCase    
    properties
    end
    
    methods (Test)
        function testVMDistribution(testCase)
            mu = 3;
            kappa = 1.5;
            vm = VMDistribution(mu,kappa);
            
            %% test pdf
            testCase.verifyEqual(vm.pdf(mu),1/2/pi/besseli(0,kappa) * exp(kappa), 'RelTol', 1E-10);
            testCase.verifyEqual(vm.pdf(mu+1.3),1/2/pi/besseli(0,kappa) * exp(kappa*cos(1.3)), 'RelTol', 1E-10);
            
            %% test cdf
            testCase.verifyEqual(vm.cdf(mu),vm.cdfNumerical(mu), 'RelTol', 1E-8);
            testCase.verifyEqual(vm.cdf(mu+1.3),vm.cdfNumerical(mu+1.3), 'RelTol', 1E-8);
            testCase.verifyEqual(vm.cdf(-10:10),vm.cdfNumerical(-10:10), 'RelTol', 1E-8);
            testCase.verifyEqual(vm.cdf(-10:10, 2),vm.cdfNumerical(-10:10, 2), 'RelTol', 1E-8);
            
            %% test integral
            testCase.verifyEqual(vm.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(vm.integral(0,pi)+vm.integral(pi,2*pi), 1, 'RelTol', 1E-6);
            
            %% test angular moments
            testCase.verifyEqual(vm.trigonometricMoment(1), vm.trigonometricMomentNumerical(1), 'RelTol', 1E-5)
            testCase.verifyEqual(vm.trigonometricMoment(2), vm.trigonometricMomentNumerical(2), 'RelTol', 1E-5)
            testCase.verifyEqual(vm.trigonometricMoment(3), vm.trigonometricMomentNumerical(3), 'RelTol', 1E-5)
            testCase.verifyEqual(vm.circularMean, mu);
            
            %% test fromMoment
            vmFromMoment = VMDistribution.fromMoment(vm.trigonometricMoment(1));
            testCase.verifyClass(vmFromMoment, 'VMDistribution');
            testCase.verifyEqual(vm.mu, vmFromMoment.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vm.kappa, vmFromMoment.kappa, 'RelTol', 1E-10);
            
            %% test conversions based on angular moment matching
            wd2 = vm.toDirac2();
            testCase.verifyClass(wd2, 'WDDistribution');
            testCase.verifyEqual(vm.trigonometricMoment(1), wd2.trigonometricMoment(1), 'RelTol', 1E-10);
            wd3 = vm.toDirac3();
            testCase.verifyClass(wd3, 'WDDistribution');
            testCase.verifyEqual(vm.trigonometricMoment(1), wd3.trigonometricMoment(1), 'RelTol', 1E-10);
            wd5 = vm.toDirac5();
            testCase.verifyClass(wd5, 'WDDistribution');
            testCase.verifyEqual(vm.trigonometricMoment(1), wd5.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyEqual(vm.trigonometricMoment(2), wd5.trigonometricMoment(2), 'RelTol', 1E-10);
            wn = vm.toWN();
            testCase.verifyClass(wn, 'WNDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-10);
            
            %% test multiplication
            vmMul = vm.multiply(vm);
            testCase.verifyClass(vmMul, 'VMDistribution');
            testCase.verifyEqual(vmMul.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(vmMul.mu, vm.mu, 'RelTol', 1E-10);
            testCase.verifyGreaterThan(vmMul.kappa, vm.kappa);
            
            % VM is closed under multiplication, i.e., new pdf is product of
            % old pdfs times a renormalization factor
            x = 0:6;
            c = vmMul.pdf(0)/(vm.pdf(0)^2);
            testCase.verifyEqual(vmMul.pdf(x), c * vm.pdf(x).^2, 'RelTol', 1E-10);
            
            %% test convolution
            vmConv = vm.convolve(vm);
            testCase.verifyClass(vmConv, 'VMDistribution');
            testCase.verifyEqual(vmConv.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(vmConv.mu, vm.mu + vm.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmConv.trigonometricMoment(1), vm.trigonometricMoment(1)*vm.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyLessThan(vmConv.kappa, vm.kappa);
            
            %% test entropy
            testCase.verifyEqual(vm.entropy(), vm.entropyNumerical(), 'RelTol', 1E-10);
            
            %% test sampling
            n = 10;
            s = vm.sample(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));
        end
        
        function testKld(testCase)
            rng default
            for i=1:20
                vm1 = VMDistribution(2*pi*rand(1), abs(randn(1)));
                vm2 = VMDistribution(2*pi*rand(1), abs(randn(1)));
                testCase.verifyEqual(vm1.kld(vm2), vm1.kldNumerical(vm2), 'RelTol', 1E-10);
                testCase.verifyEqual(vm2.kld(vm1), vm2.kldNumerical(vm1), 'RelTol', 1E-10);
            end
        end
        
        function testUniformSpecialCase(testCase)
            vm = VMDistribution(1,0);
            cu = CircularUniformDistribution();

            % test pdf
            x = 0:0.1:2*pi;
            testCase.verifyEqual(vm.pdf(x),cu.pdf(x) , 'RelTol', 1E-10);
            
            % test cdf
            testCase.verifyEqual(vm.cdf(x),vm.cdfNumerical(x), 'RelTol', 1E-8);
            testCase.verifyEqual(vm.cdf(-10:10),vm.cdfNumerical(-10:10), 'RelTol', 1E-8);
            testCase.verifyEqual(vm.cdf(-10:10, 2),vm.cdfNumerical(-10:10, 2), 'RelTol', 1E-8);
            
            % test integral
            testCase.verifyEqual(vm.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(vm.integral(0,pi)+vm.integral(pi,2*pi), 1, 'RelTol', 1E-6);
            
            % test angular moments
            testCase.verifyEqual(vm.trigonometricMoment(1), cu.trigonometricMoment(1), 'RelTol', 1E-5)
            testCase.verifyEqual(vm.trigonometricMoment(2), cu.trigonometricMoment(2), 'RelTol', 1E-5)
            testCase.verifyEqual(vm.trigonometricMoment(3), cu.trigonometricMoment(3), 'RelTol', 1E-5)
            
            % test fromMoment
            vmFromMoment = VMDistribution.fromMoment(0);
            testCase.verifyClass(vmFromMoment, 'VMDistribution');
            testCase.verifyEqual(vmFromMoment.kappa, 0, 'RelTol', 1E-10);
            
            % test conversions based on angular moment matching
            wd2 = vm.toDirac2();
            testCase.verifyClass(wd2, 'WDDistribution');
            testCase.verifyEqual(vm.trigonometricMoment(1), wd2.trigonometricMoment(1), 'AbsTol', 1E-10);
            wd3 = vm.toDirac3();
            testCase.verifyClass(wd3, 'WDDistribution');
            testCase.verifyEqual(vm.trigonometricMoment(1), wd3.trigonometricMoment(1), 'AbsTol', 1E-10);
            wd5 = vm.toDirac5();
            testCase.verifyClass(wd5, 'WDDistribution');
            testCase.verifyEqual(vm.trigonometricMoment(1), wd5.trigonometricMoment(1), 'AbsTol', 1E-10);
            testCase.verifyEqual(vm.trigonometricMoment(2), wd5.trigonometricMoment(2), 'AbsTol', 1E-10);
            wn = vm.toWN();
            testCase.verifyClass(wn, 'WNDistribution');
            testCase.verifyEqual(wn.trigonometricMoment(1), vm.trigonometricMoment(1), 'AbsTol', 1E-10);
            
            % test multiplication
            vmMul = vm.multiply(vm);
            testCase.verifyClass(vmMul, 'VMDistribution');
            testCase.verifyEqual(vmMul.integral, 1, 'RelTol', 1E-10);
            testCase.verifyGreaterThanOrEqual(vmMul.kappa, vm.kappa);
            
            % VM is closed under multiplication, i.e., new pdf is product of
            % old pdfs times a renormalization factor
            x = 0:6;
            c = vmMul.pdf(0)/(vm.pdf(0)^2);
            testCase.verifyEqual(vmMul.pdf(x), c * vm.pdf(x).^2, 'RelTol', 1E-10);
            
            % test convolution
            vmConv = vm.convolve(vm);
            testCase.verifyClass(vmConv, 'VMDistribution');
            testCase.verifyEqual(vmConv.integral, 1, 'RelTol', 1E-10);
            testCase.verifyEqual(vmConv.trigonometricMoment(1), vm.trigonometricMoment(1)*vm.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyLessThanOrEqual(vmConv.kappa, vm.kappa);
            
            % test entropy
            testCase.verifyEqual(vm.entropy(), cu.entropy(), 'RelTol', 1E-10);
            
            % test sampling
            n = 10;
            s = vm.sample(n);
            testCase.verifyEqual(size(s,1), 1);
            testCase.verifyEqual(size(s,2), n);
            testCase.verifyGreaterThanOrEqual(s,zeros(size(s)));
            testCase.verifyLessThan(s,2*pi*ones(size(s)));            
        end
        
        function testQuantization(testCase)
            vm = VMDistribution(3,1.7);
            n = 11;
            [s,w] = vm.sampleOptimalQuantization(n);
            testCase.verifyEqual(sum(w), 1, 'RelTol', 1E-10);
            testCase.verifySize(s, [1 n]);
            testCase.verifySize(w, [1 n]);
            wd = WDDistribution(s,w);
            testCase.verifyEqual(wd.trigonometricMoment(1), vm.trigonometricMoment(1), 'RelTol', 1E-2);
        end        
    end
end
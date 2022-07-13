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
            testCase.verifyEqual(vmf.dim, length(mu));
            
            %% test pdf
            x = 0:6;
            vmSame = VMDistribution(phi, kappa);
            testCase.verifyEqual(vmSame.pdf(x), vmf.pdf([cos(x);sin(x)]),  'RelTol', 1E-10);
            
            %% test mode
            testCase.verifyEqual(vmf.modeNumerical, vmf.mode, 'AbsTol',1E-5)

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
            
            vmfFittedScore =  VMFDistribution.fitScoreBased([cos(s); sin(s)]);
            testCase.verifyEqual([cos(vmFitted.mu);sin(vmFitted.mu)], vmfFittedScore.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmFitted.kappa, vmfFittedScore.kappa, 'RelTol', 1E-1);
            
            %% test deterministic sampling
            s = vmf.sampleDeterministic();
            testCase.verifyEqual(size(s,1), vmf.dim);
            testCase.verifyEqual(size(s,2), 2*vmf.dim-1);
            sVM = vmSame.toDirac3().d;
            %we need sort because VM and VMF return the samples in a different order
            testCase.verifyEqual(sort(s,2), sort([cos(sVM); sin(sVM)],2), 'RelTol',1E-10); 
            vmfFitted2 = VMFDistribution.fit(s);
            testCase.verifyEqual(vmf.mu, vmfFitted2.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.kappa, vmfFitted2.kappa, 'RelTol', 1E-10);
            
            %% test multiplication
            other = VMFDistribution([cos(2*phi); sin(2*phi)], kappa/3);
            vmfMul = vmf.multiply(other);
            vmfMul2 = other.multiply(vmf);
            c = vmfMul.pdf([1;0])/(vmf.pdf([1;0])*other.pdf([1;0]));
            testCase.verifyEqual(vmf.pdf([cos(x);sin(x)]).*other.pdf([cos(x);sin(x)])*c, vmfMul.pdf([cos(x);sin(x)]), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.pdf([cos(x);sin(x)]).*other.pdf([cos(x);sin(x)])*c, vmfMul2.pdf([cos(x);sin(x)]), 'RelTol', 1E-10);
            
            %% test colvolve
            other=VMFDistribution([0;1], kappa/3);
            vmfConv = vmf.convolve(other);
            testCase.verifyEqual(vmfConv.mu, vmf.mu, 'RelTol', 1E-10);
            d = 2;
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            
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
            testCase.verifyEqual(vmf.dim, length(mu));

            %% test mode
            testCase.verifyEqual(vmf.modeNumerical, vmf.mode, 'AbsTol',1E-5)

            %% test integral
            testCase.verifyEqual(vmf.integral(), 1, 'RelTol', 1E-5);
                        
            %% test deterministic sampling
            s = vmf.sampleDeterministic();
            testCase.verifyEqual(size(s,1), vmf.dim);
            testCase.verifyEqual(size(s,2), 2*vmf.dim-1);
            vmfFitted2 = VMFDistribution.fit(s);
            testCase.verifyEqual(vmf.mu, vmfFitted2.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmf.kappa, vmfFitted2.kappa, 'RelTol', 1E-10);
            
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
            testCase.verifyEqual(vmfConv.mu, vmf.mu, 'RelTol', 1E-10);
            d = 3;
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            
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
            testCase.verifyEqual(vmf.dim, length(mu));

            %% test mode
            testCase.verifyEqual(vmf.modeNumerical, vmf.mode, 'AbsTol',1E-5)

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
            testCase.verifyEqual(vmfConv.mu, vmf.mu, 'RelTol', 1E-10);
            d = 4;
            testCase.verifyEqual(VMFDistribution.Ad(d, vmfConv.kappa), VMFDistribution.Ad(d, vmf.kappa)*VMFDistribution.Ad(d, other.kappa), 'RelTol', 1E-10);
            
            %% test stochastic sampling
            n = 10;
            s = vmf.sample(n);
            testCase.verifySize(s, [length(mu), n])
            testCase.verifyEqual(sum(s.^2), ones(1,n), 'RelTol', 1E-10);
            
            %% test toGaussian
            g = vmf.toGaussian();
            testCase.verifyClass(g, 'GaussianDistribution');
            testCase.verifyEqual(g.mu, vmf.mu);
        end
        
        function testMeanDirection(testCase)
            mu = 1/sqrt(2)*[1;1;0];
            vmf = VMFDistribution(mu,1);
            testCase.verifyEqual(vmf.meanDirection, mu, 'AbsTol', 1e-13);
        end
        
        function testHellinger(testCase)
            % 2D
            vmf1 = VMFDistribution([1,0]', 0.9);
            vmf2 = VMFDistribution([0,1]', 1.7);
            testCase.verifyEqual(vmf1.hellingerDistance(vmf1), 0, 'AbsTol', 1E-10);
            testCase.verifyEqual(vmf2.hellingerDistance(vmf2), 0, 'AbsTol', 1E-10);
            testCase.verifyEqual(vmf1.hellingerDistance(vmf2), vmf1.hellingerDistanceNumerical(vmf2), 'RelTol', 1E-10);
            testCase.verifyEqual(vmf1.hellingerDistance(vmf2), vmf2.hellingerDistance(vmf1), 'RelTol', 1E-10);
            
            % 3D
            vmf1 = VMFDistribution([1,0,0]', 0.6);
            mu2 = [1,2,3]';
            vmf2 = VMFDistribution(mu2/norm(mu2), 2.1);
            testCase.verifyEqual(vmf1.hellingerDistance(vmf1), 0, 'AbsTol', 1E-10);
            testCase.verifyEqual(vmf2.hellingerDistance(vmf2), 0, 'AbsTol', 1E-10);
            testCase.verifyEqual(vmf1.hellingerDistance(vmf2), vmf1.hellingerDistanceNumerical(vmf2), 'RelTol', 1E-6);
            testCase.verifyEqual(vmf1.hellingerDistance(vmf2), vmf2.hellingerDistance(vmf1), 'RelTol', 1E-10);
        end
        
        function testMoment2D(testCase)
            vmf1 = VMFDistribution([sqrt(2)/2,-sqrt(2)/2]', 0.9);
            m = vmf1.moment();
            testCase.verifyEqual(m/norm(m), vmf1.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(norm(m), besseli(1,vmf1.kappa)/besseli(0,vmf1.kappa), 'RelTol', 1E-10);
            mCirc = vmf1.toVM().trigonometricMoment(1);
            testCase.verifyEqual(m, [real(mCirc), imag(mCirc)]', 'RelTol', 1E-10);
            vmf2 = VMFDistribution.fromMoment(m);
            testCase.verifyEqual(vmf1.mu, vmf2.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmf1.kappa, vmf2.kappa, 'RelTol', 1E-10);
        end
        
        function testMoment3D(testCase)    
            vmf1 = VMFDistribution([-1, 0, 0]', 1.3);
            m = vmf1.moment();
            testCase.verifyEqual(m/norm(m), vmf1.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(norm(m), besseli(1.5,vmf1.kappa)/besseli(0.5,vmf1.kappa), 'RelTol', 1E-10);
            vmf2 = VMFDistribution.fromMoment(m);
            testCase.verifyEqual(vmf1.mu, vmf2.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmf1.kappa, vmf2.kappa, 'RelTol', 1E-10);
        end
        
        function testMoment4D(testCase)    
            vmf1 = VMFDistribution([0, 1, 0, 0]', 0.5);
            m = vmf1.moment();
            testCase.verifyEqual(m/norm(m), vmf1.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(norm(m), besseli(2,vmf1.kappa)/besseli(1,vmf1.kappa), 'RelTol', 1E-10);
            vmf2 = VMFDistribution.fromMoment(m);
            testCase.verifyEqual(vmf1.mu, vmf2.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(vmf1.kappa, vmf2.kappa, 'RelTol', 1E-10);
        end
        
        function testShifting(testCase)
            dist = VMFDistribution([0; 0; 0; 1], 0.5);
            distShifted = dist.shift([1;2;3;4]/norm([1;2;3;4]));
            testCase.verifyEqual(distShifted.mu, [1;2;3;4]/norm([1;2;3;4]), 'RelTol', 1E-10);
        end
    end
end
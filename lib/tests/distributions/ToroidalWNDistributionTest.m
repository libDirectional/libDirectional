classdef ToroidalWNDistributionTest < matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        function testToroidalWNDistribution(testCase)
            mu = [1, 2]';
            C = [1.3 -0.9;
                -0.9 1.2];
            twn = ToroidalWNDistribution(mu, C);
            [X,Y] = meshgrid(1:6, 1:6);
            testpoints = [X(:)'; Y(:)'];
            
            %% sanity check
            testCase.verifyClass(twn, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn.mu, mu);
            testCase.verifyEqual(twn.C, C);

            %% test integral
            testCase.verifyEqual(twn.integral(), 1, 'RelTol', 1E-5);
            testCase.verifyEqual(twn.trigonometricMoment(0), [1;1], 'RelTol', 1E-5)
            
            %% test pdf            
            function p = matlabpdf(dist, xa, n)
                % matlab implementation
                if size(xa,2)>1
                    p = zeros(1,size(xa,2));
                    for i=1:size(xa,2)
                        p(i) = matlabpdf(dist, xa(:,i),n);
                    end
                    return
                end
                total = 0;
                Cinv = -0.5 *inv(dist.C); %#ok<MINV>
                xa = xa-dist.mu;
                for j=-n:n
                    for k=-n:n
                        dxa = [2*pi*j; 2*pi*k];
                        total = total + exp((xa+dxa)'*Cinv*(xa+dxa));
                    end
                end
                p = total/2/pi/sqrt(det(dist.C));
            end
            testCase.verifyEqual(matlabpdf(twn, testpoints, 3), twn.pdf(testpoints,3), 'RelTol', 1E-10);
            
            %% test angular moment
            m = twn.trigonometricMoment(1);
            m1 = twn.marginalizeTo1D(1).trigonometricMoment(1);
            m2 = twn.marginalizeTo1D(2).trigonometricMoment(1);
            testCase.verifyEqual(m(1), m1, 'RelTol', 1E-10);
            testCase.verifyEqual(m(2), m2, 'RelTol', 1E-10);
            testCase.verifyEqual(m, twn.trigonometricMomentNumerical(1), 'RelTol', 1E-9);
            
            %% test getMean
            mean = twn.mean4D();
            testCase.verifyEqual(mean(1), real(m1), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(m1), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(3), real(m2), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(4), imag(m2), 'RelTol', 1E-10);
            testCase.verifyEqual(twn.circularMean, mu);
            
            %% test covariance
            testCase.verifyEqual(twn.covariance4D(), twn.covariance4DNumerical(), 'RelTol', 1E-1); %very imprecise?
            
            %% test multiplication
            twnMul = twn.multiplyMomentBased(twn);
            testCase.verifyClass(twnMul, 'ToroidalWNDistribution');
            testCase.verifyEqual(twnMul.mu, twn.mu, 'RelTol', 1E-9);
            testCase.verifyLessThan(det(twnMul.C), det(twn.C));
            
            %% convolution
            twnConv = twn.convolve(twn);
            testCase.verifyClass(twnConv, 'ToroidalWNDistribution');
            testCase.verifyEqual(twnConv.integral, 1, 'RelTol', 1E-6);
            testCase.verifyEqual(twnConv.mu, twn.mu + twn.mu, 'RelTol', 1E-10);
            testCase.verifyEqual(twnConv.trigonometricMoment(1), twn.trigonometricMoment(1).*twn.trigonometricMoment(1), 'RelTol', 1E-10);
            testCase.verifyGreaterThan(det(twnConv.C), det(twn.C));
            testCase.verifyGreaterThan(eig(twnConv.C-twn.C), [0;0]); %matrix is greater in the sense that difference is pos. definite
            
            fd = ToroidalFourierDistribution.fromDistribution(twn,45);
            fdConv = fd.convolve(fd);
            testCase.verifyEqual(twnConv.pdf(testpoints), fdConv.pdf(testpoints),'AbsTol', 1E-7);
            testCase.verifyEqual(twnConv.trigonometricMoment(1), fdConv.trigonometricMoment(1), 'RelTol', 1E-7);
            
            %% test sampling
            nSamples = 5;
            s = twn.sample(nSamples);
            testCase.verifyEqual(size(s,1),2);
            testCase.verifyEqual(size(s,2),nSamples);
            testCase.verifyEqual(s, mod(s,2*pi));
            
            %% test correlation coefficients
            testCase.verifyEqual(twn.realCorrelation(), C(1,2)/sqrt(C(1,1)*C(2,2)),'RelTol', 1E-10);
            testCase.verifyEqual(twn.circularCorrelationJammalamadaka(), twn.circularCorrelationJammalamadakaNumerical(),'RelTol', 1E-6)
            
            %% test conversions
            tvm = twn.toToroidalVMSine();
            testCase.verifyClass(tvm, 'ToroidalVMSineDistribution');
            twn2 = tvm.toToroidalWN();
            testCase.verifyClass(twn2, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn.mu, twn2.mu, 'RelTol', 1E-5)
            testCase.verifyEqual(twn.C, twn2.C, 'RelTol', 1E-5)
            twd5 = twn.toToroidalWD5();
            testCase.verifyClass(twd5, 'ToroidalWDDistribution');
            testCase.verifyEqual(twd5.trigonometricMoment(1), twn.trigonometricMoment(1), 'RelTol', 1E-5)
            testCase.verifyEqual(twd5.circularCorrelationJammalamadaka(), twn.circularCorrelationJammalamadaka(), 'RelTol', 1E-5)
            twn3 = twd5.toToroidalWN();
            testCase.verifyEqual(twn.mu, twn3.mu, 'RelTol', 1E-5)
            testCase.verifyEqual(twn.C, twn3.C, 'RelTol', 1E-5)
                        
            
            %% test getMarginal
            wn1 = twn.marginalizeTo1D(1);
            wn2 = twn.marginalizeTo1D(2);
            testCase.verifyClass(wn1, 'WNDistribution');
            testCase.verifyClass(wn2, 'WNDistribution');
            testCase.verifyEqual(wn1.mu, twn.mu(1),'RelTol', 1E-10);
            testCase.verifyEqual(wn2.mu, twn.mu(2),'RelTol', 1E-10);
            testCase.verifyEqual(wn1.sigma^2, twn.C(1,1),'RelTol', 1E-10);
            testCase.verifyEqual(wn2.sigma^2, twn.C(2,2),'RelTol', 1E-10);
            
            %% test MLE
            rng default
            samples = twn.sample(2000);
            twnMleJensen = ToroidalWNDistribution.mleJensen(samples);
            twnMleNumerical = ToroidalWNDistribution.mleNumerical(samples);
            testCase.verifyClass(twnMleJensen, 'ToroidalWNDistribution');
            testCase.verifyClass(twnMleNumerical, 'ToroidalWNDistribution');
            testCase.verifyEqual(twn.mu, twnMleJensen.mu, 'RelTol', 1E-1)
            testCase.verifyEqual(twn.C, twnMleJensen.C, 'RelTol', 0.4)
            testCase.verifyEqual(twn.mu, twnMleNumerical.mu, 'RelTol', 0.02)
            testCase.verifyEqual(twn.C, twnMleNumerical.C, 'RelTol', 1E-1)
            
            %% test loglikelihood
            l = twn.logLikelihood(samples);
            l2 = sum(log(twn.pdf(samples)));
            testCase.verifyEqual(l, l2, 'RelTol', 1E-10);
        end
        
        function testShift(testCase)
            twn = ToroidalWNDistribution([3;5],eye(2));
            s = [4;2];
            twn2 = twn.shift(s);
            testCase.verifyClass(twn2, 'ToroidalWNDistribution');
            [xTest,yTest]=meshgrid(linspace(0,2*pi,20));
            testCase.verifyEqual(twn2.pdf([xTest(:)';yTest(:)']),twn.pdf([xTest(:)' - s(1);yTest(:)' - s(2)]),'AbsTol',1E-10);
        end        
    end
end
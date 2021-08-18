classdef ToroidalVMRivestDistributionTest< matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        function testToroidalVMRivestDistribution(testCase)
            mu = [1, 2]';
            kappa = [0.7, 1.4]';
            alpha = 0.3;
            beta = 0.5;
            tvm = ToroidalVMRivestDistribution(mu, kappa, alpha, beta);
            
            %% sanity check
            testCase.verifyClass(tvm, 'ToroidalVMRivestDistribution');
            testCase.verifyEqual(tvm.mu, mu);
            testCase.verifyEqual(tvm.kappa, kappa);
            testCase.verifyEqual(tvm.alpha, alpha);
            testCase.verifyEqual(tvm.beta, beta);

            %% test integral
            testCase.verifyEqual(tvm.integral(), 1, 'RelTol', 1E-5);
            testCase.verifyEqual(tvm.trigonometricMomentNumerical(0), [1;1], 'RelTol', 1E-5)
                        
            %% test pdf
            unnormalizedPdf = @(x) exp(kappa(1)*cos(x(1)-mu(1))+kappa(2)*cos(x(2)-mu(2)) + alpha * cos(x(1)-mu(1)).* cos(x(2)-mu(2)) + beta * sin(x(1)-mu(1)).* sin(x(2)-mu(2)));
            C = tvm.C;
            pdf = @(x) unnormalizedPdf(x)*C;
            testCase.verifyEqual(tvm.pdf([3;2]), pdf([3;2]), 'RelTol', 1E-10);
            testCase.verifyEqual(tvm.pdf([1;4]), pdf([1;4]), 'RelTol', 1E-10);
            testCase.verifyEqual(tvm.pdf([5;6]), pdf([5;6]), 'RelTol', 1E-10);
            testCase.verifyEqual(tvm.pdf([-3;11]), pdf([-3;11]), 'RelTol', 1E-10);
            testCase.verifyEqual(tvm.pdf([5 1;6 3]), [pdf([5;6]) pdf([1;3])] , 'RelTol', 1E-10);
            
            %% test conversion to ToroidalVMMatrixDistribution
            tvmMatrix = tvm.toToroidalVMMatrixDistribution;
            testpoints = mod([1:20;6:25],2*pi);
            testCase.verifyEqual(tvm.pdf(testpoints), tvmMatrix.pdf(testpoints), 'RelTol', 1E-3); %error is quite large due to inaccurate normalization constant
                        
            %% test sampling
            nSamples = 5;
            s = tvm.sample(nSamples);
            testCase.verifyEqual(size(s,1),2);
            testCase.verifyEqual(size(s,2),nSamples);
            testCase.verifyEqual(s, mod(s,2*pi));
                        
            %% test trigonometric moment
            m = tvm.trigonometricMoment(1);
            testCase.verifyEqual(m, tvm.trigonometricMomentNumerical(1), 'RelTol', 1E-8);
                        
            %% test circularCorrelationJammalamadaka
            testCase.verifyEqual(tvm.circularCorrelationJammalamadaka(), tvm.circularCorrelationJammalamadakaNumerical(), 'RelTol', 1E-8);
                        
            %% test mean4D
            mean = tvm.mean4D();
            testCase.verifyEqual(mean(1), real(m(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(2), imag(m(1)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(3), real(m(2)), 'RelTol', 1E-10);
            testCase.verifyEqual(mean(4), imag(m(2)), 'RelTol', 1E-10);
            testCase.verifyEqual(tvm.circularMean, mu, 'RelTol', 1E-5);         
            
            return
                                               
            %% test conversions based on angular moment matching
            tvm2 = ToroidalVMCosineDistribution(mu, [1.8, 1.6]', 0.2);
            twn = tvm2.toToroidalWN();
            testCase.verifyClass(twn, 'ToroidalWNDistribution');
            %tvm2 = twn.toToroidalVMSine();
            %testCase.verifyClass(tvm2, 'ToroidalVMSineDistribution');
            %testCase.verifyEqual(tvm.mu, tvm2.mu, 'RelTol', 1E-5)
            %testCase.verifyEqual(tvm.kappa, tvm2.kappa, 'RelTol', 1E-5)
            %testCase.verifyEqual(tvm.lambda, tvm2.lambda, 'RelTol', 1E-5)
        end
        
        function testShift(testCase)
            tvm = ToroidalVMRivestDistribution([3;5],[0.7; 1.3], 0.4, 0.3);
            s = [4;2];
            tvm2 = tvm.shift(s);
            testCase.verifyClass(tvm2, 'ToroidalVMRivestDistribution');
            [xTest,yTest]=meshgrid(linspace(0,2*pi,20));
            testCase.verifyEqual(tvm2.pdf([xTest(:)';yTest(:)']),tvm.pdf([xTest(:)' - s(1);yTest(:)' - s(2)]),'AbsTol',1E-10);
        end           
    end
end
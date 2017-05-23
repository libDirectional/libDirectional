classdef ToroidalVMMatrixDistributionTest< matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        function testToroidalVMDistribution(testCase)
            mu = [1.7, 2.3]';
            kappa = [0.8, 0.3]';
            A = 0.1*[1 2; 3 4];
            tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
            
            [X,Y] = meshgrid(1:6, 1:6);
            testpoints = [X(:)'; Y(:)'];
            
            %% sanity check
            testCase.verifyClass(tvm, 'ToroidalVMMatrixDistribution');
            testCase.verifyEqual(tvm.mu, mu);
            testCase.verifyEqual(tvm.kappa, kappa);
            testCase.verifyEqual(tvm.A, A);
            
            %tvm.plot();
            
            %% test pdf
            function p = pdf(xa, mu, kappa, A, C)
                if size(xa,2)>1
                    p = zeros(1,size(xa,2));
                    for i=1:size(xa,2)
                        p(1,i) = pdf(xa(:,i), mu, kappa, A, C);
                    end
                    return
                end
                p = C * exp ( ...
                      kappa(1) * cos(xa(1,1) - mu(1)) ...
                    + kappa(2) * cos(xa(2,1) - mu(2)) ...
                    +   [cos(xa(1,1) - mu(1)) , sin(xa(1,1) - mu(1))] ...
                      * A ...
                      * [cos(xa(2,1) - mu(2)) ; sin(xa(2,1) - mu(2))]);
            end
            
            testCase.verifyEqual(tvm.pdf(testpoints), pdf(testpoints, mu, kappa, A, tvm.C),  'RelTol', 1E-10);
            

            %% test integral
            testCase.verifyEqual(tvm.integral(), 1, 'RelTol', 1E-5);
            testCase.verifyEqual(tvm.trigonometricMomentNumerical(0), [1;1], 'RelTol', 1E-5)
            
            %% multiply
            tvm2 = ToroidalVMMatrixDistribution(mu+1, kappa+0.2, inv(A));
            tvmMul = tvm.multiply(tvm2);
            tvmMulSwapped = tvm2.multiply(tvm);
            
            C = tvm.pdf([0;0])*tvm2.pdf([0;0])/ tvmMul.pdf([0;0]); %renormalization factor
            testCase.verifyEqual(tvm.pdf(testpoints).*tvm2.pdf(testpoints), C*tvmMul.pdf(testpoints), 'RelTol', 1E-10);
            testCase.verifyEqual(tvm.pdf(testpoints).*tvm2.pdf(testpoints), C*tvmMulSwapped.pdf(testpoints), 'RelTol', 1E-10);
            
            % compare to ToroidalFourier multiplication
            n = 45; %use a lot of coefficients for high accuracy
            tf = ToroidalFourierDistribution.fromDistribution(tvm, n);
            tf2 = ToroidalFourierDistribution.fromDistribution(tvm2, n);
            tfMul = tf.multiply(tf2);
            tfMulSwapped = tf2.multiply(tf);
            %fix normconst (todo should be removed once normconst works
            %properly in this range of values)
            tvmMul.C = 1;
            tvmMul.C = 1/tvmMul.integral();
            tvmMulSwapped.C = 1;
            tvmMulSwapped.C = 1/tvmMul.integral();
            % use abstol instead of reltol, because relative accuracy is bad in areas with
            % very little probability mass
            testCase.verifyEqual(tfMul.pdf(testpoints), tvmMul.pdf(testpoints), 'AbsTol', 1E-5);
            testCase.verifyEqual(tfMulSwapped.pdf(testpoints), tvmMul.pdf(testpoints), 'AbsTol', 1E-5);
            
            %% derivatives of normalization constant
            % compare to finite differences
            CinvDiff = tvm.normConstApproxDiff();
            epsilon = 1E-9;
            tvmKappa1Plus = ToroidalVMMatrixDistribution(mu, kappa + epsilon*[1;0], A);
            testCase.verifyEqual(CinvDiff(1), (1/tvmKappa1Plus.C - 1/tvm.C)/epsilon, 'RelTol', 1E-5);
            tvmKappa2Plus = ToroidalVMMatrixDistribution(mu, kappa + epsilon*[0;1], A);
            testCase.verifyEqual(CinvDiff(2), (1/tvmKappa2Plus.C - 1/tvm.C)/epsilon, 'RelTol', 1E-5);
            tvmA11Plus = ToroidalVMMatrixDistribution(mu, kappa, A + epsilon*[1,0;0,0]);
            testCase.verifyEqual(CinvDiff(3), (1/tvmA11Plus.C - 1/tvm.C)/epsilon, 'RelTol', 1E-5);
            tvmA12Plus = ToroidalVMMatrixDistribution(mu, kappa, A + epsilon*[0,1;0,0]);
            testCase.verifyEqual(CinvDiff(4), (1/tvmA12Plus.C - 1/tvm.C)/epsilon, 'RelTol', 1E-5);
            tvmA21Plus = ToroidalVMMatrixDistribution(mu, kappa, A + epsilon*[0,0;1,0]);
            testCase.verifyEqual(CinvDiff(5), (1/tvmA21Plus.C - 1/tvm.C)/epsilon, 'RelTol', 1E-5);
            tvmA22Plus = ToroidalVMMatrixDistribution(mu, kappa, A + epsilon*[0,0;0,1]);
            testCase.verifyEqual(CinvDiff(6), (1/tvmA22Plus.C - 1/tvm.C)/epsilon, 'RelTol', 1E-5);
            
            %% derivatives of logLikelihood
            twn = ToroidalWNDistribution(tvm.mu, [1.1 0.7; 0.7 1.3]);
            rng default
            samples = twn.sample(100);
            testCase.verifyEqual( tvm.logLikelihoodDiffNumerical(samples), tvm.logLikelihoodDiff(samples), 'RelTol', 1E-7);
            testCase.verifyEqual( tvm.logLikelihoodDiffNumerical(samples(:,1)), tvm.logLikelihoodDiff(samples(:,1)), 'RelTol', 1E-7);
                        
            %% test marginals
            for dim = 1:2
                dist(dim) = tvm.marginalizeTo1D(1); %#ok<AGROW>
                testCase.verifyClass(dist(dim), 'CustomCircularDistribution');
                testCase.verifyEqual(dist(dim).integral(), 1, 'RelTol', 1E-6);
            end
            marginal1Numerical = @(x) integral(@(y) tvm.pdf([repmat(x,1,length(y));y]), 0, 2*pi);
            testCase.verifyEqual(dist(1).pdf(1), marginal1Numerical(1), 'RelTol', 1E-10);
            testCase.verifyEqual(dist(1).pdf(3), marginal1Numerical(3), 'RelTol', 1E-10);
            testCase.verifyEqual(dist(1).pdf(5), marginal1Numerical(5), 'RelTol', 1E-10);
            marginal2Numerical = @(x) integral(@(y) tvm.pdf([repmat(x,1,length(y));y]), 0, 2*pi);
            testCase.verifyEqual(dist(2).pdf(1), marginal2Numerical(1), 'RelTol', 1E-10);
            testCase.verifyEqual(dist(2).pdf(3), marginal2Numerical(3), 'RelTol', 1E-10);
            testCase.verifyEqual(dist(2).pdf(5), marginal2Numerical(5), 'RelTol', 1E-10);
        end
        
        function testMle(testCase)
            mu = [1.7, 2.3]';
            kappa = [1.4, 1.2]';
            A = 0*[1 2; 1 2];
            tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
            
            rng default
            samples = tvm.sample(500);
            
            tvmMle = ToroidalVMMatrixDistribution.mleNumerical(samples);
            testCase.verifyClass(tvmMle, 'ToroidalVMMatrixDistribution');           
            testCase.verifyEqual(tvm.mu, tvmMle.mu, 'RelTol', 0.03)
            testCase.verifyEqual(tvm.kappa, tvmMle.kappa, 'RelTol', 0.2)
            testCase.verifyEqual(tvm.A, tvmMle.A, 'AbsTol', 0.3)
            
            tvmMle = ToroidalVMMatrixDistribution.mleNumerical(samples, 'quasi-newton');
            testCase.verifyClass(tvmMle, 'ToroidalVMMatrixDistribution');           
            testCase.verifyEqual(tvm.mu, tvmMle.mu, 'RelTol', 0.03)
            testCase.verifyEqual(tvm.kappa, tvmMle.kappa, 'RelTol', 0.2)
            testCase.verifyEqual(tvm.A, tvmMle.A, 'AbsTol', 0.3)
            
            tvmMle = ToroidalVMMatrixDistribution.mleNumerical(samples, 'trust-region');
            testCase.verifyClass(tvmMle, 'ToroidalVMMatrixDistribution');           
            testCase.verifyEqual(tvm.mu, tvmMle.mu, 'RelTol', 0.03)
            testCase.verifyEqual(tvm.kappa, tvmMle.kappa, 'RelTol', 0.2)
            testCase.verifyEqual(tvm.A, tvmMle.A, 'AbsTol', 0.3)
            
            tvmMle = ToroidalVMMatrixDistribution.mleNumerical(samples, 'trust-region-reflective');
            testCase.verifyClass(tvmMle, 'ToroidalVMMatrixDistribution');           
            testCase.verifyEqual(tvm.mu, tvmMle.mu, 'RelTol', 0.03)
            testCase.verifyEqual(tvm.kappa, tvmMle.kappa, 'RelTol', 0.2)
            testCase.verifyEqual(tvm.A, tvmMle.A, 'AbsTol', 0.3)
            
            tvmMle = ToroidalVMMatrixDistribution.mleNumerical(samples, 'interior-point');
            testCase.verifyClass(tvmMle, 'ToroidalVMMatrixDistribution');           
            testCase.verifyEqual(tvm.mu, tvmMle.mu, 'RelTol', 0.03)
            testCase.verifyEqual(tvm.kappa, tvmMle.kappa, 'RelTol', 0.2)
            testCase.verifyEqual(tvm.A, tvmMle.A, 'AbsTol', 0.3)
        end
        
        function testNormconstApprox(testCase)
            mu = [0, 0]';
            
            for k1 = [0.1 0.3 0.5 1 1.5 2]
                kappa = [k1, 0.3]';
                A = 0.1*[1 2; 3 4];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-3);
                
                kappa = [k1, k1]';
                A = 0.1*[1 2; 3 4];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-2);
                
                kappa = [0.5, 0.7]';
                A = 0.1*[k1 2; 3 4];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-5);
                
                kappa = [0.5, 0.7]';
                A = 0.1*[1 k1; 3 4];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-5);
                
                kappa = [0.5, 0.7]';
                A = 0.1*[1 2; k1 4];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-5);
                
                kappa = [0.5, 0.7]';
                A = 0.1*[1 2; 3 k1];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-5);
                
                kappa = [0.5, 0.7]';
                A = 0.1*[k1 k1; k1 k1];
                tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
                tvm.C = 1;
                testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 1E-5);
            end
            
            % worst case where approximation is currently used
            kappa = [1.5, 1.5]';
            A = [1 1; 1 1];
            tvm = ToroidalVMMatrixDistribution(mu, kappa, A);
            tvm.C = 1;
            testCase.verifyEqual(tvm.integral(), 1/tvm.normConstApprox, 'RelTol', 0.02);
        end
        
        function testShift(testCase)
            tvm = ToroidalVMMatrixDistribution([3;5],[0.7; 1.3], ones(2,2));
            s = [4;2];
            tvm2 = tvm.shift(s);
            testCase.verifyClass(tvm2, 'ToroidalVMMatrixDistribution');
            [xTest,yTest]=meshgrid(linspace(0,2*pi,20));
            testCase.verifyEqual(tvm2.pdf([xTest(:)';yTest(:)']),tvm.pdf([xTest(:)' - s(1);yTest(:)' - s(2)]),'AbsTol',1E-10);
        end           
    end
end
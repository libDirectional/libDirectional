classdef HypertoroidalGridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testApproxVMMixtureT2(testCase)
            dist = HypertoroidalMixture(...
                {ToroidalWNDistribution([1;1],[1,0.5;0.5,1]),...
                ToroidalWNDistribution([3;3],[1,-0.5;-0.5,1])},[0.5,0.5]);
            
            hgd = HypertoroidalGridDistribution.fromDistribution(dist, 31);
            testCase.verifyEqual(hgd.gridValues',dist.pdf(hgd.grid));
            testCase.verifyEqual(min(hgd.grid,[],2),[0;0]);
            testCase.verifyGreaterThan(max(hgd.grid,[],2),6);
            testCase.verifyEqual(hgd.gridType,'CartesianProd');
        end
        
        function testMultiplyTWNMixtureT2(testCase)
            rng(1)% Make it reproducible.
            % Multiplication has same mean as mixture of two
            % equally-weighted ones
            dist1 = ToroidalWNDistribution(HypertoroidalUniformDistribution(2).sample(1),[1,0.5;0.5,1]);
            dist2 = ToroidalWNDistribution(HypertoroidalUniformDistribution(2).sample(1),[1,-0.5;-0.5,1]);
            distMultUnnorm = CustomHypertoroidalDistribution(...
                @(x)dist1.pdf(x).*dist2.pdf(x),2);
            int=distMultUnnorm.integral;
            distMultNorm=CustomHypertoroidalDistribution(...
                @(x)dist1.pdf(x).*dist2.pdf(x)/int,2);
            hgd1 = HypertoroidalGridDistribution.fromDistribution(dist1, 101);
            hgd2 = HypertoroidalGridDistribution.fromDistribution(dist2, 101);
            hgdMult = hgd1.multiply(hgd2);
            
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,distMultNorm,hgdMult,1e-8);
            testCase.verifyEqual(hgdMult.gridType,'CartesianProd');
        end
        
        function testFromFunction3D(testCase)
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = false;end
            if enableExpensive
                testPoints = rand(3, 1000);
            else
                testPoints = rand(3, 30);
            end
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            mu = 2 * pi * rand(3, 1);
            hwnd = HypertoroidalWNDistribution(mu, C);
            nGridPoints = [27, 27, 27];
            hfdId = HypertoroidalGridDistribution.fromFunction( ...
                @(x)hwnd.pdf(x),nGridPoints,3, 'CartesianProd');
            testCase.verifyClass(hfdId, 'HypertoroidalGridDistribution')
            testCase.verifyEqual(hfdId.pdf(testPoints), hwnd.pdf(testPoints), 'AbsTol', 1E-6);
        end

        function testFromDistribution2D(testCase)
            % Test that from Distribution and fromFunction result in equal
            % approximations
            kappa1 = 0.3;
            kappa2 = 1.5;
            lambda = 0.5;
            coeffs = [5, 7];
            tvm = ToroidalVMSineDistribution([1; 2], [kappa1; kappa2], lambda);
            hgd1 = HypertoroidalGridDistribution.fromFunction(@(x)tvm.pdf(x), coeffs,2, 'CartesianProd');
            hgd2 = HypertoroidalGridDistribution.fromDistribution(tvm, coeffs, 'CartesianProd');
            testCase.verifyClass(hgd2, 'HypertoroidalGridDistribution');
            testCase.verifySize(hgd1.gridValues, [prod(coeffs),1]);
            % Verify approximation by validating coefficients
            testCase.verifyEqual(hgd2.gridValues, hgd1.gridValues, 'AbsTol', 1E-10);
        end
        function testConvolve2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            tvm1 = ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5);
            tvm2 = ToroidalVMSineDistribution([1; 4], [0.8; 1.5], 0.2);
            hfd1 = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'identity');
            hfdt = HypertoroidalFourierDistribution.fromDistribution(tvm2, [17, 15], 'identity');
            hgd1 = HypertoroidalGridDistribution.fromDistribution(tvm1, [17, 15]);
            hgdt = HypertoroidalGridDistribution.fromDistribution(tvm2, [17, 15]);
            
            hfd2 = hfd1.convolve(hfdt);
            hgd2 = hgd1.convolve(hgdt);
            testCase.verifyClass(hgd2, 'HypertoroidalGridDistribution')
            testCase.verifySize(hgd2.gridValues,size(hgd1.gridValues));
            testCase.verifyEqual(hgd2.grid,hgd1.grid);
            
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,hgd2,hfd2,1e-6);
        end
        
        function testShift(testCase)
            twn = ToroidalWNDistribution([1;1],[1,0.5;0.5,1]);
            hgd = HypertoroidalGridDistribution.fromDistribution(twn, [41,41]);
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,hgd,twn,1e-6);
            twnShift = twn.shift([1;2]);
            hgdShift = hgd.shift([1;2]);
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,hgdShift,twnShift,1e-6);
        end
    end
    methods (Static)
        function verifyPdfEqual(testCase,dist1,dist2,tol)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            [x, y] = meshgrid(linspace(0, 2*pi, 10), linspace(0, 2*pi, 10));
            fixture1 = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fixture2 = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(dist1.pdf([x(:)'; y(:)']), dist2.pdf([x(:)'; y(:)']), 'AbsTol', tol);
            fixture1.teardown;
            fixture2.teardown;
        end
    end
end

classdef HypertoroidalGridDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testGetGrid(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            grid = [1,2,3,4;1,2,3,4];
            hgd = HypertoroidalGridDistribution(grid,[1,1,1,1]/((2*pi)^2));
            testCase.verifyWarningFree(@()hgd.getGrid);
            testCase.verifyEqual(hgd.getGrid(),grid);
            hgd.grid = [];
            testCase.verifyError(@()hgd.getGrid,'Grid:UnknownGrid');
            hgd.gridType='CartesianProd';
            testCase.verifyWarning(@()hgd.getGrid,'Grid:GenerateDuringRunTime');
            testCase.applyFixture(SuppressedWarningsFixture('Grid:GenerateDuringRunTime'))
            testCase.verifySize(hgd.getGrid(),[1,4]);
        end
        
        function testApproxVMMixtureT2(testCase)
            dist = HypertoroidalMixture(...
                {ToroidalWNDistribution([1;1],[1,0.5;0.5,1]),...
                ToroidalWNDistribution([3;3],[1,-0.5;-0.5,1])},[0.5,0.5]);
            
            hgd = HypertoroidalGridDistribution.fromDistribution(dist, 31);
            testCase.verifyEqual(hgd.gridValues',dist.pdf(hgd.getGrid()));
            testCase.verifyEqual(min(hgd.getGrid(),[],2),[0;0]);
            testCase.verifyGreaterThan(max(hgd.getGrid(),[],2),6);
            testCase.verifyEqual(hgd.gridType,'CartesianProd');
        end
        
        function testFromDistributionFourier(testCase)
            for transformation = {'identity', 'sqrt'}
                dist = HypertoroidalMixture(...
                    {ToroidalWNDistribution([1;1],[1,0.5;0.5,1]),...
                    ToroidalWNDistribution([3;3],[1,-0.5;-0.5,1])},[0.5,0.5]);
                hfd = HypertoroidalFourierDistribution.fromDistribution(dist,[31,31],[transformation{:}]);
                hgdDirect = HypertoroidalGridDistribution.fromDistribution(hfd, [31,31]);
                hgdViaFunVals = HypertoroidalGridDistribution.fromDistribution(hfd, [31,32]);

                testCase.verifyEqual(hgdDirect.gridValues',hfd.pdf(hgdDirect.getGrid()),'RelTol',1e-13);
                testCase.verifyEqual(min(hgdDirect.getGrid(),[],2),[0;0]);
                testCase.verifyEqual(hgdDirect.gridType,'CartesianProd');

                testCase.verifyEqual(hgdViaFunVals.gridValues',hfd.pdf(hgdViaFunVals.getGrid()),'RelTol',1e-13);
                testCase.verifyEqual(min(hgdViaFunVals.getGrid(),[],2),[0;0]);
                testCase.verifyEqual(hgdViaFunVals.gridType,'CartesianProd');
            end
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
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
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
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            % Do not use function defined below since this case is 3-D
            testCase.verifyEqual(hfdId.pdf(testPoints), hwnd.pdf(testPoints), 'AbsTol', 1E-6);
        end

        function testFromDistribution2D(testCase)
            % Test that fromDistribution and fromFunction result in equal
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
            testCase.verifyEqual(hgd2.getGrid(),hgd1.getGrid());
            
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,hgd2,hfd2,1e-6);
        end
        
        function testPdfUnnormalized(testCase)
            twn = ToroidalWNDistribution([1;1],[1,0.5;0.5,1]);
            hgd = HypertoroidalGridDistribution.fromDistribution(twn, [41,41]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            chd = chd.scale(2);
            hgd.gridValues = 2*hgd.gridValues;
            
            chdFromHgd = CustomHypertoroidalDistribution(@(x)hgd.pdfUnnormalized(x),2);
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,chdFromHgd,chd,1e-6);
        end
        
        function testShift(testCase)
            twn = ToroidalWNDistribution([1;1],[1,0.5;0.5,1]);
            hgd = HypertoroidalGridDistribution.fromDistribution(twn, [41,41]);
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,hgd,twn,1e-6);
            twnShift = twn.shift([1;2]);
            hgdShift = hgd.shift([1;2]);
            HypertoroidalGridDistributionTest.verifyPdfEqual(testCase,hgdShift,twnShift,1e-6);
        end
        
        function testConditioning2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            grid = linspace(0,2*pi,100);
            z1 = 6;
            z2 = 5;
            cdNormDim1 = chd.conditionOn(1,z1);
            cdNormDim2 = chd.conditionOn(2,z2);
            
            hgd = HypertoroidalGridDistribution.fromDistribution(twn,[301,301]);
            hgdNormDim1 = hgd.conditionOn(1,z1);
            hgdNormDim2 = hgd.conditionOn(2,z2);
            testCase.verifyEqual(hgdNormDim1.integral(),1,'AbsTol',5e-16);
            testCase.verifyEqual(hgdNormDim2.integral(),1,'AbsTol',5e-16);
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(hgdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'RelTol',1e-12);
            testCase.verifyEqual(hgdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'RelTol',1e-12);
        end
        
        function testSliceAt2D(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            for z=[[6;5],zeros(2,1)] % Also test case 0 because it is handled differently
                z1 = z(1);
                z2 = z(2);
                cdSlicedDim1 = chd.sliceAt(1,z1);
                cdSlicedDim2 = chd.sliceAt(2,z2);
                for useFFTN = [true,false]
                    hgd = HypertoroidalGridDistribution.fromDistribution(twn,[301,301]);
                    hgdSlicedDim1 = hgd.sliceAt(1,z1,useFFTN);
                    hgdSlicedDim2 = hgd.sliceAt(2,z2,useFFTN);
                    % Cannot use pdf because it would normalize otherwise
                    testCase.verifyEqual(hgdSlicedDim1.gridValues',cdSlicedDim1.pdf(hgdSlicedDim1.getGrid()),'RelTol',5e-13);
                    testCase.verifyEqual(hgdSlicedDim2.gridValues',cdSlicedDim2.pdf(hgdSlicedDim2.getGrid()),'RelTol',5e-13);
                end
            end
        end
        
        function testConditioning3DTo2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            
            z1 = 1;
            z2 = 3;
            z3 = 4;
            cdNormDim1 = chd.conditionOn(1,z1);
            cdNormDim2 = chd.conditionOn(2,z2);
            cdNormDim3 = chd.conditionOn(3,z3);
            
            [mesh1,mesh2] = meshgrid(linspace(0,2*pi,20));
            
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[111,111,111]);
            hgdNormDim1 = hgd.conditionOn(1,z1);
            hgdNormDim2 = hgd.conditionOn(2,z2);
            hgdNormDim3 = hgd.conditionOn(3,z3);
            testCase.verifyEqual(hgdNormDim1.integral(),1,'AbsTol',5e-15);
            testCase.verifyEqual(hgdNormDim2.integral(),1,'AbsTol',5e-15);
            testCase.verifyEqual(hgdNormDim3.integral(),1,'AbsTol',5e-15);
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(hgdNormDim1.pdf([mesh1(:)';mesh2(:)']),cdNormDim1.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
            testCase.verifyEqual(hgdNormDim2.pdf([mesh1(:)';mesh2(:)']),cdNormDim2.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
            testCase.verifyEqual(hgdNormDim3.pdf([mesh1(:)';mesh2(:)']),cdNormDim3.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
        end
        
        function testSliceAt3DTo2D(testCase)
            C = 5*[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            
            z1 = 1;
            z2 = 3;
            z3 = 4;
            cdSlicedDim1 = chd.sliceAt(1,z1);
            cdSlicedDim2 = chd.sliceAt(2,z2);
            cdSlicedDim3 = chd.sliceAt(3,z3);
            
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[41,41,41]);
            for useFFTN = [true,false]
                hgdSlicedDim1 = hgd.sliceAt(1,z1,useFFTN);
                hgdSlicedDim2 = hgd.sliceAt(2,z2,useFFTN);
                hgdSlicedDim3 = hgd.sliceAt(3,z3,useFFTN);

                testCase.verifyEqual(hgdSlicedDim1.gridValues',cdSlicedDim1.pdf(hgdSlicedDim1.getGrid()),'RelTol',5e-11);
                testCase.verifyEqual(hgdSlicedDim2.gridValues',cdSlicedDim2.pdf(hgdSlicedDim2.getGrid()),'RelTol',5e-11);
                testCase.verifyEqual(hgdSlicedDim3.gridValues',cdSlicedDim3.pdf(hgdSlicedDim3.getGrid()),'RelTol',5e-11);
            end
        end
        
        function testConditioning3DTo1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            
            z1 = [1;5];
            z2 = [2;4];
            z3 = [3;6];
            cdNormDim1 = chd.conditionOn([1,2],z1);
            cdNormDim2 = chd.conditionOn([1,3],z2);
            cdNormDim3 = chd.conditionOn([2,3],z3);
            
            grid = linspace(0,2*pi,100);
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[111,111,111]);
            hgdNormDim1 = hgd.conditionOn([1,2],z1);
            hgdNormDim2 = hgd.conditionOn([1,3],z2);
            hgdNormDim3 = hgd.conditionOn([2,3],z3);
            testCase.verifyEqual(hgdNormDim1.integral(),1,'AbsTol',1e-15);
            testCase.verifyEqual(hgdNormDim2.integral(),1,'AbsTol',1e-15);
            testCase.verifyEqual(hgdNormDim3.integral(),1,'AbsTol',1e-15);
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(hgdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'AbsTol',1e-9);
            testCase.verifyEqual(hgdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'AbsTol',1e-9);
            testCase.verifyEqual(hgdNormDim3.pdf(grid),cdNormDim3.pdf(grid),'AbsTol',1e-9);
        end
        
        function testSliceAt3DTo1D(testCase)
            C = 5*[0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            
            z1 = [1;5];
            z2 = [2;4];
            z3 = [3;6];
            cdSlicedDim1 = chd.sliceAt([1,2],z1);
            cdSlicedDim2 = chd.sliceAt([1,3],z2);
            cdSlicedDim3 = chd.sliceAt([2,3],z3);
            
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[41,41,41]);
            for useFFTN = [true,false]
                hgdSlicedDim1 = hgd.sliceAt([1,2],z1,useFFTN);
                hgdSlicedDim2 = hgd.sliceAt([1,3],z2,useFFTN);
                hgdSlicedDim3 = hgd.sliceAt([2,3],z3,useFFTN);
                testCase.verifyEqual(hgdSlicedDim1.gridValues',cdSlicedDim1.pdf(hgdSlicedDim1.getGrid()),'RelTol',5e-11);
                testCase.verifyEqual(hgdSlicedDim2.gridValues',cdSlicedDim2.pdf(hgdSlicedDim2.getGrid()),'RelTol',5e-11);
                testCase.verifyEqual(hgdSlicedDim3.gridValues',cdSlicedDim3.pdf(hgdSlicedDim3.getGrid()),'RelTol',5e-11);
            end
        end
        
        function testMarginalize2DTo1D(testCase)        
            twn = ToroidalWNDistribution([3;4],2*[1,0.8;0.8,1]);
            hgd = HypertoroidalGridDistribution.fromDistribution(twn,[31,31]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            
            gdMarg = hgd.marginalizeOut(1);
            cdMarg = chd.marginalizeOut(1);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',1e-15);
            
            gdMarg = hgd.marginalizeOut(2);
            cdMarg = chd.marginalizeOut(2);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',1e-15);
        end
        
        function testMarginalize3Dto2D(testCase)
            hwn = HypertoroidalWNDistribution([3;4;6],2*[1,0.8,0.3;0.8,1,0.5;0.3,0.5,2]);
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn, [13, 13, 13]);
            chd = CustomHypertoroidalDistribution.fromDistribution(hwn);
            
            gdMarg = hgd.marginalizeOut(1);
            cdMarg = chd.marginalizeOut(1);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',1e-15);
            
            gdMarg = hgd.marginalizeOut(2);
            cdMarg = chd.marginalizeOut(2);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',1e-15);
            
            gdMarg = hgd.marginalizeOut(3);
            cdMarg = chd.marginalizeOut(3);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',1e-15);
        end
        
        function testMarginalize3Dto1D(testCase)
            twn = HypertoroidalWNDistribution([3;4;6],2*[1,0.8,0.3;0.8,1,0.5;0.3,0.5,2]);
            hgd = HypertoroidalGridDistribution.fromDistribution(twn, [31, 31, 31]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            
            gdMarg = hgd.marginalizeOut([2,3]);
            cdMarg = chd.marginalizeOut([2,3]);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',5e-7);
            
            gdMarg = hgd.marginalizeOut([1,3]);
            cdMarg = chd.marginalizeOut([1,3]);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',5e-7);
            
            gdMarg = hgd.marginalizeOut([1,2]);
            cdMarg = chd.marginalizeOut([1,2]);
            testCase.verifyEqual(gdMarg.gridValues',cdMarg.pdf(gdMarg.getGrid()),'RelTol',5e-7);
        end
        
        function testLikelihood2DTo1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            hgd = HypertoroidalGridDistribution.fromDistribution(twn,[101,101]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            xFix = 1;
            yFix = 3;
            
            likelihoodDim1hgd = hgd.likelihood(1,xFix);
            likelihoodDim1chd = chd.likelihood(1,xFix);
            likelihoodDim2hgd = hgd.likelihood(2,yFix);
            likelihoodDim2chd = chd.likelihood(2,yFix);
            
            grid = linspace(-pi,3*pi,100);
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(likelihoodDim1hgd.pdfUnnormalized(grid),likelihoodDim1chd.pdf(grid),'RelTol',5e-6);
            testCase.verifyEqual(likelihoodDim2hgd.pdfUnnormalized(grid),likelihoodDim2chd.pdf(grid),'RelTol',5e-6);
        end
        
        function testLikelihood3DTo2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[101,101,101]);
            chd = CustomHypertoroidalDistribution.fromDistribution(hwn);
            xFix = 1;
            yFix = 3;
            zFix = 4;
            
            likelihoodDim1hgd = hgd.likelihood(1,xFix);
            likelihoodDim1chd = chd.likelihood(1,xFix);
            likelihoodDim2hgd = hgd.likelihood(2,yFix);
            likelihoodDim2chd = chd.likelihood(2,yFix);
            likelihoodDim3hgd = hgd.likelihood(3,zFix);
            likelihoodDim3chd = chd.likelihood(3,zFix);
            
            [mesh1,mesh2] = meshgrid(linspace(0,2*pi,20));
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(likelihoodDim1hgd.pdfUnnormalized([mesh1(:)';mesh2(:)']),likelihoodDim1chd.pdf([mesh1(:)';mesh2(:)']),'RelTol',5e-5);
            testCase.verifyEqual(likelihoodDim2hgd.pdfUnnormalized([mesh1(:)';mesh2(:)']),likelihoodDim2chd.pdf([mesh1(:)';mesh2(:)']),'RelTol',5e-5);
            testCase.verifyEqual(likelihoodDim3hgd.pdfUnnormalized([mesh1(:)';mesh2(:)']),likelihoodDim3chd.pdf([mesh1(:)';mesh2(:)']),'RelTol',5e-5);
        end
        
        function testLikelihood3DTo1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[101,101,101]);
            chd = CustomHypertoroidalDistribution.fromDistribution(hwn);
            xyFix = [1;5];
            xzFix = [2;4];
            yzFix = [3;6];
            
            likelihoodDim12chd = chd.likelihood([1,2],xyFix);
            likelihoodDim12hgd = hgd.likelihood([1,2],xyFix);
            likelihoodDim13chd = chd.likelihood([1,3],xzFix);
            likelihoodDim13hgd = hgd.likelihood([1,3],xzFix);
            likelihoodDim23chd = chd.likelihood([2,3],yzFix);
            likelihoodDim23hgd = hgd.likelihood([2,3],yzFix);
            
            grid = linspace(-pi,3*pi,100);
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(likelihoodDim12hgd.pdfUnnormalized(grid),likelihoodDim12chd.pdf(grid),'RelTol',5e-5);
            testCase.verifyEqual(likelihoodDim13hgd.pdfUnnormalized(grid),likelihoodDim13chd.pdf(grid),'RelTol',5e-5);
            testCase.verifyEqual(likelihoodDim23hgd.pdfUnnormalized(grid),likelihoodDim23chd.pdf(grid),'RelTol',5e-5);
        end
        
        function testPlotting(testCase)
            grid = [1,2,3,4;1,2,3,4];
            hgd = HypertoroidalGridDistribution(grid,[1,1,1,1]/((2*pi)^2));
            testCase.verifyWarningFree(@()hgd.plot());
        end
        function testPlotInterpolated(testCase)
            hgd = HypertoroidalGridDistribution.fromDistribution(ToroidalWNDistribution([1;2],[1,0.3;0.3,1]),[11,11]);
            testCase.verifyWarningFree(@()hgd.plot());
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

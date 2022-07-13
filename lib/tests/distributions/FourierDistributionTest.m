classdef FourierDistributionTest < matlab.unittest.TestCase
    properties
    end
    
    methods(Static)
        function testFourierConversion(testCase, dist, coeffs, transformation, tolerance)
            xvals = -2 * pi:0.01:3 * pi;
            fd = FourierDistribution.fromDistribution(dist, coeffs, transformation);
            testCase.verifyEqual(length(fd.a)+length(fd.b), coeffs);
            testCase.verifyEqual(fd.pdf(xvals), dist.pdf(xvals), 'AbsTol', tolerance);
        end
    end
    
    methods(Test)
        % Test conversions
        function testVMToFourierId(testCase)
            mu = 0.4;
            for kappa = .1:.1:2
                dist = VMDistribution(mu, kappa);
                FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-8);
            end
        end
        function testVMToFourierSqrt(testCase)
            mu = 0.5;
            for kappa = .1:.1:2
                dist = VMDistribution(mu, kappa);
                FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-8);
            end
        end
        function testWNToFourierId(testCase)
            mu = 0.8;
            for sigma = .2:.1:2
                dist = WNDistribution(mu, sigma);
                FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-8);
            end
        end
        function testWNToFourierSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            mu = 0.9;
            for sigma = .2:.1:2
                dist = WNDistribution(mu, sigma);
                FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-8);
            end
        end
        function testWCToFourierId(testCase)
            mu = 1.2;
            for gamma = .8:.1:3
                dist = WCDistribution(mu, gamma);
                FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-7);
            end
        end
        function testWCToFourierSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:ApproximationHypergeometric'));
            mu = 1.3;
            for gamma = .8:.1:3
                dist = WCDistribution(mu, gamma);
                FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-7);
            end
        end
        function testWEToFourierId(testCase)
            %Treat differently due to lack of continuity
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            for lambda = .1:.1:2
                xvals = -2 * pi:0.01:3 * pi;
                xvals = xvals(mod(xvals, 2*pi) > 0.5 & mod(xvals, 2*pi) < (2 * pi - 0.5));
                dist = WEDistribution(lambda);
                fd = FourierDistribution.fromDistribution(dist, 1001, 'identity');
                testCase.verifyEqual(length(fd.a)+length(fd.b), 1001);
                testCase.verifyEqual(fd.pdf(xvals), dist.pdf(xvals), 'AbsTol', 5E-3);
            end
        end
        function testWEToFourierSqrt(testCase)
            %For sqrt, same applies as to identity
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            for lambda = .1:.1:2
                xvals = -2 * pi:0.01:3 * pi;
                xvals = xvals(mod(xvals, 2*pi) > 0.5 & mod(xvals, 2*pi) < (2 * pi - 0.5));
                dist = WEDistribution(lambda);
                fd = FourierDistribution.fromDistribution(dist, 1001, 'sqrt');
                testCase.verifyEqual(length(fd.a)+length(fd.b), 1001);
                testCase.verifyEqual(fd.pdf(xvals), dist.pdf(xvals), 'AbsTol', 5E-3);
            end
        end
        function testWLToFourierId(testCase)
            %Only test parameter combinations that don't result in too
            %abrupt changes in pdf values
            for lambda = 0.1:0.2:1
                for kappa = 0.1:0.5:4
                    dist = WLDistribution(lambda, kappa);
                    FourierDistributionTest.testFourierConversion(testCase, dist, 1001, 'identity', 1E-3);
                end
            end
        end
        function testWLToFourierSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            for lambda = 0.1:0.2:1
                for kappa = 0.1:0.5:4
                    dist = WLDistribution(lambda, kappa);
                    FourierDistributionTest.testFourierConversion(testCase, dist, 1001, 'sqrt', 1E-3);
                end
            end
        end
        function testWDToFourierMoments(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:WDMatching'));
            wd = WDDistribution([1,2,3],[1,1,1]/3);
            fdId = FourierDistribution.fromDistribution(wd,19,'identity');
            fdSqrt = FourierDistribution.fromDistribution(wd,50001,'sqrt');
            
            highestMoment = 9;
            momentsWd = arrayfun(@(i)wd.trigonometricMoment(i), -highestMoment:highestMoment);
            momentsFdId = arrayfun(@(i)fdId.trigonometricMoment(i), -highestMoment:highestMoment);
            momentsFdSqrt = arrayfun(@(i)fdSqrt.trigonometricMoment(i), -highestMoment:highestMoment);
            
            testCase.verifyEqual(momentsFdId,momentsWd,'AbsTol',1e-15);
            % It is much harder for sqrt since it does not merely save them
            % but has to approximate them by convoling the coefficient
            % vector with itself
            testCase.verifyEqual(momentsFdSqrt,momentsWd,'AbsTol',2e-4);
        end
        function testWDToFourierMomentsNonequalWeights(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            wd = WDDistribution([1,2,3],[1,2,2]/5);
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:WDMatching'));
            fdId = FourierDistribution.fromDistribution(wd,19,'identity');
            fdSqrt = FourierDistribution.fromDistribution(wd,50001,'sqrt');
            
            highestMomentToCompare = 9;
            momentsWd = arrayfun(@(i)wd.trigonometricMoment(i), -highestMomentToCompare:highestMomentToCompare);
            momentsFdId = arrayfun(@(i)fdId.trigonometricMoment(i), -highestMomentToCompare:highestMomentToCompare);
            momentsFdSqrt = arrayfun(@(i)fdSqrt.trigonometricMoment(i), -highestMomentToCompare:highestMomentToCompare);
            
            testCase.verifyEqual(momentsFdId,momentsWd,'AbsTol',1e-15)
            % It is much harder for sqrt since it does not merely save them
            % but has to approximate them by convoling the coefficient
            % vector with itself
            testCase.verifyEqual(momentsFdSqrt,momentsWd,'AbsTol',2e-4)
        end
        function testWDToFourierCdf(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            wd = WDDistribution([1,2,3],[1,1,1]/3);
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:WDMatching'));
            fdId = FourierDistribution.fromDistribution(wd,10001,'identity');
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fdSqrt = FourierDistribution.fromDistribution(wd,10001,'sqrt');
            
            evalAt = linspace(0,2*pi-0.01,50);
            testCase.verifyEqual(fdId.cdf(evalAt),wd.cdf(evalAt),'AbsTol',1e-3)
            testCase.verifyEqual(fdSqrt.cdf(evalAt),wd.cdf(evalAt),'AbsTol',1e-3)
        end
        function testWDToFourierCdfNonequalWeights(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            wd = WDDistribution([1,2,3],[1,2,2]/5);
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:WDMatching'));
            fdId = FourierDistribution.fromDistribution(wd,10001,'identity');
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fdSqrt = FourierDistribution.fromDistribution(wd,10001,'sqrt');
            
            evalAt = linspace(0,2*pi-0.01,50);
            testCase.verifyEqual(fdId.cdf(evalAt),wd.cdf(evalAt),'AbsTol',1e-3)
            testCase.verifyEqual(fdSqrt.cdf(evalAt),wd.cdf(evalAt),'AbsTol',1e-3)
        end
        function testCircularUniformToFourierId(testCase)
            dist = CircularUniformDistribution();
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-8);
        end
        function testCircularUniformToFourierSqrt(testCase)
            dist = CircularUniformDistribution();
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-8);
        end
        function testGvMToFourierLog(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist1 = GvMDistribution([1; 2; 2.5], [2; 3; 4]);
            dist2 = GvMDistribution((1:50)', linspace(2, 10, 50)');
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            FourierDistributionTest.testFourierConversion(testCase, dist1, 101, 'log', 1E-8);
            testCase.verifyWarning(@()FourierDistribution.fromDistribution(dist2, 99, 'log'), 'Conversion:GvMTooFewCoefficients');
        end
        function testGCMToFourierId(testCase)
            vm = VMDistribution(1, 2);
            wn = WNDistribution(2, 1);
            dist = CircularMixture({vm, wn}, [0.3, 0.7]);
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-8);
        end
        function testGCMToFourierSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            vm = VMDistribution(1, 2);
            wn = WNDistribution(2, 1);
            dist = CircularMixture({vm, wn}, [0.3, 0.7]);
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-8);
        end
        function testCCDToFourierId(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormula'));
            vm = VMDistribution(1, 2);
            dist = CustomCircularDistribution(@(x)vm.pdf(x));
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-8);
        end
        function testCCDToFourierSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormula'));
            vm = VMDistribution(1, 2);
            dist = CustomCircularDistribution(@(x)vm.pdf(x));
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-8);
        end
        function testFIGToFourier(testCase)
            vm = VMDistribution(1,2);
            dist = FIGDistribution.fromDistribution(vm,101);
            
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-15);
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-15);
            
            dist.enforcePdfNonnegative = false;
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'identity', 1E-15);
            FourierDistributionTest.testFourierConversion(testCase, dist, 101, 'sqrt', 1E-15);
            
        end
        % Test function evaluation
        function testValueVsIFFT(testCase)
            vm = VMDistribution(1, 2);
            noCoeffs = 51;
            xVals = linspace(0, 2*pi, noCoeffs+1);
            for transformation = {'identity', 'sqrt'}
                fd = FourierDistribution.fromDistribution(vm, noCoeffs, [transformation{:}]);
                pdfValsIFFT = ifft(ifftshift(fd.c), 'symmetric') * noCoeffs;
                testCase.verifyEqual([pdfValsIFFT, pdfValsIFFT(1)], fd.value(xVals), 'AbsTol', 1E-14);
            end
        end
        % Test coefficient conversions
        function testCoefficientConversion1(testCase)
            a = [1 / pi, 4, 3, 2, 1];
            b = [4, 3, 2, 1];
            fd1 = FourierDistribution(a, b, 'identity');
            fd2 = FourierDistribution.fromComplex(fd1.c, 'identity');
            testCase.verifyEqual(fd2.a, fd1.a, 'AbsTol', 1E-10)
            testCase.verifyEqual(fd2.b, fd1.b, 'AbsTol', 1E-10)
        end
        function testCoefficientConversion2(testCase)
            kappa = 2;
            fd1 = FourierDistribution.fromDistribution(VMDistribution(0, kappa), 101, 'identity');
            fd2 = FourierDistribution.fromComplex(fd1.c, 'identity');
            testCase.verifyEqual(fd2.a, fd1.a, 'AbsTol', 1E-10)
            testCase.verifyEqual(fd2.b, fd1.b, 'AbsTol', 1E-10)
        end
        function testFromFunction(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for kappa = 0.1:0.3:4
                vm = VMDistribution(3, kappa);
                fd = FourierDistribution.fromFunction(@(x)vm.pdf(x), 101, 'sqrt');
                testCase.verifyEqual(fd.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        function testFromFunctionValues(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            vm = VMDistribution(3, 1);
            fvals = vm.pdf(linspace(0, 2*pi, 100));
            fvals(end) = [];
            fd1 = FourierDistribution.fromFunctionValues(fvals, 99, 'sqrt');
            fvals = vm.pdf(linspace(0, 2*pi, 101));
            fvals(end) = [];
            % Test truncation in fromFunction
            fd2 = FourierDistribution.fromFunctionValues(fvals, 99, 'sqrt');
            testCase.verifyEqual(fd1.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            testCase.verifyEqual(fd2.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
        end
        % Test operations for prediction and filter steps
        function testMultiplyVM(testCase)
            for transformation={'identity', 'sqrt'}
                currTrans=[transformation{:}];
                xvals = -2 * pi:0.01:3 * pi;
                for kappa1 = 0.1:0.3:4
                    for kappa2 = 0.1:0.3:4
                        dist1 = VMDistribution(0, kappa1);
                        dist2 = VMDistribution(0, kappa2);
                        f1 = FourierDistribution.fromDistribution(dist1, 101, currTrans);
                        f2 = FourierDistribution.fromDistribution(dist2, 101, currTrans);
                        fFiltered = f1.multiply(f2);
                        distResult = dist1.multiply(dist2);
                        testCase.verifyEqual(fFiltered.pdf(xvals), distResult.pdf(xvals), 'AbsTol', 1E-8);
                    end
                end
            end
        end
        function testConvolveWN(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            xvals = -2 * pi:0.01:3 * pi;
            for transformation={'identity', 'sqrt'}
                currTrans=[transformation{:}];
                for sigma1 = 0.1:0.3:4
                    for sigma2 = 0.1:0.3:4
                        dist1 = WNDistribution(0, sigma1);
                        dist2 = WNDistribution(0, sigma2);
                        f1 = FourierDistribution.fromDistribution(dist1, 101, currTrans);
                        f2 = FourierDistribution.fromDistribution(dist2, 101, currTrans);
                        fPredicted = f1.convolve(f2, 101);
                        distResult = dist1.convolve(dist2);
                        testCase.verifyEqual(fPredicted.pdf(xvals), distResult.pdf(xvals), 'AbsTol', 1E-8);
                    end
                end
            end
        end
        function testMomentsFourierRootPredictionVM(testCase)
            for kappa1 = 0.1:0.3:4
                for kappa2 = 0.1:0.3:4
                    vm1 = VMDistribution(0, kappa1);
                    vm2 = VMDistribution(0, kappa2);
                    vmRes = vm1.convolve(vm2);
                    f1 = FourierDistribution.fromDistribution(vm1, 101, 'sqrt');
                    f2 = FourierDistribution.fromDistribution(vm2, 101, 'sqrt');
                    fPredicted = f1.convolve(f2, 101);
                    testCase.verifyEqual(fPredicted.trigonometricMoment(1), vmRes.trigonometricMoment(1), 'AbsTol', 1E-8);
                end
            end
        end
        function testErrorsPredictAndFilter(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            fd1 = FourierDistribution.fromDistribution(VMDistribution(0, 1), 101, 'sqrt');
            fd2 = FourierDistribution.fromDistribution(VMDistribution(0, 1), 101, 'identity');
            testCase.verifyError(@()fd1.multiply(fd2), 'Multiply:differentTransformations');
            testCase.verifyError(@()fd1.convolve(fd2), 'Convolve:differentTransformations');
            fd3 = fd2.transformViaCoefficients('square', 101);
            testCase.verifyError(@()fd3.multiply(fd3), 'Multiply:unsupportedTransformation');
            testCase.verifyError(@()fd3.convolve(fd3), 'Convolve:unsupportedTransformation');
        end
        % Test transformations
        function testTransformViaFFT(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for kappa = 0.1:0.3:4
                vm = VMDistribution(3, kappa);
                fd1 = FourierDistribution.fromDistribution(vm, 101, 'identity');
                fd2 = fd1.transformViaFFT('sqrt', 101);
                testCase.verifyEqual(fd2.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        function testTransformViaVM(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for sigma = 0.1:0.1:4
                wntmp = WNDistribution(1, sigma);
                fd = FourierDistribution.fromDistribution(wntmp, 101, 'identity');
                vm = wntmp.toVM;
                fdtrans = fd.transformViaVM('sqrt', 101);
                testCase.verifyEqual(fdtrans.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        function testSquaring(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for kappa = 0.1:0.3:4
                dist = VMDistribution(0, kappa);
                fd1 = FourierDistribution.fromDistribution(dist, 101, 'sqrt');
                fd2 = fd1.transformViaCoefficients('square', 101);
                testCase.verifyEqual(fd2.pdf(xvals), dist.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        function testErrorRepeatedTransformation(testCase)
            vm = VMDistribution(3, 1);
            fd1 = FourierDistribution.fromFunction(@(x)vm.pdf(x), 9, 'identity');
            fd2 = fd1.transformViaFFT('sqrt', 9);
            testCase.verifyError(@()fd2.transformViaFFT('sqrt', 9), 'Transformation:cannotCombine')
        end
        function testSquareAfterSqrt(testCase)
            vm = VMDistribution(3, 1);
            fd1 = FourierDistribution.fromFunction(@(x)vm.pdf(x), 101, 'identity');
            fd2 = fd1.transformViaFFT('sqrt', 101);
            fd3 = fd2.transformViaCoefficients('square', 101);
            testCase.verifyEqual(fd3.transformation, 'identity');
            testCase.verifyEqual(fd3.a, fd1.a, 'AbsTol', 1E-8);
            testCase.verifyEqual(fd3.b, fd1.b, 'AbsTol', 1E-8);
        end
        function testCustomTransformation(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            xvals = -2 * pi:0.01:3 * pi;
            vm = VMDistribution(3, 1);
            fvals = vm.pdf(linspace(0, 2*pi, 100)).^3;
            fvals(end) = [];
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            fd = FourierDistribution.fromFunctionValues(fvals, 99, 'custom');
            testCase.verifyEqual(fd.value(xvals), vm.pdf(xvals).^3, 'AbsTol', 1E-8);
        end
        function testTransformToMoreCoeffs(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(3, 1);
            fdSqrt = FourierDistribution.fromDistribution(vm, 1001, 'sqrt');
            fdId = FourierDistribution.fromDistribution(vm, 101, 'identity');
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            fdIdSqrt = fdId.transformViaFFT('sqrt', 1001);
            testCase.verifyEqual(length(fdIdSqrt.a)+length(fdIdSqrt.b), 1001);
            testCase.verifyEqual(fdIdSqrt.a, fdSqrt.a, 'AbsTol', 1E-8);
            testCase.verifyEqual(fdIdSqrt.b, fdSqrt.b, 'AbsTol', 1E-8);
        end
        function testTransformToLessCoeffs(testCase)
            %case is easy enough to work with less coeffs
            vm = VMDistribution(3, 1);
            fdSqrt = FourierDistribution.fromDistribution(vm, 101, 'sqrt');
            fdId = FourierDistribution.fromDistribution(vm, 1001, 'identity');
            fdIdSqrt = fdId.transformViaFFT('sqrt', 101);
            testCase.verifyEqual(length(fdIdSqrt.a)+length(fdIdSqrt.b), 101, 'AbsTol', 1E-8);
            testCase.verifyEqual(fdIdSqrt.a, fdSqrt.a, 'AbsTol', 1E-8);
            testCase.verifyEqual(fdIdSqrt.b, fdSqrt.b, 'AbsTol', 1E-8);
        end
        % CDF tests
        function testCdfIdentityStartZero(testCase)
            fd = FourierDistribution.fromDistribution(VMDistribution(1, 3), 101, 'identity');
            xvals = 0:0.01:2 * pi;
            intValsNumerically = arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend), xvals);
            testCase.verifyEqual(fd.cdf(xvals, 0), intValsNumerically, 'AbsTol', 1E-8);
        end
        function testCdfSqrtStartZero(testCase)
            fd = FourierDistribution.fromDistribution(VMDistribution(1, 3), 101, 'sqrt');
            xvals = 0:0.01:2 * pi;
            intValsNumerically = arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend), xvals);
            testCase.verifyEqual(fd.cdf(xvals, 0), intValsNumerically, 'AbsTol', 1E-8);
        end
        function testCdfIdentityStartNonzero(testCase)
            fd = FourierDistribution.fromDistribution(VMDistribution(1, 3), 101, 'identity');
            xvals = 0:0.01:2 * pi;
            startingPoint = 1;
            intValsNumerically = [arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend)+integral(@(x)fd.pdf(x), startingPoint, 2*pi), xvals(xvals < startingPoint)), ...
                arrayfun(@(xend)integral(@(x)fd.pdf(x), startingPoint, xend), xvals(xvals >= startingPoint))];
            testCase.verifyEqual(fd.cdf(xvals, startingPoint), intValsNumerically, 'AbsTol', 1E-8);
        end
        %other tests
        function fromComplexTestId(testCase)
            fdOrig = FourierDistribution.fromDistribution(VMDistribution(0.5,1),13,'identity');
            fdConvertedTwice = FourierDistribution.fromComplex(fdOrig.c,'identity');
            testCase.verifyEqual(fdConvertedTwice.a,fdOrig.a);
            testCase.verifyEqual(fdConvertedTwice.b,fdOrig.b);
        end
        function fromComplexTestSqrt(testCase)
            fdOrig = FourierDistribution.fromDistribution(VMDistribution(0.5,1),13,'sqrt');
            fdConvertedTwice = FourierDistribution.fromComplex(fdOrig.c,'sqrt');
            testCase.verifyEqual(fdConvertedTwice.a,fdOrig.a);
            testCase.verifyEqual(fdConvertedTwice.b,fdOrig.b);
        end
        function testTruncation(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(3, 1);
            fd1 = FourierDistribution.fromDistribution(vm, 101, 'sqrt');
            testCase.verifyWarning(@()fd1.truncate(1001), 'Truncate:TooFewCoefficients');
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            fd2 = fd1.truncate(1001);
            testCase.verifyEqual(length(fd2.a)+length(fd2.b), 1001, 'AbsTol', 1E-8);
            fd3 = fd2.truncate(51);
            testCase.verifyEqual(length(fd3.a)+length(fd3.b), 51, 'AbsTol', 1E-8);
        end
        function testMoments(testCase)
            vm = VMDistribution(2, 1);
            fd1 = FourierDistribution.fromDistribution(vm, 15, 'identity');
            fd2 = FourierDistribution.fromDistribution(vm, 15, 'sqrt');
            for i = -2:3
                vmMoment = vm.trigonometricMoment(i);
                testCase.verifyEqual(fd1.trigonometricMoment(i), vmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(fd1.trigonometricMomentNumerical(i), vmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(fd2.trigonometricMoment(i), vmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(fd2.trigonometricMomentNumerical(i), vmMoment, 'AbsTol', 1E-6);
            end
            testCase.verifyEqual(fd1.trigonometricMoment(16), 0);
            testCase.verifyEqual(fd2.trigonometricMoment(30), 0);
        end
        function testNormalizationAfterTruncation(testCase)
            dist = VMDistribution(0, 5);
            fdSqrt = FourierDistribution.fromDistribution(dist, 1001, 'sqrt');
            testCase.verifyEqual(integral(@(x)fdSqrt.pdf(x), 0, 2*pi), 1, 'RelTol', 1E-4);
            % Perform truncation without normalization
            fdSqrtManuallyTruncated = fdSqrt;
            fdSqrtManuallyTruncated.a = fdSqrt.a(1:3);
            fdSqrtManuallyTruncated.b = fdSqrt.b(1:2);
            % Verify that this case requires normalization after truncation
            testCase.verifyTrue(abs(integral(@(x)fdSqrtManuallyTruncated.pdf(x), 0, 2*pi)-1) >= 0.01);
            fdSqrtTrunc = fdSqrt.truncate(7);
            % Verify that use of .truncate results in a normalized density
            testCase.verifyEqual(integral(@(x)fdSqrtTrunc.pdf(x), 0, 2*pi), 1, 'RelTol', 1E-4);
        end
        function testNormalizationLog(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            fdLog1 = FourierDistribution.fromDistribution(VMDistribution(1, 2), 21, 'log');
            fdLog2 = FourierDistribution.fromDistribution(VMDistribution(2, 1), 21, 'log');
            testCase.applyFixture(SuppressedWarningsFixture('Multiply:NotNormalizing'));
            fdMult = fdLog1.multiply(fdLog2);
            testCase.verifyEqual(fdMult.integralNumerical(0, 2*pi), 1, 'AbsTol', 1E-4);
            fdMultTrunc = fdMult.truncate(5);
            testCase.verifyEqual(fdMultTrunc.integralNumerical(0, 2*pi), 1, 'AbsTol', 1E-4);
        end
        function testIntegral(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            dist = VMDistribution(0, 5);
            for transformation = {'identity', 'sqrt', 'log'}
                fd = FourierDistribution.fromDistribution(dist, 15, [transformation{:}]);
                testCase.verifyEqual(fd.integral(0, 1.5), fd.integralNumerical(0, 1.5), 'RelTol', 1E-8);
                testCase.verifyEqual(fd.integral(1.5, 0), fd.integralNumerical(1.5, 0), 'RelTol', 1E-8);
                testCase.verifyEqual(fd.integral(10, -10), fd.integralNumerical(10, -10), 'RelTol', 1E-8);
            end
        end
        function testPdfOnGrid(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            for transformation = {'identity', 'sqrt', 'log'}
                dist = VMDistribution(0, 5);
                fd = FourierDistribution.fromDistribution(dist, 15, [transformation{:}]);
                [vals, xgrid] = fd.pdfOnGrid(99);
                testCase.verifyEqual(vals, fd.pdf(xgrid), 'RelTol', 1E-8);
            end
        end
        function testHellingerDistance(testCase)
            dist1 = VMDistribution(0, 5);
            dist2 = VMDistribution(0, 3);
            fd1 = FourierDistribution.fromDistribution(dist1, 51, 'sqrt');
            fd2 = FourierDistribution.fromDistribution(dist2, 51, 'sqrt');
            testCase.verifyEqual(fd1.hellingerDistanceNumerical(fd2), dist1.hellingerDistanceNumerical(dist2), 'RelTol', 1E-8);
            testCase.verifyEqual(fd1.hellingerDistanceNumerical(dist2), dist1.hellingerDistanceNumerical(dist2), 'RelTol', 1E-8);
        end
        function testTransformViaFFTForTransformed(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist = VMDistribution(2, 5);
            fdSqrt = FourierDistribution.fromDistribution(dist, 11, 'sqrt');
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            fdLog = FourierDistribution.fromDistribution(dist, 11, 'log');
            fdIdFromSqrt = fdSqrt.transformViaFFT('square');
            fdIdFromLog = fdLog.transformViaFFT('power');
            fdIdFromSqrt101 = fdSqrt.transformViaFFT('square', 511);
            fdIdFromLog101 = fdLog.transformViaFFT('power', 51);
            xvals = 0:0.01:2 * pi;
            
            % Test that approximation is bad when using few coefficients
            testCase.verifyGreaterThan(max(fdSqrt.pdf(xvals)-fdIdFromSqrt.pdf(xvals)), 0.01);
            testCase.verifyGreaterThan(max(fdLog.pdf(xvals)-fdIdFromLog.pdf(xvals)), 0.01);
            % Test that approximation is good when using many coefficients
            testCase.verifyEqual(fdSqrt.pdf(xvals), fdIdFromSqrt101.pdf(xvals), 'AbsTol', 1E-14);
            testCase.verifyEqual(fdLog.pdf(xvals), fdIdFromLog101.pdf(xvals), 'AbsTol', 1E-14);
            
        end
        function testNegativitiyId(testCase)
            % For few coefficients and distributions will little common
            % probability mass, significantly negative values occur for
            % identity transformation.
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            kappa = 10;
            dist1 = VMDistribution(1, kappa);
            dist2 = VMDistribution(3, kappa);
            fd1Id = FourierDistribution.fromDistribution(dist1, 5, 'identity');
            fd2Id = FourierDistribution.fromDistribution(dist2, 5, 'identity');
            testCase.verifyWarning(@()fd1Id.multiply(fd2Id), 'Normalization:negative');
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:negative'));
            distMultId = fd1Id.multiply(fd2Id);
            xvals = 0:0.01:2 * pi;
            testCase.verifyTrue(any(distMultId.pdf(xvals) < -0.01));
        end
        function testNonnegativitySqrtAndLog(testCase)
            % No negative values occur for sqrt and log transformations in
            % this scenario.
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            testCase.applyFixture(SuppressedWarningsFixture('Multiply:NotNormalizing'));
            kappa = 10;
            dist1 = VMDistribution(1, kappa);
            dist2 = VMDistribution(3, kappa);
            for transformation = {'sqrt', 'log'}
                fd1 = FourierDistribution.fromDistribution(dist1, 5, [transformation{:}]);
                fd2 = FourierDistribution.fromDistribution(dist2, 5, [transformation{:}]);
                distMult = fd1.multiply(fd2);
                xvals = 0:0.01:2 * pi;
                testCase.verifyGreaterThanOrEqual(distMult.pdf(xvals), 0);
            end
        end
        function testShift(testCase)
            vm = VMDistribution(1, 10);
            vmShift = VMDistribution(2, 10);
            fdId = FourierDistribution.fromDistribution(vm, 41, 'identity');
            fdSqrt = FourierDistribution.fromDistribution(vm, 41, 'sqrt');
            fdIdShift = fdId.shift(1);
            fdSqrtShift = fdSqrt.shift(1);
            
            testCase.verifyClass(fdIdShift,'FourierDistribution');
            testCase.verifyClass(fdSqrtShift,'FourierDistribution');
            
            xvals = 0:0.01:2 * pi;
            testCase.verifyEqual(fdIdShift.pdf(xvals), vmShift.pdf(xvals), 'AbsTol', 1E-6);
            testCase.verifyEqual(fdSqrtShift.pdf(xvals), vmShift.pdf(xvals), 'AbsTol', 1E-6);
        end
    end
end

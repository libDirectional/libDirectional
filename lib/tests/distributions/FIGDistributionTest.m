classdef FIGDistributionTest < matlab.unittest.TestCase
    methods(Static)
        function testGridConversion(testCase, dist, coeffs, enforceNonnegative, tolerance)
            figd = FIGDistribution.fromDistribution(dist, coeffs, enforceNonnegative);
            % Test grid values
            xvals = linspace(0, 2*pi, coeffs+1);
            xvals(end) = [];
            testCase.verifyEqual(figd.pdf(xvals), dist.pdf(xvals), 'AbsTol', tolerance);
            % Test approximation of pdf
            xvals = -2 * pi:0.01:3 * pi;
            testCase.verifyEqual(figd.pdf(xvals), dist.pdf(xvals), 'AbsTol', tolerance);
        end
    end
    
    methods(Test)
        % Test conversions
        function testVMToGridId(testCase)
            mu = 0.4;
            for kappa = .1:.1:2
                dist = VMDistribution(mu, kappa);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, false, 1E-8);
            end
        end
        function testVMToGridSqrt(testCase)
            mu = 0.5;
            for kappa = .1:.1:2
                dist = VMDistribution(mu, kappa);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, true, 1E-8);
            end
        end
        function testWNToGridId(testCase)
            mu = 0.8;
            for sigma = .2:.1:2
                dist = WNDistribution(mu, sigma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, false, 1E-8);
            end
        end
        function testWNToGridSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            mu = 0.9;
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'))
            for sigma = .2:.1:2
                dist = WNDistribution(mu, sigma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, true, 1E-8);
            end
        end
        function testWCToGridId(testCase)
            mu = 1.2;
            for gamma = .8:.1:3
                dist = WCDistribution(mu, gamma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, false, 1E-7);
            end
        end
        
        function testWCToGridSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            mu = 1.3;
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:ApproximationHypergeometric'))
            for gamma = .8:.1:3
                dist = WCDistribution(mu, gamma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, true, 1E-7);
            end
        end
        
        function testWEToGridId(testCase)
            % Treat differently due to lack of continuity
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            for lambda=.1:.1:2
                xvals=-2*pi:0.01:3*pi;
                xvals=xvals(mod(xvals,2*pi)>0.5 & mod(xvals,2*pi)<(2*pi-0.5));
                dist=WEDistribution(lambda);
                gd=FIGDistribution.fromDistribution(dist,3001, false);
                testCase.verifySize(gd.gridValues,[3001,1]);
                testCase.verifyEqual(gd.pdf(xvals),dist.pdf(xvals),'AbsTol', 5E-3);
            end
        end
        function testWEToGridSqrt(testCase)
            % For sqrt, same applies as to identity
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            for lambda=.1:.1:2
                xvals=-2*pi:0.01:3*pi;
                xvals=xvals(mod(xvals,2*pi)>0.5 & mod(xvals,2*pi)<(2*pi-0.5));
                dist=WEDistribution(lambda);
                gd=FIGDistribution.fromDistribution(dist,3001, true);
                testCase.verifySize(gd.gridValues,[3001,1]);
                testCase.verifyEqual(gd.pdf(xvals),dist.pdf(xvals),'AbsTol', 5E-3);
            end
        end
        function testWLToGridId(testCase)
            % Only test parameter combinations that don't result in too
            % abrupt changes in pdf values
            for lambda=0.1:0.2:1
                for kappa=0.1:0.5:4
                    dist=WLDistribution(lambda,kappa);
                    FIGDistributionTest.testGridConversion(testCase,dist,1001,false,1E-3);
                end
            end
        end
        function testWLToGridSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            for lambda=0.1:0.2:1
                for kappa=0.1:0.5:4
                    dist=WLDistribution(lambda,kappa);
                    FIGDistributionTest.testGridConversion(testCase,dist,1001,true,1E-3);
                end
            end
        end
        function testFourierIdToGridNoEnforcementAllValsNonnegative(testCase)
            % Case 1: Fourier identity to FIGFAN without negative (accurate)
            % Transsform from Fourier identity to FIG without enforcing
            % nonnegativity for a density that does not have negative
            % values. This should be an accurate transformation.
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(1,1);
            fd = FourierDistribution.fromDistribution(vm,3,'identity');
            % Verify grid values
            gd = FIGDistribution.fromDistribution(fd,3,false);
            testCase.verifyEqual(gd.gridValues, (ifft(ifftshift(fd.c),'symmetric')*numel(fd.c))');
            % Verify density
            testCase.verifyWarningFree(@()...
                FIGDistributionTest.testGridConversion(testCase,fd,3,false,2E-16));
            % Also the same when oversampling
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            testCase.verifyWarningFree(@()...
                FIGDistributionTest.testGridConversion(testCase,fd,11,false,1E-16));
        end
        function testFourierIdToGridNoEnforcementSomeValsNegative(testCase)
            % Case 2: Fourier Identity to FIGFAN with negative (cannot be accurate because FIG only allows nonnegative grid values)
            % Transsform from Fourier identity to FIG without enforcing
            % nonnegativity for a density that has negative
            % values. This should *not* be an accurate transformation.
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            
            vm = VMDistribution(1,10);
            fd = FourierDistribution.fromDistribution(vm,3,'identity');
            testCase.verifyWarning(@()FIGDistribution.fromDistribution(fd,3,false),'FourierToFIG:ImpreciseId');
            testCase.verifyWarning(@()FIGDistribution.fromDistribution(fd,3,false),'Normalization:notNormalized');
            
            testCase.applyFixture(SuppressedWarningsFixture('FourierToFIG:ImpreciseId'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            gd = FIGDistribution.fromDistribution(fd,3,false);
            % Verify that result deviates by at least 0.1
            testCase.verifyThat(gd.gridValues,~IsEqualTo(ifft(ifftshift(fd.c),'symmetric')*numel(fd.c),'Within',AbsoluteTolerance(0.1)))
        end
        function testFourierSqrtToGridWithEnforcementAllValsNonnegative(testCase)
            % Case 3: Fourier Square without negative values to FIGFAN (accurate)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(1,1);
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fd = FourierDistribution.fromDistribution(vm,3,'sqrt');
            % Verify grid values
            gd = FIGDistribution.fromDistribution(fd,3,true);
            testCase.verifyEqual(sqrt(gd.gridValues), (ifft(ifftshift(fd.c),'symmetric')*numel(fd.c))');
            % Verify density
            testCase.verifyWarningFree(@()...
                FIGDistributionTest.testGridConversion(testCase,fd,3,true,3E-16));
            % Also the same when oversampling
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            testCase.verifyWarningFree(@()...
                FIGDistributionTest.testGridConversion(testCase,fd,11,true,3E-16));
        end
        function testFourierSqrtToGridWithEnforcementSomeValsNegativeNoPad(testCase)
            %% Case 4: Fourier Square with negative values to FIGFDN without padding (inaccurate)
            import matlab.unittest.constraints.AbsoluteTolerance
            import matlab.unittest.constraints.IsEqualTo
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            
            vm = VMDistribution(1,10);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fd = FourierDistribution.fromDistribution(vm,3,'sqrt');
            fixture.teardown();
            testCase.verifyWarning(@()FIGDistribution.fromDistribution(fd,3,false),'FourierToFIG:OtherInterpolationWithInsufficientCoeffs');

            testCase.applyFixture(SuppressedWarningsFixture('FourierToFIG:OtherInterpolationWithInsufficientCoeffs'));
            gd = FIGDistribution.fromDistribution(fd,3,false);
            % Verify that result deviates by at least 0.1
            testCase.verifyThat(sqrt(gd.gridValues),~IsEqualTo(ifft(ifftshift(fd.c),'symmetric')*numel(fd.c),'Within',AbsoluteTolerance(0.1)))
        end
        function testFourierSqrtToGridWithEnforcementAllValsNonnegativeWithPad(testCase)
            % Case 5: Fourier Square with negative values to FIGFAN with padding (accurate)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(1,10);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            fd = FourierDistribution.fromDistribution(vm,3,'sqrt');
            fixture.teardown();
            % Verify density
            testCase.verifyWarningFree(@()...
                FIGDistributionTest.testGridConversion(testCase,fd,5,false,4E-16));
            % Also the same when oversampling
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            testCase.verifyWarningFree(@()...
                FIGDistributionTest.testGridConversion(testCase,fd,11,false,4E-16));
        end
        
        function testCircularUniformToGridId(testCase)
            dist = CircularUniformDistribution();
            FIGDistributionTest.testGridConversion(testCase, dist, 101, false, 1E-8);
        end
        function testCircularUniformToGridSqrt(testCase)
            dist = CircularUniformDistribution();
            FIGDistributionTest.testGridConversion(testCase, dist, 101, true, 1E-8);
        end
        
        function testGCMToGridId(testCase)
            vm = VMDistribution(1, 2);
            wn = WNDistribution(2, 1);
            dist = CircularMixture({vm, wn}, [0.3, 0.7]);
            FIGDistributionTest.testGridConversion(testCase, dist, 101, false, 1E-8);
        end
        function testGCMToGridSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            vm = VMDistribution(1, 2);
            wn = WNDistribution(2, 1);
            dist = CircularMixture({vm, wn}, [0.3, 0.7]);
            FIGDistributionTest.testGridConversion(testCase, dist, 101, true, 1E-8);
        end
        function testCCDToGridId(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormula'));
            vm = VMDistribution(1, 2);
            dist = CustomCircularDistribution(@(x)vm.pdf(x));
            FIGDistributionTest.testGridConversion(testCase, dist, 101, false, 1E-8);
        end
        function testCCDToGridSqrt(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormula'));
            vm = VMDistribution(1, 2);
            dist = CustomCircularDistribution(@(x)vm.pdf(x));
            FIGDistributionTest.testGridConversion(testCase, dist, 101, true, 1E-8);
        end
        function testValue(testCase)
            vm = VMDistribution(3, 1);
            fd = FIGDistribution.fromDistribution(vm, 99, true);
            testCase.verifyEqual(fd.value(0), vm.pdf(0), 'AbsTol', 1e-16);
            testCase.verifyEqual(fd.value(eps), 0, 'AbsTol', 1e-16);
        end
        function testFromFunction(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for kappa = 0.1:0.3:4
                vm = VMDistribution(3, kappa);
                fd = FIGDistribution.fromFunction(@(x)vm.pdf(x), 501, true);
                testCase.verifyEqual(fd.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        function testFromFunctionValues(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            vm = VMDistribution(3, 1);
            fvals = vm.pdf(linspace(0, 2*pi, 100))';
            fvals(end) = [];
            fd1 = FIGDistribution.fromFunctionValues(fvals, 99, true);
            % Test truncation in fromFunction
            fvals = vm.pdf(linspace(0, 2*pi, 101))';
            fvals(end) = [];
            fd2 = FIGDistribution.fromFunctionValues(fvals, 50, true);
            testCase.verifyEqual(fd1.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            testCase.verifyEqual(fd2.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
        end
        
        % Test operations for prediction and filter steps
        function testMultiplyVM(testCase)
            for currEnforcement = [false, true]
                xvals = -2 * pi:0.01:3 * pi;
                for kappa1 = 0.1:0.3:4
                    for kappa2 = 0.1:0.3:4
                        dist1 = VMDistribution(0, kappa1);
                        dist2 = VMDistribution(0, kappa2);
                        f1 = FIGDistribution.fromDistribution(dist1, 101, currEnforcement);
                        f2 = FIGDistribution.fromDistribution(dist2, 101, currEnforcement);
                        fFiltered = f1.multiply(f2);
                        distResult = dist1.multiply(dist2);
                        testCase.verifyWarningFree(@()fFiltered.pdf(xvals));
                        testCase.verifyEqual(fFiltered.pdf(xvals), distResult.pdf(xvals), 'AbsTol', 1E-8);
                    end
                end
            end
        end
        function testConvolveWN(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for currEnforcement = [false, true]
                for sigma1 = 0.1:0.3:4
                    for sigma2 = 0.1:0.3:4
                        dist1 = WNDistribution(0, sigma1);
                        dist2 = WNDistribution(0, sigma2);
                        f1 = FIGDistribution.fromDistribution(dist1, 101, currEnforcement);
                        f2 = FIGDistribution.fromDistribution(dist2, 101, currEnforcement);
                        fPredicted = f1.convolve(f2);
                        distResult = dist1.convolve(dist2);
                        testCase.verifyEqual(fPredicted.pdf(xvals), distResult.pdf(xvals), 'AbsTol', 1E-8);
                    end
                end
            end
        end
        
        function testMomentsGridSqrtPredictionVM(testCase)
            for kappa1 = [1,2]
                for kappa2 = [1,2]
                    vm1 = VMDistribution(0, kappa1);
                    vm2 = VMDistribution(0, kappa2);
                    vmRes = vm1.convolve(vm2);
                    f1 = FIGDistribution.fromDistribution(vm1, 101, true);
                    f2 = FIGDistribution.fromDistribution(vm2, 101, true);
                    fPredicted = f1.convolve(f2);
                    testCase.verifyEqual(fPredicted.trigonometricMoment(1), vmRes.trigonometricMoment(1), 'AbsTol', 1E-8);
                end
            end
        end
        
        function testNormalizationPdfSinc(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(1,50);
            % It will have to renormalize, only check interpolation below
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            gd = FIGDistribution.fromDistribution(vm,3,true);
            fixture.teardown();
            
            int = integral(@(x)reshape(gd.pdf(x(:)',true,71),size(x)),0,2*pi);
            testCase.verifyEqual(int,1,'AbsTol',1e-4);
            
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            gd = FIGDistribution.fromDistribution(vm,7,false);
            fixture.teardown();
            
            int = integral(@(x)reshape(gd.pdf(x(:)',true,71),size(x)),0,2*pi);
            testCase.verifyEqual(int,1,'AbsTol',1e-5);
        end
        % CDF tests
        function testCdfIdentityStartZero(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            fd = FIGDistribution.fromDistribution(VMDistribution(1, 3), 101, false);
            xvals = 0:0.01:2 * pi;
            testCase.applyFixture(SuppressedWarningsFixture('Grid:CanOnlyIntegrateWithLimitsNumerically'));
            intValsNumerically = arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend), xvals);
            testCase.verifyEqual(fd.cdf(xvals, 0), intValsNumerically, 'AbsTol', 1E-8);
        end
        function testCdfSqrtStartZero(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            fd = FIGDistribution.fromDistribution(VMDistribution(1, 3), 101, true);
            xvals = 0:0.01:2 * pi;
            testCase.applyFixture(SuppressedWarningsFixture('Grid:CanOnlyIntegrateWithLimitsNumerically'));
            intValsNumerically = arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend), xvals);
            testCase.verifyEqual(fd.cdf(xvals, 0), intValsNumerically, 'AbsTol', 1E-8);
        end
        function testCdfIdentityStartNonzero(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            fd = FIGDistribution.fromDistribution(VMDistribution(1, 3), 101, false);
            xvals = 0:0.01:2 * pi;
            startingPoint = 1;
            testCase.applyFixture(SuppressedWarningsFixture('Grid:CanOnlyIntegrateWithLimitsNumerically'));
            intValsNumerically = [arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend)+integral(@(x)fd.pdf(x), startingPoint, 2*pi), xvals(xvals < startingPoint)), ...
                arrayfun(@(xend)integral(@(x)fd.pdf(x), startingPoint, xend), xvals(xvals >= startingPoint))];
            testCase.verifyEqual(fd.cdf(xvals, startingPoint), intValsNumerically, 'AbsTol', 1E-8);
        end
        % Other tests
        function testTruncation(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vm = VMDistribution(3, 1);
            fd1 = FIGDistribution.fromDistribution(vm, 100, true);
            testCase.verifyWarning(@()fd1.truncate(102), 'Truncate:TooFewGridPoints');
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewGridPoints'));
            fd2 = fd1.truncate(102);
            testCase.verifySize(fd2.gridValues,[102,1]);
            fd3 = fd2.truncate(51);
            testCase.verifySize(fd3.gridValues,[51,1]);
        end
        function testMoments(testCase)
            vm = VMDistribution(2, 1);
            fd1 = FIGDistribution.fromDistribution(vm, 15, false);
            fd2 = FIGDistribution.fromDistribution(vm, 15, true);
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
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:DownsampleViaFourier'));
            dist = VMDistribution(0, 5);
            fdNonnegative = FIGDistribution.fromDistribution(dist, 1000, true);
            testCase.verifyEqual(integral(@(x)fdNonnegative.pdf(x), 0, 2*pi), 1, 'RelTol', 1E-4);
            fdSqrtTrunc = fdNonnegative.truncate(7);
            % Verify that use of .truncate results in a normalized density
            testCase.verifyEqual(integral(@(x)fdSqrtTrunc.pdf(x), 0, 2*pi), 1, 'RelTol', 1E-4);
            % Verify that no warning is given due to normalization in
            % Fourier
            testCase.verifyWarningFree(@()integral(@(x)fdSqrtTrunc.pdf(x), 0, 2*pi));
        end
        function testIntegral(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist = VMDistribution(0, 5);
            testCase.applyFixture(SuppressedWarningsFixture('Grid:CanOnlyIntegralNumericallyWithLimits'));
            for transformation = {false, true}
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
                fd = FIGDistribution.fromDistribution(dist, 15, [transformation{:}]);
                fixture.teardown();
                testCase.verifyEqual(fd.integral(0, 1.5), fd.integralNumerical(0, 1.5), 'RelTol', 1E-8);
                testCase.verifyEqual(fd.integral(1.5, 0), fd.integralNumerical(1.5, 0), 'RelTol', 1E-8);
                testCase.verifyEqual(fd.integral(10, -10), fd.integralNumerical(10, -10), 'RelTol', 1E-8);
            end
        end
        function testPdfOnGrid(testCase)
            for transformation = {false, true}
                dist = VMDistribution(0, 5);
                fd = FIGDistribution.fromDistribution(dist, 15, [transformation{:}]);
                [vals, xgrid] = fd.pdfOnGrid(15);
                testCase.verifyEqual(vals, fd.pdf(xgrid), 'RelTol', 1E-8);
            end
        end
        
        function testShift(testCase)
            vm = VMDistribution(1, 10);
            vmShift = VMDistribution(2, 10);
            gdId = FIGDistribution.fromDistribution(vm, 51, false);
            gdSqrt = FIGDistribution.fromDistribution(vm, 51, true);
            gdIdShift = gdId.shift(1);
            gdSqrtShift = gdSqrt.shift(1);
            
            testCase.verifyClass(gdIdShift, 'FIGDistribution');
            testCase.verifyClass(gdSqrtShift, 'FIGDistribution');
            
            xvals = 0:0.01:2 * pi;
            testCase.verifyEqual(gdIdShift.pdf(xvals), vmShift.pdf(xvals), 'AbsTol', 1E-6);
            testCase.verifyEqual(gdSqrtShift.pdf(xvals), vmShift.pdf(xvals), 'AbsTol', 1E-6);
        end
        function testPdfViaSinc(testCase)
            dist = WNDistribution(pi, 1);
            noPoints = 5;
            currTrans = false;
            gd = FIGDistribution.fromDistribution(dist, noPoints, currTrans);
            gd.enforcePdfNonnegative = false;

            stepSize=2*pi/noPoints;
            evalAt=linspace(0,2*pi,100);
            % Only considering sinc function at points from 0 to 2*pi (this
            % yields lowest quality)
            
            fvals = [gd.gridValues;gd.gridValues(1)];
            % Instead of 0:noPoints could use linspace(0,2*pi,noPoints+1)/stepSize
            sincFun = @(x)arrayfun(@(xCurr)sinc(xCurr/stepSize-(0:noPoints))*fvals, x);
            testCase.verifyEqual(sincFun(evalAt), gd.pdf(evalAt), 'AbsTol', 0.003);
            % Considering sinc function at points from -2*pi to 4*pi
            % (higher quality)
            xvals = linspace(-2*pi, 4*pi, noPoints*3 + 1);
            fvals = [repmat(gd.gridValues,3,1);gd.gridValues(1)];
            sincFun = @(x)arrayfun(@(xCurr)sinc(xCurr/stepSize-xvals/stepSize)*fvals, x);
            testCase.verifyEqual(sincFun(evalAt), gd.pdf(evalAt), 'AbsTol', 0.0005);
            % Considering sinc function at points from -4*pi to 6*pi
            % (even higher quality)
            xvals = linspace(-4*pi, 6*pi, noPoints*5 + 1);
            fvals = [repmat(gd.gridValues,5,1);gd.gridValues(1)];
            sincFun = @(x)arrayfun(@(xCurr)sinc(xCurr/stepSize-xvals/stepSize)*fvals, x);
            testCase.verifyEqual(sincFun(evalAt), gd.pdf(evalAt), 'AbsTol', 0.0002);
        end
        function testGetGridPoint(testCase)
            dist = VMDistribution(1, 1);
            gd = FIGDistribution.fromDistribution(dist, 37);
            entireGrid = gd.getGrid();
            testCase.verifyEqual(gd.getGridPoint([10,11]),entireGrid(:,[10,11]),'AbsTol',1e-15);
        end
        function testGetClosestPoint(testCase)
            gd = FIGDistribution.fromDistribution(VMDistribution(0,1),11);
            closestPoints = gd.getClosestPoint([3*pi/gd.noOfGridPoints-eps,3*pi/gd.noOfGridPoints+eps]);
            testCase.verifyEqual(closestPoints,2*pi/gd.noOfGridPoints*[1,2]);
        end
    end
end

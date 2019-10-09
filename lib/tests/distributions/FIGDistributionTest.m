classdef FIGDistributionTest < matlab.unittest.TestCase
    properties
    end
    
    methods(Static)
        function testGridConversion(testCase, dist, coeffs, transformation, tolerance)
            gd = FIGDistribution.fromDistribution(dist, coeffs, transformation);
            % Test grid values
            xvals = linspace(0, 2*pi, coeffs+1);
            xvals(end) = [];
            switch transformation
                case 'identity'
                    testCase.verifyEqual(gd.gridValues', dist.pdf(xvals), 'AbsTol', tolerance);
                case 'sqrt'
                    testCase.verifyEqual(gd.gridValues'.^2, dist.pdf(xvals), 'AbsTol', tolerance);
                otherwise
                    error('Transformation unsupported');
            end
            % Test approximation of pdf
            xvals = -2 * pi:0.01:3 * pi;
            testCase.verifyEqual(gd.pdf(xvals), dist.pdf(xvals), 'AbsTol', tolerance);
        end
    end
    
    methods(Test)
        % Test conversions
        function testVMToGridId(testCase)
            mu = 0.4;
            for kappa = .1:.1:2
                dist = VMDistribution(mu, kappa);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, 'identity', 1E-8);
            end
        end
        function testVMToGridSqrt(testCase)
            mu = 0.5;
            for kappa = .1:.1:2
                dist = VMDistribution(mu, kappa);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, 'sqrt', 1E-8);
            end
        end
        function testWNToGridId(testCase)
            mu = 0.8;
            for sigma = .2:.1:2
                dist = WNDistribution(mu, sigma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, 'identity', 1E-8);
            end
        end
        function testWNToGridSqrt(testCase)
            mu = 0.9;
            warningSettings = warning('off', 'Conversion:NoFormulaSqrt');
            for sigma = .2:.1:2
                dist = WNDistribution(mu, sigma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, 'sqrt', 1E-8);
            end
            warning(warningSettings);
        end
        function testWCToGridId(testCase)
            mu = 1.2;
            for gamma = .8:.1:3
                dist = WCDistribution(mu, gamma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, 'identity', 1E-7);
            end
        end
        
        function testWCToGridSqrt(testCase)
            mu = 1.3;
            warningSettings = warning('off', 'Conversion:ApproximationHypergeometric');
            for gamma = .8:.1:3
                dist = WCDistribution(mu, gamma);
                FIGDistributionTest.testGridConversion(testCase, dist, 101, 'sqrt', 1E-7);
            end
            warning(warningSettings);
        end
        
        function testWEToGridId(testCase)
            %Treat differently due to lack of continuity
            warningSettings=warning('off','Normalization:notNormalized');
            for lambda=.1:.1:2
                xvals=-2*pi:0.01:3*pi;
                xvals=xvals(mod(xvals,2*pi)>0.5 & mod(xvals,2*pi)<(2*pi-0.5));
                dist=WEDistribution(lambda);
                gd=FIGDistribution.fromDistribution(dist,3001,'identity');
                testCase.verifySize(gd.gridValues,[3001,1]);
                testCase.verifyEqual(gd.pdf(xvals),dist.pdf(xvals),'AbsTol', 5E-3);
            end
            warning(warningSettings);
        end
        function testWEToGridSqrt(testCase)
            %For sqrt, same applies as to identity
            warningSettings=warning('off','Normalization:notNormalized');
            for lambda=.1:.1:2
                xvals=-2*pi:0.01:3*pi;
                xvals=xvals(mod(xvals,2*pi)>0.5 & mod(xvals,2*pi)<(2*pi-0.5));
                dist=WEDistribution(lambda);
                gd=FIGDistribution.fromDistribution(dist,3001,'sqrt');
                testCase.verifySize(gd.gridValues,[3001,1]);
                testCase.verifyEqual(gd.pdf(xvals),dist.pdf(xvals),'AbsTol', 5E-3);
            end
            warning(warningSettings);
        end
        function testWLToGridId(testCase)
            %Only test parameter combinations that don't result in too
            %abrupt changes in pdf values
            for lambda=0.1:0.2:1
                for kappa=0.1:0.5:4
                    dist=WLDistribution(lambda,kappa);
                    FIGDistributionTest.testGridConversion(testCase,dist,1001,'identity',1E-3);
                end
            end
        end
        function testWLToFourierSqrt(testCase)
            warningSettings=warning('off','Conversion:NoFormulaSqrt');
            for lambda=0.1:0.2:1
                for kappa=0.1:0.5:4
                    dist=WLDistribution(lambda,kappa);
                    FIGDistributionTest.testGridConversion(testCase,dist,1001,'sqrt',1E-3);
                end
            end
            warning(warningSettings);
        end
        function testCircularUniformToGridId(testCase)
            dist = CircularUniformDistribution();
            FIGDistributionTest.testGridConversion(testCase, dist, 101, 'identity', 1E-8);
        end
        function testCircularUniformToGridSqrt(testCase)
            dist = CircularUniformDistribution();
            FIGDistributionTest.testGridConversion(testCase, dist, 101, 'sqrt', 1E-8);
        end
        
        function testGCMToGridId(testCase)
            vm = VMDistribution(1, 2);
            wn = WNDistribution(2, 1);
            dist = CircularMixture({vm, wn}, [0.3, 0.7]);
            FIGDistributionTest.testGridConversion(testCase, dist, 101, 'identity', 1E-8);
        end
        function testGCMToGridSqrt(testCase)
            warningSettings = warning('off', 'Conversion:NoFormulaSqrt');
            vm = VMDistribution(1, 2);
            wn = WNDistribution(2, 1);
            dist = CircularMixture({vm, wn}, [0.3, 0.7]);
            FIGDistributionTest.testGridConversion(testCase, dist, 101, 'sqrt', 1E-8);
            warning(warningSettings);
        end
        function testCCDToGridId(testCase)
            warningSettings = warning('off', 'Conversion:NoFormula');
            vm = VMDistribution(1, 2);
            dist = CustomCircularDistribution(@(x)vm.pdf(x));
            FIGDistributionTest.testGridConversion(testCase, dist, 101, 'identity', 1E-8);
            warning(warningSettings);
        end
        function testCCDToGridSqrt(testCase)
            warningSettings = warning('off', 'Conversion:NoFormula');
            vm = VMDistribution(1, 2);
            dist = CustomCircularDistribution(@(x)vm.pdf(x));
            FIGDistributionTest.testGridConversion(testCase, dist, 101, 'sqrt', 1E-8);
            warning(warningSettings);
        end
        
        function testFromFunction(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for kappa = 0.1:0.3:4
                vm = VMDistribution(3, kappa);
                fd = FIGDistribution.fromFunction(@(x)vm.pdf(x), 501, 'sqrt');
                testCase.verifyEqual(fd.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        function testFromFunctionValues(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            vm = VMDistribution(3, 1);
            fvals = vm.pdf(linspace(0, 2*pi, 100))';
            fvals(end) = [];
            fd1 = FIGDistribution.fromFunctionValues(fvals, 99, 'sqrt');
            % Test truncation in fromFunction
            fvals = vm.pdf(linspace(0, 2*pi, 101))';
            fvals(end) = [];
            fd2 = FIGDistribution.fromFunctionValues(fvals, 50, 'sqrt');
            testCase.verifyEqual(fd1.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
            testCase.verifyEqual(fd2.pdf(xvals), vm.pdf(xvals), 'AbsTol', 1E-8);
        end
        
        % Test oeprations for prediction and filter steps
        function testMultiplyVM(testCase)
            for transformation={'identity', 'sqrt'}
                currTrans=[transformation{:}];
                xvals = -2 * pi:0.01:3 * pi;
                for kappa1 = 0.1:0.3:4
                    for kappa2 = 0.1:0.3:4
                        dist1 = VMDistribution(0, kappa1);
                        dist2 = VMDistribution(0, kappa2);
                        f1 = FIGDistribution.fromDistribution(dist1, 101, currTrans);
                        f2 = FIGDistribution.fromDistribution(dist2, 101, currTrans);
                        fFiltered = f1.multiply(f2);
                        distResult = dist1.multiply(dist2);
                        testCase.verifyEqual(fFiltered.pdf(xvals), distResult.pdf(xvals), 'AbsTol', 1E-8);
                    end
                end
            end
        end
        function testConvolveWN(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for transformation={'identity', 'sqrt'}
                currTrans=[transformation{:}];
                for sigma1 = 0.1:0.3:4
                    for sigma2 = 0.1:0.3:4
                        dist1 = WNDistribution(0, sigma1);
                        dist2 = WNDistribution(0, sigma2);
                        f1 = FIGDistribution.fromDistribution(dist1, 101, currTrans);
                        f2 = FIGDistribution.fromDistribution(dist2, 101, currTrans);
                        fPredicted = f1.convolve(f2);
                        distResult = dist1.convolve(dist2);
                        testCase.verifyEqual(fPredicted.pdf(xvals), distResult.pdf(xvals), 'AbsTol', 1E-8);
                    end
                end
            end
        end
        
        function testMomentsGridSqrtPredictionVM(testCase)
            for kappa1 = [1,2] %0.1:0.3:4
                for kappa2 = [1,2] %0.1:0.3:4
                    vm1 = VMDistribution(0, kappa1);
                    vm2 = VMDistribution(0, kappa2);
                    vmRes = vm1.convolve(vm2);
                    f1 = FIGDistribution.fromDistribution(vm1, 101, 'sqrt');
                    f2 = FIGDistribution.fromDistribution(vm2, 101, 'sqrt');
                    fPredicted = f1.convolve(f2);
                    testCase.verifyEqual(fPredicted.trigonometricMoment(1), vmRes.trigonometricMoment(1), 'AbsTol', 1E-8);
                end
            end
        end
        
        function testErrorsPredictAndFilter(testCase)
            warningSettings = warning('off', 'Normalization:cannotTest');
            fd1 = FIGDistribution.fromDistribution(VMDistribution(0, 1), 101, 'sqrt');
            fd2 = FIGDistribution.fromDistribution(VMDistribution(0, 1), 101, 'identity');
            testCase.verifyError(@()fd1.multiply(fd2), 'Multiply:differentTransformations');
            testCase.verifyError(@()fd1.convolve(fd2), 'Convolve:differentTransformations');
            fd3 = fd2.transformViaValues('square');
            testCase.verifyError(@()fd3.multiply(fd3), 'Multiply:unsupportedTransformation');
            testCase.verifyError(@()fd3.convolve(fd3), 'Convolve:unsupportedTransformation');
            warning(warningSettings);
        end
        
        % Test transformations
        function testSquaring(testCase)
            xvals = -2 * pi:0.01:3 * pi;
            for kappa = 0.1:0.3:4
                dist = VMDistribution(0, kappa);
                fd1 = FIGDistribution.fromDistribution(dist, 101, 'sqrt');
                fd2 = fd1.transformViaValues('square');
                testCase.verifyEqual(fd2.pdf(xvals), dist.pdf(xvals), 'AbsTol', 1E-8);
            end
        end
        
        function testSquareAfterSqrt(testCase)
            vm = VMDistribution(3, 1);
            fd1 = FIGDistribution.fromFunction(@(x)vm.pdf(x), 101, 'identity');
            fd2 = fd1.transformViaValues('sqrt');
            fd3 = fd2.transformViaValues('square');
            testCase.verifyEqual(fd3.transformation, 'identity');
            testCase.verifyEqual(fd3.gridValues, fd1.gridValues, 'AbsTol', 1E-8);
        end
        %{
        function testTransformToMoreCoeffs(testCase)
            vm = VMDistribution(3, 1);
            fdSqrt = FIGDistribution.fromDistribution(vm, 1001, 'sqrt');
            fdId = FIGDistribution.fromDistribution(vm, 101, 'identity');
            warningSettings = warning('off', 'Truncate:TooFewCoefficients');
            fdIdSqrt = fdId.transformViaFFT('sqrt', 1001);
            warning(warningSettings);
            testCase.verifyEqual(length(fdIdSqrt.a)+length(fdIdSqrt.b), 1001);
            testCase.verifyEqual(fdIdSqrt.a, fdSqrt.a, 'AbsTol', 1E-8);
            testCase.verifyEqual(fdIdSqrt.b, fdSqrt.b, 'AbsTol', 1E-8);
        end
        function testTransformToLessCoeffs(testCase)
            %case is easy enough to work with less coeffs
            vm = VMDistribution(3, 1);
            fdSqrt = FIGDistribution.fromDistribution(vm, 101, 'sqrt');
            fdId = FIGDistribution.fromDistribution(vm, 1001, 'identity');
            fdIdSqrt = fdId.transformViaFFT('sqrt', 101);
            testCase.verifyEqual(length(fdIdSqrt.a)+length(fdIdSqrt.b), 101, 'AbsTol', 1E-8);
            testCase.verifyEqual(fdIdSqrt.a, fdSqrt.a, 'AbsTol', 1E-8);
            testCase.verifyEqual(fdIdSqrt.b, fdSqrt.b, 'AbsTol', 1E-8);
        end
        %}
        % CDF tests
        function testCdfIdentityStartZero(testCase)
            fd = FIGDistribution.fromDistribution(VMDistribution(1, 3), 101, 'identity');
            xvals = 0:0.01:2 * pi;
            intValsNumerically = arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend), xvals);
            testCase.verifyEqual(fd.cdf(xvals, 0), intValsNumerically, 'AbsTol', 1E-8);
        end
        function testCdfSqrtStartZero(testCase)
            fd = FIGDistribution.fromDistribution(VMDistribution(1, 3), 101, 'sqrt');
            xvals = 0:0.01:2 * pi;
            intValsNumerically = arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend), xvals);
            testCase.verifyEqual(fd.cdf(xvals, 0), intValsNumerically, 'AbsTol', 1E-8);
        end
        function testCdfIdentityStartNonzero(testCase)
            fd = FIGDistribution.fromDistribution(VMDistribution(1, 3), 101, 'identity');
            xvals = 0:0.01:2 * pi;
            startingPoint = 1;
            intValsNumerically = [arrayfun(@(xend)integral(@(x)fd.pdf(x), 0, xend)+integral(@(x)fd.pdf(x), startingPoint, 2*pi), xvals(xvals < startingPoint)), ...
                arrayfun(@(xend)integral(@(x)fd.pdf(x), startingPoint, xend), xvals(xvals >= startingPoint))];
            testCase.verifyEqual(fd.cdf(xvals, startingPoint), intValsNumerically, 'AbsTol', 1E-8);
        end
        %other tests
        %{
        function testTruncation(testCase)
            vm = VMDistribution(3, 1);
            fd1 = FIGDistribution.fromDistribution(vm, 101, 'sqrt');
            testCase.verifyWarning(@()fd1.truncate(1001), 'Truncate:TooFewCoefficients');
            warningSettings = warning('off', 'Truncate:TooFewCoefficients');
            fd2 = fd1.truncate(1001);
            warning(warningSettings);
            testCase.verifyEqual(length(fd2.a)+length(fd2.b), 1001, 'AbsTol', 1E-8);
            fd3 = fd2.truncate(51);
            testCase.verifyEqual(length(fd3.a)+length(fd3.b), 51, 'AbsTol', 1E-8);
        end
        %}
        function testMoments(testCase)
            vm = VMDistribution(2, 1);
            fd1 = FIGDistribution.fromDistribution(vm, 15, 'identity');
            fd2 = FIGDistribution.fromDistribution(vm, 15, 'sqrt');
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
        %{
        function testNormalizationAfterTruncation(testCase)
            dist = VMDistribution(0, 5);
            fdSqrt = FIGDistribution.fromDistribution(dist, 1001, 'sqrt');
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
        %}
        function testIntegral(testCase)
            dist = VMDistribution(0, 5);
            for transformation = {'identity', 'sqrt'}
                warningSetting = warning('off', 'Normalization:cannotTest');
                fd = FIGDistribution.fromDistribution(dist, 15, [transformation{:}]);
                warning(warningSetting);
                testCase.verifyEqual(fd.integral(0, 1.5), fd.integralNumerical(0, 1.5), 'RelTol', 1E-8);
                testCase.verifyEqual(fd.integral(1.5, 0), fd.integralNumerical(1.5, 0), 'RelTol', 1E-8);
                testCase.verifyEqual(fd.integral(10, -10), fd.integralNumerical(10, -10), 'RelTol', 1E-8);
            end
        end
        function testPdfOnGrid(testCase)
            for transformation = {'identity', 'sqrt'}
                dist = VMDistribution(0, 5);
                fd = FIGDistribution.fromDistribution(dist, 15, [transformation{:}]);
                [vals, xgrid] = fd.pdfOnGrid(15);
                testCase.verifyEqual(vals, fd.pdf(xgrid), 'RelTol', 1E-8);
            end
        end
        
        function testShift(testCase)
            vm = VMDistribution(1, 10);
            vmShift = VMDistribution(2, 10);
            gdId = FIGDistribution.fromDistribution(vm, 51, 'identity');
            gdSqrt = FIGDistribution.fromDistribution(vm, 51, 'sqrt');
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
            currTrans = 'identity';
            gd = FIGDistribution.fromDistribution(dist, noPoints, currTrans);
            gd.enforcePdfNonnegative = false;

            stepSize=2*pi/noPoints;
            evalAt=linspace(0,2*pi,100);
            % Only considering sinc function at points from 0 to 2*pi (this
            % yields lowest quality)
            xvals = linspace(0,2*pi,noPoints+1);
            fvals = [gd.gridValues;gd.gridValues(1)];
            sincFun = @(x)arrayfun(@(xCurr)sinc(xCurr/stepSize-xvals/stepSize)*fvals, x);
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
    end
end

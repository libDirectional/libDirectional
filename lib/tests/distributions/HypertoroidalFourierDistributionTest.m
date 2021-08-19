classdef HypertoroidalFourierDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testConstructor(testCase)
            hfd = HypertoroidalFourierDistribution([0, 0, 0; 0, 1 / sqrt(2*pi)^2, 0; 0, 0, 0]);
            testCase.verifySize(hfd.C, [3, 3]);
            testCase.verifyEqual(hfd.transformation, 'sqrt');
            testCase.verifyWarning(@()HypertoroidalFourierDistribution(1/sqrt(2*pi)), 'fourierCoefficients:singleCoefficient');
            testCase.verifyError(@()HypertoroidalFourierDistribution(WNDistribution(0, 1)), 'fourierCoefficients:invalidCoefficientMatrix');
        end
        
        function testNormalization2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            unnormalizedCoeffs2D = fftshift(fftn(rand(3, 7)+0.5));
            unnormalizedCoeffs2D(2, 4) = 1;
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdId = HypertoroidalFourierDistribution(unnormalizedCoeffs2D, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution(unnormalizedCoeffs2D, 'sqrt');
            fixture.teardown;
            testCase.verifyEqual(integral2(@(x, y)reshape(hfdId.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi), 1, 'RelTol', 1E-4);
            testCase.verifyEqual(integral2(@(x, y)reshape(hfdSqrt.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi), 1, 'RelTol', 1E-4);
            % Test warnings
            testCase.verifyWarning(@()HypertoroidalFourierDistribution([0, 0, 0; 0, -1, 0; 0, 0, 0], 'identity'), 'Normalization:negative');
            testCase.verifyError(@()HypertoroidalFourierDistribution([0, 0, 0; 0, 1E-201, 0; 0, 0, 0], 'identity'), 'Normalization:almostZero');
        end
        
        function testNormalization3D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                unnormalizedCoeffs3D = fftshift(fftn(rand(3, 11, 7)+0.5));
                unnormalizedCoeffs3D(2, 6, 4) = 1;
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
                hfdId = HypertoroidalFourierDistribution(unnormalizedCoeffs3D, 'identity');
                hfdSqrt = HypertoroidalFourierDistribution(unnormalizedCoeffs3D, 'sqrt');
                fixture.teardown;
                testCase.verifyEqual(integral3(@(x, y, z)reshape(hfdId.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, 2*pi, 0, 2*pi, 0, 2*pi), 1, 'RelTol', 1E-4);
                testCase.verifyEqual(integral3(@(x, y, z)reshape(hfdSqrt.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, 2*pi, 0, 2*pi, 0, 2*pi), 1, 'RelTol', 1E-4);
            end
        end
        
        function testTruncation1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.applyFixture(SuppressedWarningsFixture('Conversion:NoFormulaSqrt'));
            hfd1 = HypertoroidalFourierDistribution.fromDistribution(WNDistribution(1, 1), 101);
            hfd2 = hfd1.truncate(51);
            testCase.verifySize(hfd2.C, [51, 1]);
            xvals = linspace(0, 2*pi, 100);
            testCase.verifyEqual(hfd1.pdf(xvals), hfd2.pdf(xvals), 'AbsTol', 1E-8);
        end
        
        function testTruncation3D(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            
            coeffs = [15, 15, 15];
            mu = [0; 0; 0];
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            hwnd = HypertoroidalWNDistribution(mu, C);
            hfdId = HypertoroidalFourierDistribution.fromDistribution(hwnd, coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromDistribution(hwnd, coeffs, 'sqrt');
            hfdIdTruncWithoutNormalization = hfdId;
            hfdIdTruncWithoutNormalization.C = hfdId.C(7:9, 7:9, 7:9);
            hfdSqrtTruncWithoutNormalization = hfdSqrt;
            hfdSqrtTruncWithoutNormalization.C = hfdSqrt.C(7:9, 7:9, 7:9);
            % Test warnings
            testCase.verifyWarning(@()hfdId.truncate([21, 21, 21]), 'Truncate:TooFewCoefficients');
            testCase.verifyWarning(@()hfdSqrt.truncate([21, 21, 21]), 'Truncate:TooFewCoefficients');
            
            if enableExpensive
                % No problems exist when using identity transformation
                testCase.verifyEqual(integral3(@(x, y, z)reshape(hfdIdTruncWithoutNormalization.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, 2*pi, 0, 2*pi, 0, 2*pi), 1, 'RelTol', 1E-4);
                % Show that in this example, naive truncation does not work for
                % sqrt case
                testCase.verifyTrue(abs(integral3(@(x, y, z)reshape(hfdSqrtTruncWithoutNormalization.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, 2*pi, 0, 2*pi, 0, 2*pi)-1) > 0.01);
                % Show that .truncate handles it correctly
                hfdSqrtTrunc = hfdSqrt.truncate([3, 3, 3]);
                testCase.verifyEqual(integral3(@(x, y, z)reshape(hfdSqrtTrunc.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, 2*pi, 0, 2*pi, 0, 2*pi), 1, 'RelTol', 1E-4);
            end
        end
        
        function testFromFunction2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                [xTest, yTest] = meshgrid(-pi:0.1:3*pi);
            else
                [xTest, yTest] = meshgrid(0:1:2*pi);
            end
            coeffs = [13, 15];
            mu = [1; 0];
            for kappa1 = [0.2, 1]
                for kappa2 = [0.5, 2]
                    for lambda = [0, 1]
                        tvm = ToroidalVMSineDistribution(mu, [kappa1; kappa2], lambda);
                        hfdId = HypertoroidalFourierDistribution.fromFunction( ...
                            @(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'identity');
                        hfdSqrt = HypertoroidalFourierDistribution.fromFunction( ...
                            @(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'sqrt');
                        fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
                        hfdLog = HypertoroidalFourierDistribution.fromFunction( ...
                            @(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'log');
                        fixture.teardown;
                        testCase.verifyClass(hfdId, 'HypertoroidalFourierDistribution');
                        testCase.verifyClass(hfdSqrt, 'HypertoroidalFourierDistribution');
                        testCase.verifyClass(hfdLog, 'HypertoroidalFourierDistribution');
                        testCase.verifySize(hfdId.C, coeffs);
                        testCase.verifySize(hfdSqrt.C, coeffs);
                        testCase.verifySize(hfdLog.C, coeffs);
                        testCase.verifyEqual(tvm.pdf([xTest(:)'; yTest(:)']), hfdId.pdf([xTest(:)'; yTest(:)']), 'AbsTol', 1E-5);
                        testCase.verifyEqual(tvm.pdf([xTest(:)'; yTest(:)']), hfdSqrt.pdf([xTest(:)'; yTest(:)']), 'AbsTol', 1E-6);
                        fixture = testCase.applyFixture(SuppressedWarningsFixture('pdf:mayNotBeNormalized'));
                        testCase.verifyEqual(tvm.pdf([xTest(:)'; yTest(:)']), hfdLog.pdf([xTest(:)'; yTest(:)']), 'AbsTol', 1E-6);
                        fixture.teardown;
                    end
                end
            end
            testCase.verifyError(@()HypertoroidalFourierDistribution.fromFunction( ...
                @(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'abc'), 'fromFunctionValues:unrecognizedTranformation')
        end
        
        function testFromFunction3D(testCase)
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
            coeffs = [25, 27, 23];
            hfdId = HypertoroidalFourierDistribution.fromFunction( ...
                @(x, y, z)reshape(hwnd.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromFunction( ...
                @(x, y, z)reshape(hwnd.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'sqrt');
            testCase.verifyClass(hfdId, 'HypertoroidalFourierDistribution')
            testCase.verifySize(hfdId.C, coeffs);
            testCase.verifyEqual(hfdId.pdf(testPoints), hwnd.pdf(testPoints), 'AbsTol', 1E-6);
            testCase.verifyEqual(hfdSqrt.pdf(testPoints), hwnd.pdf(testPoints), 'AbsTol', 1E-5);
        end
        
        function testFromFunction4D(testCase)
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                C = [0.7, 0.4, 0.2, -0.5; 0.4, 0.6, 0.1, 0; 0.2, 0.1, 1, -0.3; -0.5, 0, -0.3, 0.9] * 2;
                mu = 2 * pi * rand(4, 1);
                hwnd = HypertoroidalWNDistribution(mu, C);
                coeffs = [19, 17, 15, 21];
                hfdId = HypertoroidalFourierDistribution.fromFunction( ...
                    @(x, y, z, w)reshape(hwnd.pdf([x(:)'; y(:)'; z(:)'; w(:)']), size(x)), coeffs, 'identity');
                hfdSqrt = HypertoroidalFourierDistribution.fromFunction( ...
                    @(x, y, z, w)reshape(hwnd.pdf([x(:)'; y(:)'; z(:)'; w(:)']), size(x)), coeffs, 'sqrt');
                testCase.verifyClass(hfdId, 'HypertoroidalFourierDistribution')
                testCase.verifyClass(hfdSqrt, 'HypertoroidalFourierDistribution')
                testCase.verifySize(hfdId.C, coeffs);
                testCase.verifySize(hfdSqrt.C, coeffs);
                testPoints = rand(4, 100);
                testCase.verifyEqual(hfdId.pdf(testPoints), hwnd.pdf(testPoints, 4), 'AbsTol', 1E-6);
                testCase.verifyEqual(hfdSqrt.pdf(testPoints), hwnd.pdf(testPoints, 4), 'AbsTol', 1E-5);
            end
        end
        
        function testFromDistribution2D(testCase)
            % Test that from Distribution and fromFunction result in equal
            % approximations
            kappa1 = 0.3;
            kappa2 = 1.5;
            lambda = 0.5;
            coeffs = [5, 7];
            tvm = ToroidalVMSineDistribution([1; 2], [kappa1; kappa2], lambda);
            hfd1id = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'identity');
            hfd1sqrt = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'sqrt');
            hfd2id = HypertoroidalFourierDistribution.fromDistribution(tvm, coeffs, 'identity');
            hfd2sqrt = HypertoroidalFourierDistribution.fromDistribution(tvm, coeffs, 'sqrt');
            testCase.verifyClass(hfd2id, 'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfd2sqrt, 'HypertoroidalFourierDistribution');
            testCase.verifySize(hfd2id.C, coeffs);
            testCase.verifySize(hfd2sqrt.C, coeffs);
            % Verify approximation by validating coefficients
            testCase.verifyEqual(hfd2id.C, hfd1id.C, 'AbsTol', 1E-10);
            testCase.verifyEqual(hfd2sqrt.C, hfd1sqrt.C, 'AbsTol', 1E-10);
        end
        
        function testFromDistribution3D(testCase)
            % Test that error is thrown if dimensionality is wrong
            testCase.verifyError( ...
                @()HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution([1; 1; 1], eye(3)), [15, 15]), ...
                'fromDistribution:invalidObject');
            % Test that it yields same results as fromFunction
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            mu = [3; 5; 2];
            coeffs = [21, 21, 21];
            hwnd = HypertoroidalWNDistribution(mu, C);
            hfd1id = HypertoroidalFourierDistribution.fromFunction(@(x, y, z)reshape(hwnd.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'identity');
            hfd1sqrt = HypertoroidalFourierDistribution.fromFunction(@(x, y, z)reshape(hwnd.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'sqrt');
            hfd2id = HypertoroidalFourierDistribution.fromDistribution(hwnd, coeffs, 'identity');
            hfd2sqrt = HypertoroidalFourierDistribution.fromDistribution(hwnd, coeffs, 'sqrt');
            testCase.verifyClass(hfd2id, 'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfd2sqrt, 'HypertoroidalFourierDistribution');
            testCase.verifySize(hfd2id.C, coeffs);
            testCase.verifySize(hfd2sqrt.C, coeffs);
            testCase.verifyEqual(hfd2id.C, hfd1id.C, 'AbsTol', 1E-10);
            testCase.verifyEqual(hfd2sqrt.C, hfd1sqrt.C, 'AbsTol', 1E-10);
        end
        
        function testFromDistributionHWN(testCase)
            % Verify closed form solution against FFT version
            % 2D case
            mu = [1; 2];
            C = 2 * [1, 0.5; 0.5, 1.3];
            hwn = HypertoroidalWNDistribution(mu, C);
            hfd1id = HypertoroidalFourierDistribution.fromDistribution(hwn, [35, 35], 'identity');
            hfd2id = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(hwn.pdf([x(:)'; y(:)']), size(x)), [35, 35], 'identity');
            testCase.verifyEqual(hfd1id.C, hfd2id.C, 'AbsTol', 1E-7);
            % 3D case
            mu = [1; 2; 4];
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            hwn = HypertoroidalWNDistribution(mu, C);
            hfd1id = HypertoroidalFourierDistribution.fromDistribution(hwn, [19, 19, 19], 'identity');
            hfd2id = HypertoroidalFourierDistribution.fromFunction(@(x, y, z)reshape(hwn.pdf([x(:)'; y(:)'; z(:)']), size(x)), [19, 19, 19], 'identity');
            testCase.verifyEqual(hfd1id.C, hfd2id.C, 'AbsTol', 1E-7);
        end
        
        function testMultiply2DWithoutTruncation(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                [xTest, yTest] = meshgrid(-pi:0.1:3*pi);
            else
                [xTest, yTest] = meshgrid(0:1:2*pi);
            end
            testCase.applyFixture(SuppressedWarningsFixture('pdf:mayNotBeNormalized'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            
            tvm1 = ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5);
            tvm2 = ToroidalVMSineDistribution([1; 4], [0.8; 1.5], 0.2);
            
            hfd1id = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'identity');
            hfd2id = HypertoroidalFourierDistribution.fromDistribution(tvm2, [15, 17], 'identity');
            hfd1sqrt = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'sqrt');
            hfd2sqrt = HypertoroidalFourierDistribution.fromDistribution(tvm2, [15, 17], 'sqrt');
            hfd1log = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'log');
            hfd2log = HypertoroidalFourierDistribution.fromDistribution(tvm2, [15, 17], 'log');
            hfdMultid = hfd1id.multiply(hfd2id, [17, 15]+[15, 17]-1);
            hfdMultsqrt = hfd1sqrt.multiply(hfd2sqrt, [17, 15]+[15, 17]-1);
            
            fixture1 = testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            fixture2 = testCase.applyFixture(SuppressedWarningsFixture('Multiply:NotNormalizing'));
            hfdMultlog = hfd1log.multiply(hfd2log, [17, 17]);
            fixture1.teardown;
            fixture2.teardown;
            
            testCase.verifyClass(hfdMultid, 'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfdMultsqrt, 'HypertoroidalFourierDistribution');
            testCase.verifyClass(hfdMultlog, 'HypertoroidalFourierDistribution');
            % Verify that they integrate to 1 (except log)
            testCase.verifyEqual( ...
                integral2(@(x, y)reshape(hfdMultid.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi), 1, 'AbsTol', 1E-6);
            testCase.verifyEqual( ...
                integral2(@(x, y)reshape(hfdMultsqrt.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi), 1, 'AbsTol', 1E-6);
            % Normalize multiplication using integral to obtain approximate
            % ground truth
            normConst = integral2(@(x, y)reshape(tvm1.pdf([x(:)'; y(:)']).*tvm2.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi);
            valTrueApprox = reshape(tvm1.pdf([xTest(:)'; yTest(:)']).*tvm2.pdf([xTest(:)'; yTest(:)']), size(xTest)) / normConst;
            valFourierId = reshape(hfdMultid.pdf([xTest(:)'; yTest(:)']), size(xTest));
            valFourierSqrt = reshape(hfdMultsqrt.pdf([xTest(:)'; yTest(:)']), size(xTest));
            valFourierLog = reshape(hfdMultlog.pdf([xTest(:)'; yTest(:)']), size(xTest));
            % Verify correct function values
            testCase.verifyEqual(valFourierId, valTrueApprox, 'AbsTol', 1E-6);
            testCase.verifyEqual(valFourierSqrt, valTrueApprox, 'AbsTol', 1E-7);
            normConstLog = 1 / integral2(@(x, y)reshape(hfdMultlog.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi);
            testCase.verifyEqual(valFourierLog*normConstLog, valTrueApprox, 'AbsTol', 1E-6);
        end
        
        function testMultiply3DWithoutTruncation(testCase)
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            mu = [1; 2; 3];
            hwnd1 = HypertoroidalWNDistribution(mu, C);
            C = [1.6, 0.8, 1.3; 0.8, 0.7, 0.6; 1.3, 0.6, 1.1];
            hwnd2 = HypertoroidalWNDistribution(mu, C);
            coeffs = 19 * ones(1, 3);
            
            hfd1id = HypertoroidalFourierDistribution.fromFunction( ...
                @(x, y, z)reshape(hwnd1.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'identity');
            hfd2id = HypertoroidalFourierDistribution.fromFunction( ...
                @(x, y, z)reshape(hwnd2.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'identity');
            hfdResultId = hfd1id.multiply(hfd2id, 2*coeffs-1);
            
            [x, y, z] = meshgrid(linspace(-4*pi, pi, 5));
            testPoints = [x(:)'; y(:)'; z(:)'];
            % Calculate approximate normalization factor because integral3
            % is expensive.
            approxNormFactorId = hfdResultId.pdf([1; 1; 1]) / (hfd1id.pdf([1; 1; 1]) * hfd2id.pdf([1; 1; 1]));
            testCase.verifyEqual(approxNormFactorId*(hfd1id.pdf(testPoints) .* hfd2id.pdf(testPoints)), hfdResultId.pdf(testPoints), 'AbsTol', 1E-10);
            
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                hfd1sqrt = HypertoroidalFourierDistribution.fromFunction( ...
                    @(x, y, z)reshape(hwnd1.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'sqrt');
                hfd2sqrt = HypertoroidalFourierDistribution.fromFunction( ...
                    @(x, y, z)reshape(hwnd2.pdf([x(:)'; y(:)'; z(:)']), size(x)), coeffs, 'sqrt');
                hfdResultSqrt = hfd1sqrt.multiply(hfd2sqrt, 2*coeffs-1);
                approxNormFactorSqrt = hfdResultSqrt.pdf([1; 2; 1]) / (hfd1sqrt.pdf([1; 2; 1]) * hfd2sqrt.pdf([1; 2; 1]));
                testCase.verifyEqual(approxNormFactorSqrt*(hfd1sqrt.pdf(testPoints) .* hfd2sqrt.pdf(testPoints)), hfdResultSqrt.pdf(testPoints), 'AbsTol', 1E-10);
            end
        end
        
        function testConvolve2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            tvm1 = ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5);
            tvm2 = ToroidalVMSineDistribution([1; 4], [0.8; 1.5], 0.2);
            hfd1id = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'identity');
            hfdtid = HypertoroidalFourierDistribution.fromDistribution(tvm2, [15, 17], 'identity');
            hfd1sqrt = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'sqrt');
            hfdtsqrt = HypertoroidalFourierDistribution.fromDistribution(tvm2, [15, 17], 'sqrt');
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:cannotTest'));
            hfd1log = HypertoroidalFourierDistribution.fromDistribution(tvm1, [17, 15], 'log');
            hfdtlog = HypertoroidalFourierDistribution.fromDistribution(tvm2, [15, 17], 'log');
            testCase.verifyError(@()hfd1log.convolve(hfdtlog, [15, 13]), 'transformation:unrecognizedTransformation');
            fixture.teardown;
            
            hfd2id = hfd1id.convolve(hfdtid, [15, 13]);
            hfd2sqrt = hfd1sqrt.convolve(hfdtsqrt, [15, 13]);
            testCase.verifyClass(hfd2id, 'HypertoroidalFourierDistribution')
            testCase.verifyClass(hfd2sqrt, 'HypertoroidalFourierDistribution')
            testCase.verifySize(hfd2id.C, [15, 13])
            testCase.verifySize(hfd2sqrt.C, [15, 13])
            hfd2id = hfd1id.convolve(hfdtid);
            hfd2sqrt = hfd1sqrt.convolve(hfdtsqrt);
            testCase.verifySize(hfd2id.C, [17, 15]);
            testCase.verifySize(hfd2sqrt.C, [17, 15]);
            % Verify that they integrate to 1
            testCase.verifyEqual( ...
                integral2(@(x, y)reshape(hfd2id.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi), 1, 'AbsTol', 1E-6);
            testCase.verifyEqual( ...
                integral2(@(x, y)reshape(hfd2sqrt.pdf([x(:)'; y(:)']), size(x)), 0, 2*pi, 0, 2*pi), 1, 'AbsTol', 1E-6);
            [xTest, yTest] = meshgrid(0:2*pi/100:2*pi-2*pi/100);
            valFourierId = reshape(hfd2id.pdf([xTest(:)'; yTest(:)']), size(xTest));
            valFourierSqrt = reshape(hfd2sqrt.pdf([xTest(:)'; yTest(:)']), size(xTest));
            % Calculate cyclic discrete convolution, the result is not
            % normalized, so use values of Fourier for normalization
            % Cyclic convolution would be easier with FFT & IFFT but since
            % is approach is used in FourierDistribution, a different
            % approach is used for validation
            tvm1vals = reshape(tvm1.pdf([xTest(:)'; yTest(:)']), size(xTest));
            tvm2valspadded = padarray(reshape(tvm2.pdf([xTest(:)'; yTest(:)']), size(xTest)), size(tvm1vals), 'circular');
            convolutionTmp = conv2(tvm1vals, tvm2valspadded);
            valConvUnnorm = convolutionTmp(size(tvm1vals, 1)+1:2*size(tvm1vals, 1), size(tvm1vals, 2)+1:2*size(tvm1vals, 2));
            testCase.verifyEqual(valConvUnnorm/sum(valConvUnnorm,[1,2])*sum(valFourierId,[1,2]), valFourierId, 'AbsTol', 1E-6);
            testCase.verifyEqual(valConvUnnorm/sum(valConvUnnorm,[1,2])*sum(valFourierSqrt,[1,2]), valFourierSqrt, 'AbsTol', 1E-6);
        end
        function testShift(testCase)
            Ccell = {2 * [1, 0.5; 0.5, 1], [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1], ...
                [0.7, 0.4, 0.2, -0.5; 0.4, 0.6, 0.1, 0; 0.2, 0.1, 1, -0.3; -0.5, 0, -0.3, 0.9] * 2};
            offsets = {4 * pi * rand(2, 1) - pi, 4 * pi * rand(3, 1) - pi, 4 * pi * rand(4, 1) - pi};
            coeffsPerDim = 13;
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                maxDim = 4;
            else
                maxDim = 3;
            end
            for dim = 2:maxDim
                hfdId = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(zeros(dim, 1), Ccell{dim-1}), coeffsPerDim*ones(1, dim), 'identity');
                hfdSqrt = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(zeros(dim, 1), Ccell{dim-1}), coeffsPerDim*ones(1, dim), 'sqrt');
                hfdIdShiftWN = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(offsets{dim-1}, Ccell{dim-1}), coeffsPerDim*ones(1, dim), 'identity');
                hfdSqrtShiftWN = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(offsets{dim-1}, Ccell{dim-1}), coeffsPerDim*ones(1, dim), 'sqrt');
                hfdIdShiftFD = hfdId.shift(offsets{dim-1});
                hfdSqrtShiftFD = hfdSqrt.shift(offsets{dim-1});
                testCase.verifyEqual(hfdIdShiftFD.C, hfdIdShiftWN.C, 'AbsTol', abs(max(hfdIdShiftWN.C(:))/1000));
                testCase.verifyEqual(hfdSqrtShiftFD.C, hfdSqrtShiftWN.C, 'AbsTol', abs(max(hfdSqrtShiftWN.C(:))/1000));
            end
        end
        function testIntegral2D(testCase)
            % Test against implementation in toroidal (test case that this
            % works correctly exists in ToroidalFourierDistributionTest)
            kappa1 = 0.3;
            kappa2 = 1.5;
            lambda = 0.5;
            coeffs = [5, 7];
            tvm = ToroidalVMSineDistribution([1; 2], [kappa1; kappa2], lambda);
            hfdId = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'sqrt');
            
            tfdId = ToroidalFourierDistribution(hfdId.C, hfdId.transformation);
            tfdSqrt = ToroidalFourierDistribution(hfdSqrt.C, hfdSqrt.transformation);
            testCase.verifyEqual(hfdId.integral([0; 0], [pi; pi]), tfdId.integral([0; 0], [pi; pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0; 0], [pi; pi]), tfdSqrt.integral([0; 0], [pi; pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdId.integral([0; 0], [pi; 2 * pi]), tfdId.integral([0; 0], [pi; 2 * pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0; 0], [pi; 2 * pi]), tfdSqrt.integral([0; 0], [pi; 2 * pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdId.integral([0; -1], [3 * pi; 5 * pi]), tfdId.integral([0; -1], [3 * pi; 5 * pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0; -1], [3 * pi; 5 * pi]), tfdSqrt.integral([0; -1], [3 * pi; 5 * pi]), 'AbsTol', 1E-4);
        end
        function testIntegral2DUnnorm(testCase)
            % Test against implementation in toroidal (test case that this
            % works correctly exists in ToroidalFourierDistributionTest)
            kappa1 = 0.3;
            kappa2 = 1.5;
            lambda = 0.5;
            coeffs = [5, 7];
            tvm = ToroidalVMSineDistribution([1; 2], [kappa1; kappa2], lambda);
            hfdId = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), coeffs, 'sqrt');
            
            
            tfdId = ToroidalFourierDistribution(hfdId.C, hfdId.transformation);
            tfdSqrt = ToroidalFourierDistribution(hfdSqrt.C, hfdSqrt.transformation);
            testCase.verifyEqual(hfdId.integral([0; 0], [pi; pi]), tfdId.integral([0; 0], [pi; pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0; 0], [pi; pi]), tfdSqrt.integral([0; 0], [pi; pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdId.integral([0; 0], [pi; 2 * pi]), tfdId.integral([0; 0], [pi; 2 * pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0; 0], [pi; 2 * pi]), tfdSqrt.integral([0; 0], [pi; 2 * pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdId.integral([0; -1], [3 * pi; 5 * pi]), tfdId.integral([0; -1], [3 * pi; 5 * pi]), 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.integral([0; -1], [3 * pi; 5 * pi]), tfdSqrt.integral([0; -1], [3 * pi; 5 * pi]), 'AbsTol', 1E-4);
        end
        function testIntegral3D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test against integral
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            mu = [1; 2; 5];
            coeffs = [5, 9, 7];
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized')); % Very coarse approximation, lack of normalization is to be expected
            hfdId = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu, C), coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu, C), coeffs, 'sqrt');
            fixture.teardown;
            
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true;end
            if enableExpensive
                testCase.verifyEqual(hfdId.integral([0; 0; 0], [pi; pi; pi]), integral3(@(x, y, z)reshape(hfdId.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, pi, 0, pi, 0, pi), 'AbsTol', 1E-4);
                testCase.verifyEqual(hfdSqrt.integral([0; 0; 0], [pi; pi; pi]), integral3(@(x, y, z)reshape(hfdSqrt.pdf([x(:)'; y(:)'; z(:)']), size(x)), 0, pi, 0, pi, 0, pi), 'AbsTol', 1E-4);
            end
        end
        function testCovariance2dimd2D(testCase)
            % Test against integral
            C = [0.7, 0.4; 0.4, 0.6];
            mu = [1; 2];
            coeffs = [15, 17];
            twn = ToroidalWNDistribution(mu, C);
            hfdId = HypertoroidalFourierDistribution.fromDistribution(twn, coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromDistribution(twn, coeffs, 'sqrt');
            
            testCase.verifyEqual(hfdId.covariance2dimD, twn.covariance4D, 'AbsTol', 1E-4);
            testCase.verifyEqual(hfdSqrt.covariance2dimD, twn.covariance4D, 'AbsTol', 1E-4);
        end
        function testCovariance2dimd3D(testCase)
            C = [0.7, 0.4, 0.2; 0.4, 0.6, 0.1; 0.2, 0.1, 1];
            mu = [1; 2; 5];
            coeffs = [15, 19, 17];
            hfdId = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu, C), coeffs, 'identity');
            hfdSqrt = HypertoroidalFourierDistribution.fromDistribution(HypertoroidalWNDistribution(mu, C), coeffs, 'sqrt');
            
            % Test by taking only two correlated dimensions and comparing
            % to that result for each pair of dimensions
            Cov6DId = hfdId.covariance2dimD;
            Cov6DSqrt = hfdSqrt.covariance2dimD;
            twn = ToroidalWNDistribution(mu(1:2), C(1:2, 1:2));
            testCase.verifyEqual(Cov6DId(1:4, 1:4), twn.covariance4D, 'AbsTol', 1E-4);
            testCase.verifyEqual(Cov6DSqrt(1:4, 1:4), twn.covariance4D, 'AbsTol', 1E-4);
            twn = ToroidalWNDistribution(mu([1, 3]), C([1, 3], [1, 3]));
            testCase.verifyEqual(Cov6DId([1, 2, 5, 6], [1, 2, 5, 6]), twn.covariance4D, 'AbsTol', 1E-4);
            testCase.verifyEqual(Cov6DSqrt([1, 2, 5, 6], [1, 2, 5, 6]), twn.covariance4D, 'AbsTol', 1E-4);
            twn = ToroidalWNDistribution(mu(2:3), C(2:3, 2:3));
            testCase.verifyEqual(Cov6DId(3:6, 3:6), twn.covariance4D, 'AbsTol', 1E-4);
            testCase.verifyEqual(Cov6DSqrt(3:6, 3:6), twn.covariance4D, 'AbsTol', 1E-4);
        end
        function testMarginalize2Dto1Dtvm(testCase)
            tvm = ToroidalVMSineDistribution([1; 2], [2; 3], 3);
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(tvm, [53, 53], [transformation{:}]);
                for d = 1:2
                    vm = tvm.marginalizeTo1D(d);
                    fd1 = hfd.marginalizeTo1D(d);
                    fd2 = hfd.marginalizeOut(1+(d == 1));
                    testCase.verifyEqual(vm.hellingerDistanceNumerical(fd1), 0, 'AbsTol', 1E-3);
                    testCase.verifyEqual(vm.hellingerDistanceNumerical(fd2), 0, 'AbsTol', 1E-3);
                end
            end
        end
        function testMarginalize2Dto1Dtwn(testCase)
            twn = ToroidalWNDistribution([3;4],2*[1,0.8;0.8,1]);
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(twn, [71, 53], [transformation{:}]);
                grid = linspace(-pi,3*pi,300);

                hfdMarginalized = hfd.marginalizeOut(2);
                marginlized1DViaIntegral=@(x)arrayfun(@(xCurr)integral(@(y)reshape(twn.pdf([xCurr*ones(1,size(y,2));y(:)']),size(y)),0,2*pi),x);
                testCase.verifyEqual(hfdMarginalized.pdf(grid), marginlized1DViaIntegral(grid), 'AbsTol', 1E-13);

                hfdMarginalized = hfd.marginalizeOut(1);
                marginlized1DViaIntegral=@(y)arrayfun(@(yCurr)integral(@(x)reshape(twn.pdf([x(:)';yCurr*ones(1,size(x,2))]),size(x)),0,2*pi),y);
                testCase.verifyEqual(hfdMarginalized.pdf(grid), marginlized1DViaIntegral(grid), 'AbsTol', 1E-13);
            end
        end
        function testMarginalize3Dto2D(testCase)
            twn = HypertoroidalWNDistribution([3;4;6],2*[1,0.8,0.3;0.8,1,0.5;0.3,0.5,2]);
            for transformation = {'identity', 'sqrt'}
                if strcmp([transformation{:}],'identity')
                    tol = 1e-15;
                else
                    tol = 1e-5;
                end

                hfd = HypertoroidalFourierDistribution.fromDistribution(twn, [41, 53, 53], [transformation{:}]);
                [mesh1,mesh2] = meshgrid(linspace(-pi,3*pi,25));

                hfdMarginalized = hfd.marginalizeOut(3);            
                marginlized1DViaIntegral=@(x,y)arrayfun(@(xCurr,yCurr)integral(@(z)reshape(twn.pdf([xCurr*ones(1,size(z,2));yCurr*ones(1,size(z,2));z(:)']),size(z)),0,2*pi),x,y);
                testCase.verifyEqual(hfdMarginalized.pdf([mesh1(:)';mesh2(:)']), marginlized1DViaIntegral(mesh1(:)',mesh2(:)'), 'AbsTol', tol);

                hfdMarginalized = hfd.marginalizeOut(2);
                marginlized1DViaIntegral=@(x,z)arrayfun(@(xCurr,zCurr)integral(@(y)reshape(twn.pdf([xCurr*ones(1,size(y,2));y(:)';zCurr*ones(1,size(y,2))]),size(y)),0,2*pi),x,z);
                testCase.verifyEqual(hfdMarginalized.pdf([mesh1(:)';mesh2(:)']), marginlized1DViaIntegral(mesh1(:)',mesh2(:)'), 'AbsTol', tol);

                hfdMarginalized = hfd.marginalizeOut(1);
                marginlized1DViaIntegral=@(y,z)arrayfun(@(yCurr,zCurr)integral(@(x)reshape(twn.pdf([x(:)';yCurr*ones(1,size(x,2));zCurr*ones(1,size(x,2))]),size(x)),0,2*pi),y,z);
                testCase.verifyEqual(hfdMarginalized.pdf([mesh1(:)';mesh2(:)']), marginlized1DViaIntegral(mesh1(:)',mesh2(:)'), 'AbsTol', tol);
            end
        end
        function testMarginalize3Dto1D(testCase)
            twn = HypertoroidalWNDistribution([3;4;6],2*[1,0.8,0.3;0.8,1,0.5;0.3,0.5,2]);
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(twn, [41, 53, 53], [transformation{:}]);
                grid = linspace(-pi,3*pi,100);

                hfdMarginalized = hfd.marginalizeOut([2,3]); 
                marginlized1DViaIntegral=@(x)arrayfun(@(xCurr)integral2(@(y,z)reshape(twn.pdf([xCurr*ones(1,numel(y));y(:)';z(:)']),size(y)),0,2*pi,0,2*pi),x);
                testCase.verifyEqual(hfdMarginalized.pdf(grid), marginlized1DViaIntegral(grid), 'RelTol', 5E-7);

                hfdMarginalized = hfd.marginalizeOut([1,3]);
                marginlized1DViaIntegral=@(y)arrayfun(@(yCurr)integral2(@(x,z)reshape(twn.pdf([x(:)';yCurr*ones(1,numel(x));z(:)']),size(x)),0,2*pi,0,2*pi),y);
                testCase.verifyEqual(hfdMarginalized.pdf(grid), marginlized1DViaIntegral(grid), 'RelTol', 5E-7);

                hfdMarginalized = hfd.marginalizeOut([1,2]);
                marginlized1DViaIntegral=@(z)arrayfun(@(zCurr)integral2(@(x,y)reshape(twn.pdf([x(:)';y(:)';zCurr*ones(1,numel(x))]),size(x)),0,2*pi,0,2*pi),z);
                testCase.verifyEqual(hfdMarginalized.pdf(grid), marginlized1DViaIntegral(grid), 'RelTol', 5E-7);
            end
        end
        function testPlotting(testCase)
            dist = ToroidalVMSineDistribution([1; 2], [2; 3], 3);
            tfd = ToroidalFourierDistribution.fromDistribution(dist, [51, 201], 'identity');
            figure(987)
            h = tfd.plot;
            noPoints = 102;
            [alpha, beta] = meshgrid(linspace(0, 2*pi, noPoints));
            fvals = reshape(tfd.pdf([alpha(:)'; beta(:)']), size(alpha));
            testCase.verifyEqual(get(h, 'ZData'), fvals, 'AbsTol', 1E-15);
            close(987)
        end
        function testConv1D(testCase)
            for transformation = {'identity', 'sqrt'}
                fd1 = FourierDistribution.fromDistribution(VMDistribution(0, 1), 13, [transformation{:}]);
                fd2 = FourierDistribution.fromDistribution(VMDistribution(2, 1), 13, [transformation{:}]);
                hfd1 = HypertoroidalFourierDistribution(fd1.c', [transformation{:}]);
                hfd2 = HypertoroidalFourierDistribution(fd2.c', [transformation{:}]);
                fdConv = fd1.convolve(fd2);
                hfdConv = hfd1.convolve(hfd2);
                testCase.verifyEqual(hfdConv.C, fdConv.c', 'AbsTol', 1E-8);
            end
        end
        function testMult1D(testCase)
            for transformation = {'identity', 'sqrt'}
                fd1 = FourierDistribution.fromDistribution(VMDistribution(0, 1), 13, [transformation{:}]);
                fd2 = FourierDistribution.fromDistribution(VMDistribution(2, 1), 13, [transformation{:}]);
                hfd1 = HypertoroidalFourierDistribution(fd1.c', [transformation{:}]);
                hfd2 = HypertoroidalFourierDistribution(fd2.c', [transformation{:}]);
                fdMult = fd1.multiply(fd2, 2*(numel(fd1.a) + numel(fd1.b))-1);
                hfdMult = hfd1.multiply(hfd2, 2*size(hfd1.C, 1)-1);
                testCase.verifyEqual(hfdMult.C, fdMult.c', 'AbsTol', 1E-10);
            end
        end
        
        function testApproxFromEvenGrid1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            hwn = WNDistribution(1,0.5);
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized')); % Very few coefficients, expect unnormalized
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,4);
            %% Without square root, we get the same values agian. The reason is that the Fourier series goes through these points and appending it more coefficients will not change it
            hfdId = HypertoroidalFourierDistribution.fromFunctionValues(reshape(hgd.gridValues,[hgd.noOfGridPoints,1]),5,'identity');
            testCase.verifyEqual(hfdId.pdf(hgd.getGrid())',hgd.gridValues,'AbsTol',1E-15);
            %% Values are not precisely matched!
            % The reason is that while the non-padded Fourier series would
            % have gone through these points, the padded one is changed
            % because the normalization.
            hfdSqrt = HypertoroidalFourierDistribution.fromFunctionValues(reshape(hgd.gridValues,[hgd.noOfGridPoints,1]),5,'sqrt');
            testCase.verifyEqual(hfdSqrt.pdf(hgd.getGrid())',hgd.gridValues,'AbsTol',0.05);
            
        end
        
        function testApproxFromEvenGrid2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            hwn = HypertoroidalWNDistribution([1;1],0.1*eye(2));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized')); % Very few coefficients, expect unnormalized
            hgd = HypertoroidalGridDistribution.fromDistribution(hwn,[4,4]);
            %% Without square root, we get the same values agian. See above.
            hfdId = HypertoroidalFourierDistribution.fromFunctionValues(reshape(hgd.gridValues,hgd.noOfGridPoints),[5,5],'identity');
            testCase.verifyEqual(hfdId.pdf(hgd.getGrid())',hgd.gridValues,'AbsTol',1E-16);
            %% Values are not precisely matched! See above.
            hfdSqrt = HypertoroidalFourierDistribution.fromFunctionValues(reshape(hgd.gridValues,hgd.noOfGridPoints),[5,5],'sqrt');
            testCase.verifyEqual(hfdSqrt.pdf(hgd.getGrid())',hgd.gridValues,'AbsTol',0.1);
        end
        
        function testMatchFunctionOdd1D(testCase)
            fNeeds6 = @(x)1/(2*pi)+0.05*cos(x)+0.05*cos(2*x)+0.05*cos(3*x)+0.05*sin(x)-0.05*sin(2*x);
            fNeeds7 = @(x)1/(2*pi)+0.05*cos(x)+0.05*cos(2*x)+0.05*cos(3*x)+0.05*sin(x)-0.05*sin(2*x)+0.05*sin(3*x);
            cdNeeds6 = CustomCircularDistribution(fNeeds6);
            cdNeeds7 = CustomCircularDistribution(fNeeds7);

            hgdNeeds6Has6 = HypertoroidalGridDistribution.fromDistribution(cdNeeds6,6);
            hgdNeeds6Has7 = HypertoroidalGridDistribution.fromDistribution(cdNeeds6,7);
            hgdNeeds7Has6 = HypertoroidalGridDistribution.fromDistribution(cdNeeds7,6);
            hgdNeeds7Has7 = HypertoroidalGridDistribution.fromDistribution(cdNeeds7,7);

            hfdNeeds6BasedOn6 = HypertoroidalFourierDistribution.fromDistribution(hgdNeeds6Has6,7,'identity');
            hfdNeeds6BasedOn7 = HypertoroidalFourierDistribution.fromDistribution(hgdNeeds6Has7,7,'identity');
            hfdNeeds7BasedOn6 = HypertoroidalFourierDistribution.fromDistribution(hgdNeeds7Has6,7,'identity');
            hfdNeeds7BasedOn7 = HypertoroidalFourierDistribution.fromDistribution(hgdNeeds7Has7,7,'identity');
            
            % If only 6 coefficients are required and we convert a
            % GridDistribution with 6, it should be perfect
            testCase.verifyEqual(cdNeeds6.l2distanceCdfNumerical(hfdNeeds6BasedOn6),0,'AbsTol',1E-30);
            % Same for 7
            testCase.verifyEqual(cdNeeds7.l2distanceCdfNumerical(hfdNeeds7BasedOn7),0,'AbsTol',1E-30);
            
            % If we only need 6, then it must not matter if we use 6 or 7
            % grid points
            testCase.verifyEqual(hfdNeeds6BasedOn6.C,hfdNeeds6BasedOn7.C,'AbsTol',1E-16);
            
            % If we need 7 and only use 6, then the last coefficient only
            % influences the complex part and all real parts should be fine
            testCase.verifyEqual(real(hfdNeeds7BasedOn6.C),real(hfdNeeds7BasedOn7.C),'AbsTol',1E-16);
            % Further, the middle ones should be perfect also in the
            % complex part
            testCase.verifyEqual(hfdNeeds7BasedOn6.C(2:end-1),hfdNeeds7BasedOn7.C(2:end-1),'AbsTol',1E-16);
            % There should be no imaginary part in the outermost
            % coefficients if only 6 grid points were used
            testCase.verifyEqual(imag(hfdNeeds7BasedOn6.C([1,end])),[0;0],'AbsTol',1E-30);
        end
        
        function testConditioning2D(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            grid = linspace(0,2*pi,100);
            z1 = 6;
            z2 = 5;
            cdNormDim1 = chd.conditionOn(1,z1);
            cdNormDim2 = chd.conditionOn(2,z2);
            
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(twn,[301,301],[transformation{:}]);
                % Use IFFT
                ffdNormDim1 = hfd.conditionOn(1,z1,false);
                ffdNormDim2 = hfd.conditionOn(2,z2,false);
                testCase.verifyEqual(ffdNormDim1.integral(),1);
                testCase.verifyEqual(ffdNormDim2.integral(),1);
                testCase.verifyEqual(ffdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'RelTol',8e-11);
                testCase.verifyEqual(ffdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'RelTol',8e-11);
                % Use FFTN
                ffdNormDim1 = hfd.conditionOn(1,z1,true);
                ffdNormDim2 = hfd.conditionOn(2,z2,true);
                testCase.verifyEqual(ffdNormDim1.integral(),1);
                testCase.verifyEqual(ffdNormDim2.integral(),1);
                testCase.verifyEqual(ffdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'RelTol',8e-11);
                testCase.verifyEqual(ffdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'RelTol',8e-11);
            end
        end
        
        function testSliceAt2D(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            grid = linspace(0,2*pi,100);
            z1 = 6;
            z2 = 5;
            cdSliceDim1 = chd.sliceAt(1,z1);
            cdSliceDim2 = chd.sliceAt(2,z2);
            
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(twn,[301,301],[transformation{:}]);
                ffdSliceDim1IFFT = hfd.sliceAt(1,z1,false);
                ffdSliceDim2IFTT = hfd.sliceAt(2,z2,false);
                ffdSliceDim1FFTN = hfd.sliceAt(1,z1,true);
                ffdSliceDim2FFTN = hfd.sliceAt(2,z2,true);
                
                testCase.verifyEqual(ffdSliceDim1IFFT.pdf(grid),cdSliceDim1.pdf(grid),'RelTol',8e-11);
                testCase.verifyEqual(ffdSliceDim2IFTT.pdf(grid),cdSliceDim2.pdf(grid),'RelTol',8e-11);
                testCase.verifyEqual(ffdSliceDim1FFTN.pdf(grid),cdSliceDim1.pdf(grid),'RelTol',8e-11);
                testCase.verifyEqual(ffdSliceDim2FFTN.pdf(grid),cdSliceDim2.pdf(grid),'RelTol',8e-11);
            end
        end
        
        function testConditioning3DTo2D(testCase)
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
            
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(hwn,[101,101,101],[transformation{:}]);
                % Use IFFT
                ffdNormDim1 = hfd.conditionOn(1,z1,false);
                ffdNormDim2 = hfd.conditionOn(2,z2,false);
                ffdNormDim3 = hfd.conditionOn(3,z3,false);
                testCase.verifyEqual(ffdNormDim1.integral(),1);
                testCase.verifyEqual(ffdNormDim2.integral(),1);
                testCase.verifyEqual(ffdNormDim3.integral(),1);
                testCase.verifyEqual(ffdNormDim1.pdf([mesh1(:)';mesh2(:)']),cdNormDim1.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdNormDim2.pdf([mesh1(:)';mesh2(:)']),cdNormDim2.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdNormDim3.pdf([mesh1(:)';mesh2(:)']),cdNormDim3.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                % Use FFTN
                ffdNormDim1 = hfd.conditionOn(1,z1,true);
                ffdNormDim2 = hfd.conditionOn(2,z2,true);
                ffdNormDim3 = hfd.conditionOn(3,z3,true);
                testCase.verifyEqual(ffdNormDim1.integral(),1);
                testCase.verifyEqual(ffdNormDim2.integral(),1);
                testCase.verifyEqual(ffdNormDim3.integral(),1);
                testCase.verifyEqual(ffdNormDim1.pdf([mesh1(:)';mesh2(:)']),cdNormDim1.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdNormDim2.pdf([mesh1(:)';mesh2(:)']),cdNormDim2.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdNormDim3.pdf([mesh1(:)';mesh2(:)']),cdNormDim3.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
            end
        end
        
        function testSliceAt3DTo2D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            
            z1 = 1;
            z2 = 3;
            z3 = 4;
            cdNormDim1 = chd.sliceAt(1,z1);
            cdNormDim2 = chd.sliceAt(2,z2);
            cdNormDim3 = chd.sliceAt(3,z3);
            
            [mesh1,mesh2] = meshgrid(linspace(0,2*pi,20));
            
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(hwn,[101,101,101],[transformation{:}]);
                % Use IFFT
                ffdSliceDim1 = hfd.sliceAt(1,z1,false);
                ffdSliceDim2 = hfd.sliceAt(2,z2,false);
                ffdSliceDim3 = hfd.sliceAt(3,z3,false);
                testCase.verifyEqual(ffdSliceDim1.pdf([mesh1(:)';mesh2(:)']),cdNormDim1.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdSliceDim2.pdf([mesh1(:)';mesh2(:)']),cdNormDim2.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdSliceDim3.pdf([mesh1(:)';mesh2(:)']),cdNormDim3.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                % Use FFTN
                ffdSliceDim1 = hfd.sliceAt(1,z1,true);
                ffdSliceDim2 = hfd.sliceAt(2,z2,true);
                ffdSliceDim3 = hfd.sliceAt(3,z3,true);
                testCase.verifyEqual(ffdSliceDim1.pdf([mesh1(:)';mesh2(:)']),cdNormDim1.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdSliceDim2.pdf([mesh1(:)';mesh2(:)']),cdNormDim2.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
                testCase.verifyEqual(ffdSliceDim3.pdf([mesh1(:)';mesh2(:)']),cdNormDim3.pdf([mesh1(:)';mesh2(:)']),'AbsTol',1e-9);
            end
        end
        
        function testConditioning3DTo1D(testCase)
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
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(hwn,[101,101,101],[transformation{:}]);
                if strcmp([transformation{:}],'identity')
                    tol=5e-15;
                else
                    tol=5e-9;
                end
                % Use IFFT
                ffdNormDim1 = hfd.conditionOn([1,2],z1,false);
                ffdNormDim2 = hfd.conditionOn([1,3],z2,false);
                ffdNormDim3 = hfd.conditionOn([2,3],z3,false);

                testCase.verifyEqual(ffdNormDim1.integral(),1);
                testCase.verifyEqual(ffdNormDim2.integral(),1);
                testCase.verifyEqual(ffdNormDim3.integral(),1);
                testCase.verifyEqual(ffdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim3.pdf(grid),cdNormDim3.pdf(grid),'AbsTol',tol);
                % Use FFTN
                ffdNormDim1 = hfd.conditionOn([1,2],z1,true);
                ffdNormDim2 = hfd.conditionOn([1,3],z2,true);
                ffdNormDim3 = hfd.conditionOn([2,3],z3,true);

                testCase.verifyEqual(ffdNormDim1.integral(),1);
                testCase.verifyEqual(ffdNormDim2.integral(),1);
                testCase.verifyEqual(ffdNormDim3.integral(),1);
                testCase.verifyEqual(ffdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim3.pdf(grid),cdNormDim3.pdf(grid),'AbsTol',tol);
            end
        end
        
        function testSliceAt3DTo1D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution(@(xs)hwn.pdf(xs),hwn.dim);
            
            z1 = [1;5];
            z2 = [2;4];
            z3 = [3;6];
            cdNormDim1 = chd.sliceAt([1,2],z1);
            cdNormDim2 = chd.sliceAt([1,3],z2);
            cdNormDim3 = chd.sliceAt([2,3],z3);
            
            grid = linspace(0,2*pi,100);
            for transformation = {'identity', 'sqrt'}
                hfd = HypertoroidalFourierDistribution.fromDistribution(hwn,[101,101,101],[transformation{:}]);
                if strcmp([transformation{:}],'identity')
                    tol=5e-15;
                else
                    tol=5e-9;
                end
                % Use IFFT
                ffdNormDim1 = hfd.sliceAt([1,2],z1,false);
                ffdNormDim2 = hfd.sliceAt([1,3],z2,false);
                ffdNormDim3 = hfd.sliceAt([2,3],z3,false);

                testCase.verifyEqual(ffdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim3.pdf(grid),cdNormDim3.pdf(grid),'AbsTol',tol);
                % Use FFTN
                ffdNormDim1 = hfd.sliceAt([1,2],z1,true);
                ffdNormDim2 = hfd.sliceAt([1,3],z2,true);
                ffdNormDim3 = hfd.sliceAt([2,3],z3,true);

                testCase.verifyEqual(ffdNormDim1.pdf(grid),cdNormDim1.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim2.pdf(grid),cdNormDim2.pdf(grid),'AbsTol',tol);
                testCase.verifyEqual(ffdNormDim3.pdf(grid),cdNormDim3.pdf(grid),'AbsTol',tol);
            end
        end
        
        function testLikelihood2DTo1D(testCase)
            twn = ToroidalWNDistribution([3;4],[1,0.8;0.8,1]);
            chd = CustomHypertoroidalDistribution.fromDistribution(twn);
            xFix = 1;
            yFix = 3;
            
            likelihoodDim1chd = chd.likelihood(1,xFix);
            likelihoodDim2chd = chd.likelihood(2,yFix);
            
            for transformation = {'identity', 'sqrt'}
                hgd = HypertoroidalFourierDistribution.fromDistribution(twn,[101,101],[transformation{:}]);
                likelihoodDim1hfd = hgd.likelihood(1,xFix);
                likelihoodDim2hfd = hgd.likelihood(2,yFix);
                grid = linspace(-pi,3*pi,100);
                testCase.verifyEqual(likelihoodDim1hfd.pdf(grid),likelihoodDim1chd.pdf(grid),'RelTol',5e-6);
                testCase.verifyEqual(likelihoodDim2hfd.pdf(grid),likelihoodDim2chd.pdf(grid),'RelTol',5e-6);
            end
        end
        
        function testLikelihood3DTo2D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            chd = CustomHypertoroidalDistribution.fromDistribution(hwn);
            xFix = 1;
            yFix = 3;
            zFix = 4;
            
            likelihoodDim1chd = chd.likelihood(1,xFix);
            likelihoodDim2chd = chd.likelihood(2,yFix);
            likelihoodDim3chd = chd.likelihood(3,zFix);
            
            [mesh1,mesh2] = meshgrid(linspace(0,2*pi,20));
            for transformation = {'identity', 'sqrt'}
                hgd = HypertoroidalFourierDistribution.fromDistribution(hwn,[101,101,101]);
                likelihoodDim1hfd = hgd.likelihood(1,xFix);
                likelihoodDim2hfd = hgd.likelihood(2,yFix);
                likelihoodDim3hfd = hgd.likelihood(3,zFix);
                testCase.verifyEqual(likelihoodDim1hfd.pdf([mesh1(:)';mesh2(:)']),likelihoodDim1chd.pdf([mesh1(:)';mesh2(:)']),'RelTol',5e-5);
                testCase.verifyEqual(likelihoodDim2hfd.pdf([mesh1(:)';mesh2(:)']),likelihoodDim2chd.pdf([mesh1(:)';mesh2(:)']),'RelTol',5e-5);
                testCase.verifyEqual(likelihoodDim3hfd.pdf([mesh1(:)';mesh2(:)']),likelihoodDim3chd.pdf([mesh1(:)';mesh2(:)']),'RelTol',5e-5);
            end
        end
        
        function testLikelihood3DTo1D(testCase)
            C = [0.7,0.4,0.2;0.4,0.6,0.1;0.2,0.1,1];
            hwn = HypertoroidalWNDistribution([2;4;6],C);
            hgd = HypertoroidalFourierDistribution.fromDistribution(hwn,[101,101,101]);
            chd = CustomHypertoroidalDistribution.fromDistribution(hwn);
            xyFix = [1;5];
            xzFix = [2;4];
            yzFix = [3;6];
            
            grid = linspace(-pi,3*pi,100);
            
            likelihoodDim12chd = chd.likelihood([1,2],xyFix);
            likelihoodDim13chd = chd.likelihood([1,3],xzFix);
            likelihoodDim23chd = chd.likelihood([2,3],yzFix);
            
            for transformation = {'identity', 'sqrt'}
                likelihoodDim12hfd = hgd.likelihood([1,2],xyFix);
                likelihoodDim13hfd = hgd.likelihood([1,3],xzFix);
                likelihoodDim23hfd = hgd.likelihood([2,3],yzFix);

                testCase.verifyEqual(likelihoodDim12hfd.pdf(grid),likelihoodDim12chd.pdf(grid),'RelTol',5e-5);
                testCase.verifyEqual(likelihoodDim13hfd.pdf(grid),likelihoodDim13chd.pdf(grid),'RelTol',5e-5);
                testCase.verifyEqual(likelihoodDim23hfd.pdf(grid),likelihoodDim23chd.pdf(grid),'RelTol',5e-5);
            end
        end
    end
end

classdef SphericalHarmonicsDistributionComplexTest < matlab.unittest.TestCase
    methods(Test)
        function testNormalization(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            testCase.verifyWarning(@()SphericalHarmonicsDistributionComplex(1), 'Normalization:notNormalized');
            testCase.verifyError(@()SphericalHarmonicsDistributionComplex(0), 'Normalization:almostZero');
            coeffRand = rand(1, 9);
            unnormalizedCoeffs = [coeffRand(1), NaN, NaN, NaN, NaN; coeffRand(2) + 1i * coeffRand(3), coeffRand(4), -coeffRand(2) + 1i * coeffRand(3), NaN, NaN; ...
                coeffRand(5) + 1i * coeffRand(6), coeffRand(7) + 1i * coeffRand(8), coeffRand(9), -coeffRand(7) + 1i * coeffRand(8), coeffRand(5) - 1i * coeffRand(6)];
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            shd = SphericalHarmonicsDistributionComplex(unnormalizedCoeffs);
            fixture.teardown;
            testCase.verifyEqual(shd.integral, 1, 'AbsTol', 1E-5);
            % Enforce unnormalized coefficients and compare ratio
            [phi, theta] = deal(rand(1, 10)*2*pi, rand(1, 10)*pi-pi/2);
            [x, y, z] = sph2cart(phi, theta, 1);
            valsNormalized = shd.pdf([x; y; z]);
            shd.coeffMat = unnormalizedCoeffs;
            valsUnnormalized = shd.pdf([x; y; z]);
            testCase.verifyEqual(diff(valsNormalized./valsUnnormalized), zeros(1, size(x, 2)-1), 'AbsTol', 1E-6);
        end
        function testIntegralAnalytical(testCase)
            rng(1);
            coeffRand = rand(1, 9);
            unnormalizedCoeffs = [coeffRand(1), NaN, NaN, NaN, NaN; coeffRand(2) + 1i * coeffRand(3), coeffRand(4), -coeffRand(2) + 1i * coeffRand(3), NaN, NaN; ...
                coeffRand(5) + 1i * coeffRand(6), coeffRand(7) + 1i * coeffRand(8), coeffRand(9), -coeffRand(7) + 1i * coeffRand(8), coeffRand(5) - 1i * coeffRand(6)];
            % Initialize with uniform and overwrite to prevent normalization
            shd = SphericalHarmonicsDistributionComplex.fromDistributionViaIntegral(HypersphericalUniformDistribution(3), 1);
            for transformation = {'identity', 'sqrt'}
                shd.coeffMat = unnormalizedCoeffs;
                shd.transformation = [transformation{:}];
                intValNum = shd.integralNumerical;
                intValAna = shd.integralAnalytical;
                testCase.verifyEqual(intValAna, intValNum, 'AbsTol', 1E-5);
            end
        end
        function testTruncation(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            coeffRand = rand(1, 9);
            unnormalizedCoeffs = [coeffRand(1), NaN, NaN, NaN, NaN; coeffRand(2) + 1i * coeffRand(3), coeffRand(4), -coeffRand(2) + 1i * coeffRand(3), NaN, NaN; ...
                coeffRand(5) + 1i * coeffRand(6), coeffRand(7) + 1i * coeffRand(8), coeffRand(9), -coeffRand(7) + 1i * coeffRand(8), coeffRand(5) - 1i * coeffRand(6)];
            cshd = SphericalHarmonicsDistributionComplex(unnormalizedCoeffs);
            fixture.teardown;
            testCase.verifyWarning(@()cshd.truncate(4), 'Truncate:TooFewCoefficients');
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            cshd2 = cshd.truncate(4);
            testCase.verifySize(cshd2.coeffMat, [5, 9]);
            testCase.verifyTrue(all(isnan(cshd2.coeffMat(5, :)) | cshd2.coeffMat(5, :) == 0));
            cshd3 = cshd.truncate(5);
            fixture.teardown;
            testCase.verifySize(cshd3.coeffMat, [6, 11]);
            testCase.verifyTrue(all(isnan(cshd3.coeffMat(5:6, :)) | cshd3.coeffMat(5:6, :) == 0,[1,2]));
            cshd4 = cshd2.truncate(3);
            testCase.verifySize(cshd4.coeffMat, [4, 7]);
            cshd5 = cshd3.truncate(3);
            testCase.verifySize(cshd5.coeffMat, [4, 7]);
            
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi, theta, 1);
            testCase.verifyEqual(cshd2.pdf([x; y; z]), cshd.pdf([x; y; z]), 'AbsTol', 1E-6);
            testCase.verifyEqual(cshd3.pdf([x; y; z]), cshd.pdf([x; y; z]), 'AbsTol', 1E-6);
            testCase.verifyEqual(cshd4.pdf([x; y; z]), cshd.pdf([x; y; z]), 'AbsTol', 1E-6);
            testCase.verifyEqual(cshd5.pdf([x; y; z]), cshd.pdf([x; y; z]), 'AbsTol', 1E-6);
        end
        % Test some basis functions by comparing them with the table given
        % on https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
        function testl0m0(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi)); % Initialize with 1/sqrt(4*pi) and overwrite to prevent normalization
            shd.coeffMat = [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, 0, 0, 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), ones(1, length(x))*sqrt(1/(4 * pi)), 'AbsTol', 1E-6);
        end
        function testl1m0(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 0, 1, 0, NaN, NaN; 0, 0, 0, 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), sqrt(3/(4 * pi))*z, 'AbsTol', 1E-6);
        end
        function testl2m0(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, 0, 1, 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(5/pi)*(2 * z.^2 - x.^2 - y.^2), 'AbsTol', 1E-6);
        end
        function testl3m0(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 0, 1, 0, 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(7/pi)*(z .* (2 * z.^2 - 3 * x.^2 - 3 * y.^2)), 'AbsTol', 1E-6);
        end
        % Test by contructing a combination that yields real values that
        % can be verified using the formulae for the real coefficients
        % (since the code does not allow for imaginary values as a pdf
        % result, basis functions cannot be tested individually).
        % We use the basis functions of real spherical harmonics that are
        % given as a combination of complex spherical harmonics in
        % https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
        function testl1mneg1real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 1i * sqrt(1/2), 0, 1i * sqrt(1/2), NaN, NaN; 0, 0, 0, 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), sqrt(3/(4 * pi))*y, 'AbsTol', 1E-6);
        end
        function testl1m1real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; sqrt(1/2), 0, -sqrt(1/2), NaN, NaN; 0, 0, 0, 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), sqrt(3/(4 * pi))*x, 'AbsTol', 1E-6);
        end
        function testl2mneg2real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 1i * sqrt(1/2), 0, 0, 0, -1i * sqrt(1/2)];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/2*sqrt(15/pi)*x.*y, 'AbsTol', 1E-6);
        end
        function testl2mneg1real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, 1i * sqrt(1/2), 0, 1i * sqrt(1/2), 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/2*sqrt(15/pi)*y.*z, 'AbsTol', 1E-6);
        end
        function testl2m1real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, sqrt(1/2), 0, -sqrt(1/2), 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/2*sqrt(15/pi)*x.*z, 'AbsTol', 1E-6);
        end
        function testl2m2real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; sqrt(1/2), 0, 0, 0, sqrt(1/2)];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(15/pi)*(x.^2 - y.^2), 'AbsTol', 1E-6);
        end
        function testl3mneg3real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1i / sqrt(2), 0, 0, 0, 0, 0, 1i / sqrt(2)];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(35/(2 * pi))*y.*(3 * x.^2 - y.^2), 'AbsTol', 1E-6);
        end
        function testl3mneg2real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 1i / sqrt(2), 0, 0, 0, -1i / sqrt(2), 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/2*sqrt(105/pi)*x.*y.*z, 'AbsTol', 1E-6);
        end
        function testl3mneg1real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 1i / sqrt(2), 0, 1i / sqrt(2), 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(21/(2 * pi))*y.*(4 * z.^2 - x.^2 - y.^2), 'AbsTol', 1E-6);
        end
        function testl3m1real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 1 / sqrt(2), 0, -1 / sqrt(2), 0, 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(21/(2 * pi))*x.*(4 * z.^2 - x.^2 - y.^2), 'AbsTol', 1E-6);
        end
        function testl3m2real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 1 / sqrt(2), 0, 0, 0, 1 / sqrt(2), 0];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(105/pi)*z.*(x.^2 - y.^2), 'AbsTol', 1E-6);
        end
        function testl3m3real(testCase)
            shd = SphericalHarmonicsDistributionComplex(1/sqrt(4*pi));
            shd.coeffMat = [0, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1 / sqrt(2), 0, 0, 0, 0, 0, -1 / sqrt(2)];
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shd.pdf([x; y; z]), 1/4*sqrt(35/(2 * pi))*x.*(x.^2 - 3 * y.^2), 'AbsTol', 1E-6);
        end
        % Test conversion to RealSphericalHarmonic
        function testl0m0conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, 0, 0, 0, 0]);
        end
        function testl1mneg1conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 1i * sqrt(1/2), 0, 1i * sqrt(1/2), NaN, NaN; 0, 0, 0, 0, 0]);
        end
        function testl1m0conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 1, 0, NaN, NaN; 0, 0, 0, 0, 0]);
        end
        function testl1m1conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; sqrt(1/2), 0, -sqrt(1/2), NaN, NaN; 0, 0, 0, 0, 0]);
        end
        function testl2mneg2conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 1i * sqrt(1/2), 0, 0, 0, -1i * sqrt(1/2)]);
        end
        function testl2mneg1conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, 1i * sqrt(1/2), 0, 1i * sqrt(1/2), 0]);
        end
        function testl2m0conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, 0, 1, 0, 0]);
        end
        function testl2m1conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; 0, sqrt(1/2), 0, -sqrt(1/2), 0]);
        end
        function testl2m2conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN; sqrt(1/2), 0, 0, 0, sqrt(1/2)]);
        end
        function testl3mneg3conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1i / sqrt(2), 0, 0, 0, 0, 0, 1i / sqrt(2)]);
        end
        function testl3mneg2conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 1i / sqrt(2), 0, 0, 0, -1i / sqrt(2), 0]);
        end
        function testl3mneg1conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 1i / sqrt(2), 0, 1i / sqrt(2), 0, 0]);
        end
        function testl3m0conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 0, 1, 0, 0, 0]);
        end
        function testl3m1conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 1 / sqrt(2), 0, -1 / sqrt(2), 0, 0]);
        end
        function testl3m2conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 1 / sqrt(2), 0, 0, 0, 1 / sqrt(2), 0]);
        end
        function testl3m3conversion(testCase)
            testConversion(testCase, [1, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1 / sqrt(2), 0, 0, 0, 0, 0, -1 / sqrt(2)]);
        end
        function testTransformationViaIntegral(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            global enableExpensive
            if ~islogical(enableExpensive) || ~enableExpensive, return; end
            % Test approximating a VMF
            dist = VMFDistribution([0; -1; 0], 10);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            shd = SphericalHarmonicsDistributionComplex.fromFunctionViaIntegralCart(@(x)dist.pdf(x), 11);
            fixture.teardown;
            points = rand(3, 10000);
            testCase.verifyEqual(shd.pdf(points), dist.pdf(points), 'AbsTol', 2E-3);
            
            % Test approximating a spherical harmonic distrubution
            dist = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1 / sqrt(2), 0, 0, 0, 0, 0, -1 / sqrt(2)]);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            shd = SphericalHarmonicsDistributionComplex.fromFunctionViaIntegralCart(@(x)dist.pdf(x), 3);
            fixture.teardown;
            testCase.verifyEqual(shd.pdf(points), dist.pdf(points), 'AbsTol', 1E-6);
            testCase.verifyEqual(shd.coeffMat, dist.coeffMat, 'AbsTol', 1E-6);
        end
        function testCovergence(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            if enableExpensive
                noDiffs = 10;
            else
                noDiffs = 3;
            end
            dist = VMFDistribution([0; -1; 0], 10);
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            diffs = NaN(noDiffs, 1);
            for i = 2:noDiffs + 1
                shd = SphericalHarmonicsDistributionComplex.fromFunctionViaIntegralCart(@(x)dist.pdf(x), i);
                diffs(i-1) = shd.totalVariationDistanceNumerical(dist);
            end
            testCase.verifyLessThan(diff(diffs), 0);
        end
        function testMultiplicationFixedCoeffs(testCase)
            [phi, theta] = meshgrid(linspace(0, 2*pi, 100), linspace(-pi/2, pi/2, 100));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            
            cshd1 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN; 0, 1, 0]);
            cshd2 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1i / sqrt(2), 0, 0, 0, 0, 0, 1i / sqrt(2)]);
            testCase.verifyWarning(@()cshd1.multiply(cshd2, size(cshd1.coeffMat, 1)+size(cshd2.coeffMat, 1)-1), 'Multiplication:degreeTooHigh');
            cshdMult = cshd1.multiply(cshd2, size(cshd1.coeffMat, 1)+size(cshd2.coeffMat, 1)-2);
            testCase.verifySize(cshdMult.coeffMat, [size(cshd1.coeffMat, 1) + size(cshd2.coeffMat, 1) - 1, 2 * (size(cshd1.coeffMat, 1) + size(cshd2.coeffMat, 1) - 1) - 1]);
            % Only test if ratio is constant because cshdMult is
            % normalized and the multiplied pdf values aren't
            testCase.verifyEqual(diff(cshdMult.pdf([x; y; z])./(cshd1.pdf([x; y; z]) .* cshd2.pdf([x; y; z]))), zeros(1, size(x, 2)-1), 'AbsTol', 1E-6);
            
            cshd1 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN; sqrt(1/2), 0, -sqrt(1/2)]);
            cshd2 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1i / sqrt(2), 0, 0, 0, 0, 0, 1i / sqrt(2)]);
            cshdMult = cshd1.multiply(cshd2, size(cshd1.coeffMat, 1)+size(cshd2.coeffMat, 1)-2);
            testCase.verifyEqual(diff(cshdMult.pdf([x; y; z])./(cshd1.pdf([x; y; z]) .* cshd2.pdf([x; y; z]))), zeros(1, size(x, 2)-1), 'AbsTol', 1E-6);
            
            cshd1 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN; 1i * sqrt(1/2), 0, 1i * sqrt(1/2), NaN, NaN; 0, 0, 0, 0, 0]);
            cshd2 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN; 0, 1, 0, NaN, NaN; 0, 0, 0, 0, 0]);
            cshdMult = cshd1.multiply(cshd2, size(cshd1.coeffMat, 1)+size(cshd2.coeffMat, 1)-2);
            testCase.verifyEqual(diff(cshdMult.pdf([x; y; z])./(cshd1.pdf([x; y; z]) .* cshd2.pdf([x; y; z]))), zeros(1, size(x, 2)-1), 'AbsTol', 1E-6);
            
            cshd1 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 1i / sqrt(2), 0, 0, 0, 0, 0, 1i / sqrt(2)]);
            cshd2 = SphericalHarmonicsDistributionComplex([sqrt(1/pi) / 2, NaN, NaN, NaN, NaN, NaN, NaN; 0, 0, 0, NaN, NaN, NaN, NaN; 0, 0, 0, 0, 0, NaN, NaN; 0, 0, 0, 1, 0, 0, 0]);
            cshdMult = cshd1.multiply(cshd2, size(cshd1.coeffMat, 1)+size(cshd2.coeffMat, 1)-2);
            testCase.verifyEqual(diff(cshdMult.pdf([x; y; z])./(cshd1.pdf([x; y; z]) .* cshd2.pdf([x; y; z]))), zeros(1, size(x, 2)-1), 'AbsTol', 1E-6);
        end
        function testMultiplicationVsVmf(testCase)
            vmf = VMFDistribution([0; -1; 0], 0.5);
            shd = SphericalHarmonicsDistributionComplex.fromFunctionViaIntegralCart(@(x)vmf.pdf(x), 15);
            vmfCurr = vmf;
            shdCurr = shd;
            for i = 1:5
                vmfCurr = vmfCurr.multiply(vmf);
                shdCurr = shdCurr.multiply(shd);
            end
            [phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(shdCurr.pdf([x; y; z]), vmfCurr.pdf([x; y; z]), 'AbsTol', 1E-6);
        end
        function testTransformViaCoefficients(testCase)
            degree = 21;
            % Test identity -> sqrt
            shd = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast( ...
                VMFDistribution([1; 1; 0]/sqrt(2), 2), degree, 'identity');
            shdSqrt = shd.transformViaCoefficients('sqrt', degree);
            testCase.verifyEqual(shdSqrt.totalVariationDistanceNumerical(shd), 0, 'AbsTol', 1E-14);
            testCase.verifySize(shdSqrt.coeffMat, [degree + 1, 2 * degree + 1]);
            % Test sqrt -> identity
            shd = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast( ...
                VMFDistribution([1; 1; 0]/sqrt(2), 2), degree, 'sqrt');
            shdId = shd.transformViaCoefficients('square', degree-1);
            testCase.verifySize(shdId.coeffMat, [(degree - 1) + 1, 2 * (degree - 1) + 1]);
            testCase.verifyEqual(shdId.totalVariationDistanceNumerical(shd), 0, 'AbsTol', 1E-14);
            shdId = shd.transformViaCoefficients('square', degree+2);
            testCase.verifySize(shdId.coeffMat, [(degree + 2) + 1, 2 * (degree + 2) + 1]);
            testCase.verifyEqual(shdId.totalVariationDistanceNumerical(shd), 0, 'AbsTol', 1E-14);
        end
        function testMultiplicationViaTransform(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            wignercycle(10);
            for transformation = {'identity', 'sqrt'}
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
                shd1 = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([1; 1; 1]/sqrt(3), 2), 5, [transformation{:}]);
                shd2 = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([1; 1; -1]/sqrt(3), 2), 5, [transformation{:}]);
                fixture.teardown;
                shdMultFast = shd1.multiply(shd2);
                shdMultCoeff = shd1.multiplyViaCoefficients(shd2);
                
                shdMultFastDeg10 = shd1.multiply(shd2, 10);
                shdMultCoeffDeg10 = shd1.multiplyViaCoefficients(shd2, 10);
                
                [phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));
                [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
                testCase.verifyEqual(shdMultFast.pdf([x; y; z]), shdMultCoeff.pdf([x; y; z]), 'AbsTol', 1E-6);
                testCase.verifyEqual(shdMultFastDeg10.pdf([x; y; z]), shdMultCoeffDeg10.pdf([x; y; z]), 'AbsTol', 1E-6);
                % Verify that the above test actually is of use (i.e., there is
                % a difference when respecting degree up to 10 instead)
                testCase.verifyTrue(any(shdMultFast.pdf([x; y; z])-shdMultCoeffDeg10.pdf([x; y; z]) > 1E-6));
            end
        end
        function testConvolutionForSqrt(testCase)
            % Convolution for identity stransformation can be found in
            % AbstractSphericalHarmonicsDistributionTest. We here test that
            % convolution result approximately same as for identity
            % transformation, given a sufficiently high number of
            % coefficients.
            degree = 21;
            thisId = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([1; 1; 0]/sqrt(2), 2), degree, 'identity');
            otherId = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([0; 0; 1], 1), degree, 'identity');
            
            this = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([1; 1; 0]/sqrt(2), 2), degree, 'sqrt');
            other = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(VMFDistribution([0; 0; 1], 1), degree, 'sqrt');
            
            convolutionResultId = thisId.convolve(otherId);
            convolutionResultSqrt = this.convolve(other);
            
            testCase.verifyEqual(convolutionResultSqrt.totalVariationDistanceNumerical(convolutionResultId), 0, 'AbsTol', 1.5E-15);
        end
        function testRotationGrid(testCase)
            % Test rotation by verifying function values at identically
            % transformed points.
            rng default
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            if enableExpensive
                gridpoints = 7;
            else
                gridpoints = 3;
            end
            [phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            Rz = @(alpha)[cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];
            Ry = @(beta)[cos(beta), 0, sin(beta); 0, 1, 0; -sin(beta), 0, cos(beta)];
            
            [alphas, betas, gammas] = ndgrid(linspace(0, 2*pi, gridpoints));
            for degree = 7:10
                for i = 1:numel(alphas)
                    coeffMat = rand(7, 13);
                    coeffMat(1) = 1 / sqrt(4*pi);
                    rshd = SphericalHarmonicsDistributionReal(coeffMat);
                    shd = rshd.toSphericalHarmonicsDistributionComplex;
                    shd2 = shd.rotate(alphas(i), betas(i), gammas(i));
                    rotm = Rz(gammas(i)) * Ry(betas(i)) * Rz(alphas(i));
                    xyzRotated = rotm * [x; y; z];
                    testCase.verifyEqual(shd.pdf([x; y; z]), shd2.pdf(xyzRotated), 'AbsTol', 1E-6);
                end
            end
        end
        function testRotationSequentialRotation(testCase)
            % The angles are sequentially increased in the order they are
            % used. This is interesting not only as a unit test but also
            % for plotting and debugging purposes.
            global enableExpensive
            if ~islogical(enableExpensive) || ~enableExpensive, return; end
            [phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            Rz = @(alpha)[cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];
            Ry = @(beta)[cos(beta), 0, sin(beta); 0, 1, 0; -sin(beta), 0, cos(beta)];
            coeffMat = rand(7, 13);
            coeffMat(1) = 1 / sqrt(4*pi);
            rshd = SphericalHarmonicsDistributionReal(coeffMat);
            shd = rshd.toSphericalHarmonicsDistributionComplex;
            
            
            [alphamin, betamin, gammamin] = deal(0);
            [alphamax, betamax, gammamax] = deal(9/5*pi);
            beta = betamin;
            gamma = gammamin;
            for alpha = alphamin:pi / 20:alphamax
                shd2 = shd.rotate(alpha, beta, gamma);
                rotm = Rz(gamma) * Ry(beta) * Rz(alpha);
                xyzRotated = rotm * [x; y; z];
                testCase.verifyEqual(shd.pdf([x; y; z]), shd2.pdf(xyzRotated), 'AbsTol', 1E-6);
                if (alphamax - alpha) < pi / 40
                    for beta = betamin:pi / 40:betamax
                        shd2 = shd.rotate(alpha, beta, gamma);
                        rotm = Rz(gamma) * Ry(beta) * Rz(alpha);
                        xyzRotated = rotm * [x; y; z];
                        testCase.verifyEqual(shd.pdf([x; y; z]), shd2.pdf(xyzRotated), 'AbsTol', 1E-6);
                        if (betamax - beta) < pi / 40
                            for gamma = gammamin:pi / 40:gammamax
                                shd2 = shd.rotate(alpha, beta, gamma);
                                rotm = Rz(gamma) * Ry(beta) * Rz(alpha);
                                xyzRotated = rotm * [x; y; z];
                                testCase.verifyEqual(shd.pdf([x; y; z]), shd2.pdf(xyzRotated), 'AbsTol', 1E-6);
                            end
                        end
                    end
                end
            end
        end
        function testMean(testCase)
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; 0, 1, 0]);
            testCase.verifyEqual(shd.meanDirection, [0; 0; 1], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; 0, -1, 0]);
            testCase.verifyEqual(shd.meanDirection, [0; 0; -1], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; 1i * sqrt(1/2), 0, 1i * sqrt(1/2)]);
            testCase.verifyEqual(shd.meanDirection, [0; 1; 0], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; -1i * sqrt(1/2), 0, -1i * sqrt(1/2)]);
            testCase.verifyEqual(shd.meanDirection, [0; -1; 0], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; sqrt(1/2), 0, -sqrt(1/2)]);
            testCase.verifyEqual(shd.meanDirection, [1; 0; 0], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; -sqrt(1/2), 0, sqrt(1/2)]);
            testCase.verifyEqual(shd.meanDirection, [-1; 0; 0], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; -1i * sqrt(1/2) - sqrt(1/2), 0, -1i * sqrt(1/2) + sqrt(1/2)]);
            testCase.verifyEqual(shd.meanDirection, 1/sqrt(2)*[-1; -1; 0], 'AbsTol', 1E-6);
            shd = SphericalHarmonicsDistributionComplex([1 / sqrt(4*pi), NaN, NaN; 1i * sqrt(1/2) + sqrt(1/2), 1, 1i * sqrt(1/2) - sqrt(1/2)]);
            testCase.verifyEqual(shd.meanDirection, 1/sqrt(3)*[1; 1; 1], 'AbsTol', 1E-6);
            testCase.verifyTrue(isreal(shd.meanDirection));
        end
        function integralAnalytical(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            coeffRand = rand(1, 9);
            unnormalizedCoeffs = [coeffRand(1), NaN, NaN, NaN, NaN; coeffRand(2) + 1i * coeffRand(3), coeffRand(4), -coeffRand(2) + 1i * coeffRand(3), NaN, NaN; ...
                coeffRand(5) + 1i * coeffRand(6), coeffRand(7) + 1i * coeffRand(8), coeffRand(9), -coeffRand(7) + 1i * coeffRand(8), coeffRand(5) - 1i * coeffRand(6)];
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            shd = SphericalHarmonicsDistributionComplex(unnormalizedCoeffs);
            fixture.teardown;
            testCase.verifyEqual(shd.integralAnalytical, shd.integral, 'AbsTol', 1E-6);
        end
        function testfromDistributionNumericalFast(testCase)
            [phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            dists = {VMFDistribution([1; 0; 0], 1), VMFDistribution([0; 1; 0], 2), VMFDistribution([0; 0; 1], 3), ...
                HypersphericalMixture({VMFDistribution([1; 1; 1]/sqrt(3), 3), VMFDistribution([1; 0; 0], 2)}, [0.5, 0.5])};
            for transformation = {'identity', 'sqrt'}
                for i = 1:numel(dists)
                    shd = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(dists{i}, 21, [transformation{:}]);
                    testCase.verifyEqual(dists{i}.pdf([x; y; z]), shd.pdf([x; y; z]), 'AbsTol', 1E-6)
                end
            end
        end
        function testFromGrid(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            dist = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Transformation:notEq_Point_set'));
            sgd = SphericalGridDistribution.fromDistribution(dist, 1012, 'sh_grid');
            fixture.teardown();
            
            % Test without providing grid
            deg = (-6+sqrt(36-8*(4-numel(sgd.gridValues))))/4;
            shd1Id = SphericalHarmonicsDistributionComplex.fromGrid(reshape(sgd.gridValues,deg+2,2*deg+2));
            shd1Sqrt = SphericalHarmonicsDistributionComplex.fromGrid(reshape(sgd.gridValues,deg+2,2*deg+2), [], 'sqrt');
            % Test when providing grid
            shd2Id = SphericalHarmonicsDistributionComplex.fromGrid(sgd.gridValues, sgd.getGrid());
            shd2Sqrt = SphericalHarmonicsDistributionComplex.fromGrid(sgd.gridValues, sgd.getGrid(), 'sqrt');
            
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            
            testCase.verifyEqual(shd1Id.pdf([x; y; z]), dist.pdf([x; y; z]), 'AbsTol', 1E-11);
            testCase.verifyEqual(shd1Sqrt.pdf([x; y; z]), dist.pdf([x; y; z]), 'AbsTol', 1E-11);
            testCase.verifyEqual(shd2Id.pdf([x; y; z]), dist.pdf([x; y; z]), 'AbsTol', 1E-11);
            testCase.verifyEqual(shd2Sqrt.pdf([x; y; z]), dist.pdf([x; y; z]), 'AbsTol', 1E-11);
        end
    end
    methods
        function testConversion(testCase, coeffMat)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            cshd = SphericalHarmonicsDistributionComplex(coeffMat);
            fixture.teardown;
            rshd = cshd.toSphericalHarmonicsDistributionReal;
            [phi, theta] = meshgrid(linspace(0, 2*pi, 10), linspace(-pi/2, pi/2, 10));
            [x, y, z] = sph2cart(phi(:)', theta(:)', 1);
            testCase.verifyEqual(rshd.pdf([x; y; z]), cshd.pdf([x; y; z]), 'AbsTol', 1E-6);
        end
    end
end

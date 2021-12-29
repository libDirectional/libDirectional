classdef ToroidalFourierDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testAgainstHypertoroidal(testCase)
            tvm = ToroidalVMSineDistribution([1; 4], [0.3; 0.7], 0.5);
            tfd = ToroidalFourierDistribution.fromDistribution(tvm, [15, 15], 'sqrt');
            hfd = HypertoroidalFourierDistribution.fromDistribution(tvm, [15, 15], 'sqrt');
            testPoints = rand(2, 100) * 4 * pi - pi;
            testCase.verifyEqual(tfd.pdf(testPoints), hfd.pdf(testPoints), 'AbsTol', 1E-6);
        end
        function testSingleElementVector(testCase)
            testCase.verifyWarning(@(x)ToroidalFourierDistribution([0; 1; 0]), 'ToroidalFourierDistribution:VectorGiven');
        end
        function testFromFunction(testCase)
            tvm = ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5);
            tfd1Fun = ToroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), [15, 15], 'identity');
            tfd2Fun = ToroidalFourierDistribution.fromFunction(@(x, y)reshape(tvm.pdf([x(:)'; y(:)']), size(x)), [15, 15], 'sqrt');
            for i = -2:3
                tvmMoment = tvm.trigonometricMoment(i);
                testCase.verifyEqual(tfd1Fun.trigonometricMoment(i), tvmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(tfd1Fun.trigonometricMomentNumerical(i), tvmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(tfd2Fun.trigonometricMoment(i), tvmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(tfd2Fun.trigonometricMomentNumerical(i), tvmMoment, 'AbsTol', 1E-6);
            end
        end
        function testFromFunctionAndMoments(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            tvm = ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5);
            tfd1 = ToroidalFourierDistribution.fromDistribution(tvm, [15, 15], 'identity');
            tfd2 = ToroidalFourierDistribution.fromDistribution(tvm, [15, 15], 'sqrt');
            for i = -2:3
                tvmMoment = tvm.trigonometricMoment(i);
                testCase.verifyEqual(tfd1.trigonometricMoment(i), tvmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(tfd1.trigonometricMomentNumerical(i), tvmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(tfd2.trigonometricMoment(i), tvmMoment, 'AbsTol', 1E-6);
                testCase.verifyEqual(tfd2.trigonometricMomentNumerical(i), tvmMoment, 'AbsTol', 1E-6);
            end
            testCase.applyFixture(SuppressedWarningsFixture('Truncate:TooFewCoefficients'));
            testCase.verifyEqual(tfd1.trigonometricMoment(16), [0; 0]);
            testCase.verifyEqual(tfd2.trigonometricMoment(30), [0; 0]);
        end
        function testCorrelation(testCase)
            tvm = ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5);
            tfd1 = ToroidalFourierDistribution.fromDistribution(tvm, [101, 101], 'identity');
            tfd2 = ToroidalFourierDistribution.fromDistribution(tvm, [101, 101], 'sqrt');
            testCase.verifyEqual(tfd1.circularCorrelationJammalamadaka, tvm.circularCorrelationJammalamadakaNumerical, 'AbsTol', 1E-6);
            testCase.verifyEqual(tfd2.circularCorrelationJammalamadaka, tvm.circularCorrelationJammalamadakaNumerical, 'AbsTol', 1E-6);
            testCase.verifyEqual(tfd1.covariance4D, tvm.covariance4D, 'AbsTol', 1E-6);
            testCase.verifyEqual(tfd2.covariance4D, tvm.covariance4D, 'AbsTol', 1E-6);
        end
        function testIntegral(testCase)
            tvm = ToroidalVMSineDistribution([1; 2], [0.6; 1], 0.2);
            tfdid = ToroidalFourierDistribution.fromDistribution(tvm, [45, 49], 'identity');
            tfdsqrt = ToroidalFourierDistribution.fromDistribution(tvm, [45, 49], 'sqrt');
            l = [0.3; 0.3];
            r = [1.5; 1.5];
            %l1=0.3;r1=1.5;
            %l2=0.3;r2=1.5;
            testCase.verifyEqual(tfdid.integral(l, r), tvm.integral(l, r), 'AbsTol', 1E-6);
            testCase.verifyEqual(tfdid.integralNumerical(l, r), tvm.integral(l, r), 'AbsTol', 1E-6);
            testCase.verifyEqual(tfdsqrt.integral(l, r), tvm.integral(l, r), 'AbsTol', 1E-6);
            testCase.verifyEqual(tfdsqrt.integralNumerical(l, r), tvm.integral(l, r), 'AbsTol', 1E-6);
            % Test special case case
            tfdsimple = ToroidalFourierDistribution(blkdiag(0, 1/(4 * pi^2), 0), 'identity');
            testCase.verifyEqual(tfdsimple.integral(), 1, 'AbsTol', 1E-6);
        end
        function testToTWN(testCase)
            mu = [1; 3];
            C = [9, 0.3; 0.3, 2];
            twn = ToroidalWNDistribution(mu, C);
            tfd = ToroidalFourierDistribution.fromDistribution(twn, 9, 'identity');
            twnConv = tfd.toTWN;
            testCase.verifyEqual(twnConv.mu, twn.mu, 'AbsTol', 1E-4);
            testCase.verifyEqual(twnConv.C, twn.C, 'AbsTol', 1E-4);
            tfd = ToroidalFourierDistribution.fromDistribution(twn, 9, 'sqrt');
            twnConv = tfd.toTWN;
            testCase.verifyEqual(twnConv.mu, twn.mu, 'AbsTol', 1E-4);
            testCase.verifyEqual(twnConv.C, twn.C, 'AbsTol', 1E-4);
        end
    end
end

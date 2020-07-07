classdef FourierFilterTest < matlab.unittest.TestCase
    
    methods(Test)
        function testInitialization(testCase)
            filter = FourierFilter(101);
            testCase.verifyClass(filter.fd, 'FourierDistribution');
            testCase.verifyLength(filter.fd.a, 51);
            testCase.verifyLength(filter.fd.b, 50);
            testCase.verifyEqual(filter.fd.transformation, 'sqrt');
            testCase.verifyEqual(filter.fd.a(1), sqrt(2/pi), 'AbsTol', 1E-10);
        end
        
        function testSetState(testCase)
            filter = FourierFilter(101);
            fd = FourierDistribution.fromDistribution(VMDistribution(2, 3), 101);
            filter.setState(fd);
            testCase.verifyEqual(filter.fd, fd);
            testCase.verifyWarning(@()filter.setState(VMDistribution(2, 3)), 'setState:nonFourier');
            fdtmp = fd;
            fdtmp.a = fd.a(1:end-1);
            fdtmp.b = fd.b(1:end-1);
            testCase.verifyWarning(@()filter.setState(fdtmp), 'setState:noOfCoeffsDiffer');
            fdtmp = fd;
            fdtmp.transformation = 'identity';
            testCase.verifyWarning(@()filter.setState(fdtmp), 'setState:transDiffer');
            testCase.verifyWarning(@()filter.setState(VMDistribution(2, 3)), 'setState:nonFourier');
        end
        
        function testPrediction(testCase)
            % Validity of operation is tested in Fourier class, only test
            % filter interface
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            for transformation = {'identity', 'sqrt'}
                filter = FourierFilter(101, [transformation{:}]);
                vm = VMDistribution(3, 2);
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Predict:automaticConversion'));
                filter.predictIdentity(vm);
                fixture.teardown;
                fd1 = FourierDistribution.fromDistribution(CircularUniformDistribution(), 101, [transformation{:}]);
                fd2 = fd1.convolve(FourierDistribution.fromDistribution(vm, 101, [transformation{:}]));
                testCase.verifyLength(filter.fd.a, 51);
                testCase.verifyLength(filter.fd.b, 50);
                testCase.verifyEqual(filter.fd.a, fd2.a, 'AbsTol', 1E-8);
                testCase.verifyEqual(filter.fd.b, fd2.b, 'AbsTol', 1E-8);
            end
        end
        
        function testUpdateIdentity(testCase)
            % As above, only test filter interface
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            for transformation = {'identity', 'sqrt'}
                filter = FourierFilter(101, [transformation{:}]);
                vmMultiply = VMDistribution(3, 2);
                vmFilter = VMDistribution(0, 2);
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Predict:automaticConversion'));
                filter.updateIdentity(vmFilter, 3);
                fixture.teardown;
                fd1 = FourierDistribution.fromDistribution(CircularUniformDistribution(), 101, [transformation{:}]);
                fd2 = fd1.multiply(FourierDistribution.fromDistribution(vmMultiply, 101, [transformation{:}]));
                fd3 = fd2.truncate(101);
                testCase.verifyLength(filter.fd.a, 51);
                testCase.verifyLength(filter.fd.b, 50);
                testCase.verifyEqual(filter.fd.a, fd3.a, 'AbsTol', 1E-8);
                testCase.verifyEqual(filter.fd.b, fd3.b, 'AbsTol', 1E-8);
            end
            testCase.verifyWarning(@()filter.updateIdentity(WNDistribution(1, 2), 3), 'Update:automaticConversion');
        end
        
        function testUpdateNonlinear(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            z = 3;
            likelihood = @(z, x)1 ./ (abs(x-z) + .1);
            for transformation = {'identity', 'sqrt'}
                filter = FourierFilter(23, [transformation{:}]);
                hfd1 = FourierDistribution.fromDistribution(VMDistribution(3, 2), 23, [transformation{:}]);
                filter.setState(hfd1);
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
                filter.updateNonlinear(likelihood, z);
                fixture.teardown;
                testCase.verifyEqual(filter.getEstimateMean, z, 'AbsTol', 0.2); % The mean should stay approximately equal, but not entirely
            end
        end
        
        function testUpdateNonlinearViaIFFT(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            z = 3;
            likelihood = @(z, x)1 ./ (abs(x-z) + .1);
            noCoeffs = 501;
            for transformation = {'identity', 'sqrt'}
                filter1 = FourierFilter(noCoeffs, [transformation{:}]);
                filter2 = FourierFilter(noCoeffs, [transformation{:}]);
                hfd1 = FourierDistribution.fromDistribution(VMDistribution(3, 2), noCoeffs, [transformation{:}]);
                filter1.setState(hfd1);
                filter2.setState(hfd1);
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
                filter1.updateNonlinear(likelihood, z);
                filter2.updateNonlinearViaIFFT(likelihood, z);
                fixture.teardown;
                % Given a high number of coefficients, the results should
                % be similar
                testCase.verifyEqual(filter1.fd.totalVariationDistanceNumerical(filter2.fd), ...
                    0, 'AbsTol', 1E-3);
            end
        end
        
        function testPredictNonlinearForLinear(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            densityInit = VMDistribution(3, 5);
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 31;
            for transformation = {'identity', 'sqrt'}
                fourierFilterLin = FourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterLin.setState(FourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                fourierFilterNl = FourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterNl.setState(FourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                
                testCase.applyFixture(SuppressedWarningsFixture('Predict:automaticConversion'));
                fourierFilterLin.predictIdentity(fNoiseDist)
                fourierFilterNl.predictNonlinear(@(x)x, fNoiseDist, true)
                testCase.verifyEqual(fourierFilterLin.getEstimate.kldNumerical(fourierFilterNl.getEstimate), 0, 'AbsTol', 1E-8);
                
                fNoiseDistShifted = fNoiseDist.shift(1);
                fourierFilterLin.predictIdentity(fNoiseDistShifted)
                fourierFilterNl.predictNonlinear(@(x)x+1, fNoiseDist, false)
                testCase.verifyEqual(fourierFilterLin.getEstimate.kldNumerical(fourierFilterNl.getEstimate), 0, 'AbsTol', 1E-8);
            end
        end
        
        function testPredictIdentityGrid(testCase)
            noCoeffs = 101;
            fd = FourierDistribution.fromDistribution(VMDistribution(1, 1), noCoeffs, 'sqrt');
            sysNoise = VMDistribution(2, 2);
            fdPred = fd.convolve(FourierDistribution.fromDistribution(sysNoise, noCoeffs, 'sqrt'));
            
            ff = FourierFilter(noCoeffs, 'sqrt');
            ff.setState(fd);
            ff.predictIdentity(sysNoise.pdf((0:2 * pi / noCoeffs:2 * pi - 2 * pi / noCoeffs)));
            
            testCase.verifyEqual(fdPred.totalVariationDistanceNumerical(ff.getEstimate), 0, 'AbsTol', 1E-13);
        end
        
        function testTransitionDensity(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test that transition density is correct if .C is multiply
            % with 2*pi for identity and with sqrt(2*pi) for sqrt
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 401;
            % Igor's function
            aGen = @(a)@(te) pi * (sin(sign(te-pi)/2.*abs(te-pi).^a/pi^(a - 1)) + 1);
            a = aGen(4);
            
            fTrans = @(xkk, xk)reshape(fNoiseDist.pdf(xkk(:)'-a(xk(:)')), size(xk));
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdTransId = ToroidalFourierDistribution.fromFunction(fTrans, noCoeffs*[1, 1], 'identity');
            hfdTransId.C = 2 * pi * hfdTransId.C;
            hfdTransSqrt = ToroidalFourierDistribution.fromFunction(fTrans, noCoeffs*[1, 1], 'sqrt');
            hfdTransSqrt.C = sqrt(2*pi) * hfdTransSqrt.C;
            fixture.teardown;
            
            [xVals, yVals] = meshgrid(linspace(0, 2*pi, 50));
            fTransVals = reshape(fTrans(xVals(:)', yVals(:)'), size(xVals));
            hfdTransIdVals = reshape(hfdTransId.pdf([xVals(:)'; yVals(:)']), size(xVals));
            hfdTransSqrtVals = reshape(hfdTransSqrt.pdf([xVals(:)'; yVals(:)']), size(xVals));
            
            testCase.verifyEqual(fTransVals, hfdTransIdVals, 'AbsTol', 3E-6);
            testCase.verifyEqual(fTransVals, hfdTransSqrtVals, 'AbsTol', 3E-6);
        end
        
        function testJointDensity(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test that transition density is correct if .C is multiply
            % with 2*pi for identity and with sqrt(2*pi) for sqrt
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 401;
            % Igor's function
            aGen = @(a)@(te) pi * (sin(sign(te-pi)/2.*abs(te-pi).^a/pi^(a - 1)) + 1);
            a = aGen(4);
            
            densityInit = VMDistribution(3, 5);
            fdInitId = FourierDistribution.fromDistribution(densityInit, noCoeffs, 'identity');
            fdInitSqrt = FourierDistribution.fromDistribution(densityInit, noCoeffs, 'sqrt');
            
            fTrans = @(xkk, xk)reshape(fNoiseDist.pdf(xkk(:)'-a(xk(:)')), size(xk));
            % Not correctly scaled here.
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdTransId = ToroidalFourierDistribution.fromFunction(fTrans, noCoeffs*[1, 1], 'identity');
            hfdTransSqrt = ToroidalFourierDistribution.fromFunction(fTrans, noCoeffs*[1, 1], 'sqrt');
            fixture.teardown;
            
            fJoint = @(xkk, xk)fTrans(xkk, xk) .* reshape(densityInit.pdf(xk(:)'), size(xk));
            % Initially first and overwrite C (this ensures no additional
            % normalization takes place)
            hfdJointId = ToroidalFourierDistribution(diag([0, 1 / (2 * pi)^2, 0]), 'identity');
            hfdJointSqrt = ToroidalFourierDistribution(diag([0, sqrt(1/(2 * pi)^2), 0]), 'sqrt');
            % This intermediate step is not performed in FourierFilter for
            % the identity transformation. Next test case shows it is equal
            % for identity transformation. Scale up fTrans after calculation of
            % joint.
            hfdJointId.C = (2 * pi) * conv2(hfdTransId.C, fdInitId.c, 'same');
            % Keep calculations as close to the code as possible. C of
            % hfdTrans is scaled up before convolution.
            hfdJointSqrt.C = conv2(sqrt(2*pi)*hfdTransSqrt.C, fdInitSqrt.c, 'same');
            
            [xVals, yVals] = meshgrid(linspace(0, 2*pi, 50));
            fJointVals = reshape(fJoint(xVals(:)', yVals(:)'), size(xVals));
            hfdJointIdVals = reshape(hfdJointId.pdf([xVals(:)'; yVals(:)']), size(xVals));
            hfdJointSqrtVals = reshape(hfdJointSqrt.pdf([xVals(:)'; yVals(:)']), size(xVals));
            
            testCase.verifyEqual(fJointVals, hfdJointIdVals, 'AbsTol', 2E-9);
            testCase.verifyEqual(fJointVals, hfdJointSqrtVals, 'AbsTol', 1E-9);
        end
        
        function testTruncatePredictionJointId(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test that using truncation with .valid when using the
            % identity transformation is valid
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 401;
            % Igor's function
            aGen = @(a)@(te) pi * (sin(sign(te-pi)/2.*abs(te-pi).^a/pi^(a - 1)) + 1);
            a = aGen(4);
            
            densityInit = VMDistribution(3, 5);
            fdInitId = FourierDistribution.fromDistribution(densityInit, noCoeffs, 'identity');
            
            fTrans = @(xkk, xk)reshape(fNoiseDist.pdf(xkk(:)'-a(xk(:)')), size(xk));
            % Not correctly scaled here.
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdTransId = ToroidalFourierDistribution.fromFunction(fTrans, noCoeffs*[1, 1], 'identity');
            fixture.teardown;
            
            % Initially first and overwrite C (this ensures no additional
            % normalization takes place)
            hfdJointIdSame = ToroidalFourierDistribution(diag([0, 1 / (2 * pi)^2, 0]), 'identity');
            hfdJointIdFull = ToroidalFourierDistribution(diag([0, 1 / (2 * pi)^2, 0]), 'identity');
            % Obtain joint (truncated to 'same')
            hfdJointIdSame.C = conv2((2 * pi)*hfdTransId.C, fdInitId.c, 'same');
            hfdJointIdFull.C = conv2((2 * pi)*hfdTransId.C, fdInitId.c, 'full');
            % Marginalization: take central row and multiply with 2*pi
            cPredictedIdUsingSame = (2 * pi) * hfdJointIdSame.C(:, (size(hfdJointIdSame.C, 2) + 1)/2);
            cPredictedIdUsingFull = (2 * pi) * hfdJointIdFull.C(:, (size(hfdJointIdFull.C, 2) + 1)/2);
            
            % Using FourierFilter class
            fourierFilterNl = FourierFilter(noCoeffs, 'identity');
            fourierFilterNl.setState(FourierDistribution.fromDistribution(densityInit, noCoeffs, 'identity'));
            fourierFilterNl.predictNonlinear(a, fNoiseDist);
            
            cPredictedIdUsingValid = fourierFilterNl.fd.c;
            
            testCase.verifyEqual(cPredictedIdUsingSame.', cPredictedIdUsingValid, 'AbsTol', 2E-16);
            testCase.verifyEqual(cPredictedIdUsingFull.', cPredictedIdUsingValid, 'AbsTol', 2E-16);
        end
        
        function testPredictNonlinear(testCase)
            densityInit = VMDistribution(3, 5);
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 101;
            % Igor's function
            aGen = @(a)@(te) pi * (sin(sign(te-pi)/2.*abs(te-pi).^a/pi^(a - 1)) + 1);
            a = aGen(4);
            
            fTrans = @(xkk, xk)reshape(fNoiseDist.pdf(xkk(:)'-a(xk(:)')), size(xk));
            fPredUsingInt = @(xkk)arrayfun(@(xkkCurr) ...
                integral(@(xkCurr)fTrans(xkkCurr, xkCurr).*reshape(densityInit.pdf(xkCurr(:)'), size(xkCurr)), 0, 2*pi), xkk);
            xVals = linspace(0, 2*pi, 100);
            
            for transformation = {'identity', 'sqrt'}
                fourierFilterNl = FourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterNl.setState(FourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                fourierFilterNl.predictNonlinear(a, fNoiseDist);
                testCase.verifyEqual(fourierFilterNl.fd.pdf(xVals), fPredUsingInt(xVals), 'AbsTol', 1E-9);
            end
            
        end
        
        function showUpdateViaIFFTCouldBeIdentical(testCase)
            % Show that if we padded before the IFFT, we could obtain the
            % same multiplication result.
            noCoeffs = 3;
            prior = FourierDistribution.fromDistribution(VMDistribution(1, 1), noCoeffs, 'identity');
            likelihood = FourierDistribution.fromDistribution(VMDistribution(2, 2), noCoeffs, 'identity');
            
            cUnnormViaConv = conv(prior.c, likelihood.c, 'full');
            
            xVals = linspace(0, 2*pi, 6);
            xVals(end) = [];
            cUnnormViaFFT = fftshift(fft(ifft(5*ifftshift([0, prior.c, 0])).*likelihood.pdf(xVals)*1/5));
            
            testCase.verifyEqual(cUnnormViaFFT, cUnnormViaConv, 'AbsTol', 1E-17);
        end
        
        function testGetEstimateMean(testCase)
            vm = VMDistribution(3, 0.1);
            for transformation = {'identity', 'sqrt'}
                filter = FourierFilter(31, [transformation{:}]);
                filter.setState(FourierDistribution.fromDistribution(vm, 2*length(filter.fd.a)-1, [transformation{:}]));
                testCase.verifyEqual(filter.getEstimateMean, vm.mu, 'AbsTol', 1E-9);
            end
        end
        
        function testAssociationLikelihood(testCase)
            vm1 = VMDistribution(3, 3);
            vm2 = VMDistribution(1, 10);
            
            ffid = FourierFilter(101, 'identity');
            ffsqrt = FourierFilter(101, 'sqrt');
            ffid.setState(FourierDistribution.fromDistribution(vm1, 101, 'identity'));
            ffsqrt.setState(FourierDistribution.fromDistribution(vm1, 101, 'sqrt'));
            
            fdLikelihoodId = FourierDistribution.fromDistribution(vm2, 101, 'identity');
            fdLikelihoodSqrt = FourierDistribution.fromDistribution(vm2, 101, 'sqrt');
            
            testCase.verifyEqual(ffid.associationLikelihood(fdLikelihoodId), ...
                ffid.associationLikelihoodNumerical(fdLikelihoodId), 'AbsTol', 1E-10);
            testCase.verifyEqual(ffsqrt.associationLikelihood(fdLikelihoodSqrt), ...
                ffsqrt.associationLikelihoodNumerical(fdLikelihoodSqrt), 'AbsTol', 1E-10);
        end
    end
end

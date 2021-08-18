classdef HypertoroidalFourierFilterTest < matlab.unittest.TestCase
    methods(Test)
        function testInitialization(testCase)
            filter = HypertoroidalFourierFilter([23, 23]);
            testCase.verifyClass(filter.hfd, 'HypertoroidalFourierDistribution');
            testCase.verifySize(filter.hfd.C, [23, 23]);
            testCase.verifyEqual(filter.hfd.transformation, 'sqrt');
            testCase.verifyEqual(filter.hfd.C(12, 12), 1/sqrt(4*pi^2), 'AbsTol', 1E-10);
            filter = HypertoroidalFourierFilter([23, 25, 27]);
            testCase.verifyClass(filter.hfd, 'HypertoroidalFourierDistribution');
            testCase.verifySize(filter.hfd.C, [23, 25, 27]);
            testCase.verifyEqual(filter.hfd.C(12, 13, 14), 1/sqrt(8*pi^3), 'AbsTol', 1E-10);
            filter = HypertoroidalFourierFilter([23, 25, 27], 'identity');
            testCase.verifyEqual(filter.hfd.transformation, 'identity');
        end
        
        function testSetState(testCase)
            % Test warnings
            filter = HypertoroidalFourierFilter([21, 21, 5]);
            testCase.verifyWarning(@()filter.setState(HypertoroidalWNDistribution([1; 1; 1], eye(3))), 'setState:nonFourier');
            hfdtmp = filter.hfd;
            hfdtmp.C = hfdtmp.C(:, :, 3);
            testCase.verifyWarning(@()filter.setState(hfdtmp), 'setState:noOfDimsDiffer');
            hfdtmp.C = hfdtmp.C(2:20, 2:20);
            testCase.verifyWarning(@()filter.setState(hfdtmp), 'setState:noOfCoeffsDiffer');
            hfdtmp.transformation = 'identity';
            testCase.verifyWarning(@()filter.setState(hfdtmp), 'setState:transDiffer');
            % Test correct functionality
            filter = HypertoroidalFourierFilter([23, 25]);
            filter.setState(HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5), [23, 25]));
            testCase.verifyEqual(filter.hfd, HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([1; 2], [0.3; 0.5], 0.5), [23, 25]));
        end
        
        function testPredictIdentity(testCase)
            % Only test filter interface, validity of operation is tested
            % in HypertoroidalFourierDistribution
            for transformation = {'identity', 'sqrt'}
                filter = HypertoroidalFourierFilter([23, 25], [transformation{:}]);
                hfd1 = HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5), [23, 25], [transformation{:}]);
                hfd2 = HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([0; 0], [0.8; 1.5], 0.2), [23, 25], [transformation{:}]);
                hfdRes = hfd1.convolve(hfd2, size(hfd1.C));
                filter.setState(hfd1);
                filter.predictIdentity(hfd2);
                testCase.verifyClass(filter.hfd, 'HypertoroidalFourierDistribution');
                testCase.verifySize(filter.hfd.C, [23, 25]);
                testCase.verifyEqual(filter.hfd.C, hfdRes.C, 'Rel', 1E-5);
            end
            testCase.verifyWarning(@()filter.predictIdentity(HypertoroidalWNDistribution([1; 1], eye(2))), 'PredictIdentity:automaticConversion');
        end
        
        function testUpdateIdentity(testCase)
            % As above, only test filter interface
            z = [3; 4];
            for transformation = {'identity', 'sqrt'}
                filter = HypertoroidalFourierFilter([23, 25], [transformation{:}]);
                hfd1 = HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5), [23, 25], [transformation{:}]);
                hfd2 = HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution(z, [0.8; 1.5], 0.2), [23, 25], [transformation{:}]);
                hfd2zeroMean = HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([0; 0], [0.8; 1.5], 0.2), [23, 25], [transformation{:}]);
                hfdRes = hfd1.multiply(hfd2);
                hfdResTrunc = hfdRes.truncate([23, 25]);
                filter.setState(hfd1);
                filter.updateIdentity(hfd2zeroMean, z);
                testCase.verifyClass(filter.hfd, 'HypertoroidalFourierDistribution');
                testCase.verifySize(filter.hfd.C, [23, 25]);
                testCase.verifyEqual(filter.hfd.C, hfdResTrunc.C, 'AbsTol', 1E-6);
            end
            testCase.verifyWarning(@()filter.updateIdentity(HypertoroidalWNDistribution([1; 1], eye(2)), z), 'Update:automaticConversion');
        end
        
        function testUpdateNonlinear2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            z = [1; 3];
            for transformation = {'identity', 'sqrt'}
                filter = HypertoroidalFourierFilter([23, 25], [transformation{:}]);
                hfd1 = HypertoroidalFourierDistribution.fromDistribution(ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5), [23, 25], [transformation{:}]);
                likelihood = @(z, x)1 ./ (sum(abs(x-z)) + .1);
                filter.setState(hfd1);
                fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
                filter.updateNonlinear(likelihood, z);
                fixture.teardown;
                testCase.verifyEqual(filter.getEstimateMean, z, 'AbsTol', 0.2); % The mean should stay approximately equal, but not entirely
            end
        end
        
        function testPredictNonlinearForLinear1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture 
            densityInit = VMDistribution(3, 5);
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 31;
            for transformation = {'identity', 'sqrt'}
                fourierFilterLin = HypertoroidalFourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterLin.setState(HypertoroidalFourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                fourierFilterNl = HypertoroidalFourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterNl.setState(HypertoroidalFourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                
                testCase.applyFixture(SuppressedWarningsFixture('PredictIdentity:automaticConversion'));
                fourierFilterLin.predictIdentity(fNoiseDist)
                fourierFilterNl.predictNonlinear(@(x)x, fNoiseDist, true)
                testCase.verifyEqual(fourierFilterLin.getEstimate.kldNumerical(fourierFilterNl.getEstimate), 0, 'AbsTol', 1E-8);
                
                fNoiseDistShifted = fNoiseDist.shift(1);
                fourierFilterLin.predictIdentity(fNoiseDistShifted)
                fourierFilterNl.predictNonlinear(@(x)x+1, fNoiseDist, false)
                testCase.verifyEqual(fourierFilterLin.getEstimate.kldNumerical(fourierFilterNl.getEstimate), 0, 'AbsTol', 1E-8);
            end
        end
        
        function testPredictNonlinearForLinear2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            densityInit = HypertoroidalWNDistribution([3; 3], eye(2));
            fNoiseDist = HypertoroidalWNDistribution([0; 0], [1, 0.1; 0.1, 1]);
            noCoeffs = [15, 15];
            for transformation = {'identity', 'sqrt'}
                if strcmp([transformation{:}], 'identity')
                    tol = 1E-8;
                else
                    tol = 1E-6;
                end
                fourierFilterLin = HypertoroidalFourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterLin.setState(HypertoroidalFourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                fourierFilterNl = HypertoroidalFourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterNl.setState(HypertoroidalFourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                
                testCase.applyFixture(SuppressedWarningsFixture('PredictIdentity:automaticConversion'));
                fourierFilterLin.predictIdentity(fNoiseDist);
                fourierFilterNl.predictNonlinear(@(varargin)varargin{:}, fNoiseDist);
                testCase.verifyEqual(fourierFilterNl.hfd.transformation, [transformation{:}]);
                testCase.verifyEqual(fourierFilterLin.getEstimate.C, fourierFilterNl.getEstimate.C, 'AbsTol', tol);
                
                fNoiseDistShifted = fNoiseDist.shift([1; 1]);
                fourierFilterLin.predictIdentity(fNoiseDistShifted);
                fourierFilterNl.predictNonlinear(@(x, y)deal(x+1, y+1), fNoiseDist);
                testCase.verifyEqual(fourierFilterNl.hfd.transformation, [transformation{:}]);
                testCase.verifyEqual(fourierFilterLin.getEstimate.C, fourierFilterNl.getEstimate.C, 'AbsTol', tol);
            end
        end
        
        function testPredictNonlinear1D(testCase)
            densityInit = VMDistribution(3, 5);
            fNoiseDist = VMDistribution(0.5, 10);
            noCoeffs = 31;
            % Igor's function
            aGen = @(a)@(te) pi * (sin(sign(te-pi)/2.*abs(te-pi).^a/pi^(a - 1)) + 1);
            a = aGen(4);
            
            fTrans = @(xkk, xk)reshape(fNoiseDist.pdf(xkk(:)'-a(xk(:)')), size(xk));
            fPredUsingInt = @(xkk)arrayfun(@(xkkCurr) ...
                integral(@(xkCurr)fTrans(xkkCurr, xkCurr).*reshape(densityInit.pdf(xkCurr(:)'), size(xkCurr)), 0, 2*pi), xkk);
            xvals = linspace(0, 2*pi, 100);
            for transformation = {'identity', 'sqrt'}
                fourierFilterNl = HypertoroidalFourierFilter(noCoeffs, [transformation{:}]);
                fourierFilterNl.setState(HypertoroidalFourierDistribution.fromDistribution(densityInit, noCoeffs, [transformation{:}]));
                fourierFilterNl.predictNonlinear(a, fNoiseDist);
                testCase.verifyEqual(fourierFilterNl.hfd.pdf(xvals), fPredUsingInt(xvals), 'AbsTol', 2E-5);
            end
        end
        
        function testWarningFreePredictNonlinearViaTransitionDensity(testCase)
            fourierFilterNl = HypertoroidalFourierFilter([11, 11]);
            testCase.verifyError(@()fourierFilterNl.predictNonlinearViaTransitionDensity(HypertoroidalWNDistribution([1; 1; 1], eye(3))), 'predictNonlinear:fTransInvalid');
            testCase.verifyWarningFree(@()fourierFilterNl.predictNonlinearViaTransitionDensity(HypertoroidalFourierDistribution.fromDistribution( ...
                HypertoroidalWNDistribution([1; 1; 1; 1], eye(4)), [11, 11, 11, 11])));
        end
        
        function testTransitionDensity2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test that transition density is correct if .C is multiply
            % with 2*pi for identity and with sqrt(2*pi) for sqrt
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            noiseDistribution = ToroidalWNDistribution([1; 2], [1, -0.3; -0.3, 4]);
            coeffsPerDim = 71;
            f = @a2D;
            % Initialize variable like this to make sure that fTrans can be
            % used as in HypertoroidalFourierFilter
            this = struct('hfd', []);
            this.hfd = HypertoroidalFourierDistribution(diag([zeros(1, 50), 1 / (2 * pi)^2, zeros(1, 50)]), 'identity');
            
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdTransId = HypertoroidalFourierDistribution.fromFunction( ...
                @fTrans, coeffsPerDim*ones(1, 2*noiseDistribution.dim), 'identity');
            hfdTransSqrt = HypertoroidalFourierDistribution.fromFunction( ...
                @fTrans, coeffsPerDim*ones(1, 2*noiseDistribution.dim), 'sqrt');
            fixture.teardown;
            
            hfdTransId.C = (2 * pi)^noiseDistribution.dim * hfdTransId.C;
            hfdTransSqrt.C = sqrt(2*pi)^noiseDistribution.dim * hfdTransSqrt.C;
            
            % Denser grid too expensive
            [a, b, c, d] = ndgrid([0, pi], [0, pi], [0, pi], [0, pi]);
            fTransVals = reshape(fTrans(a, b, c, d), 1, []);
            % Flatten tensor to dim x (numel(varargin{1}) tensor for
            % evaluation in pdf of the noise. Because the values for
            % the different dimensions are iterated through last, use
            % reshape into numel(a) x dim first and
            % transpose. Make sure to make the concatenation this
            % dependent on ndims(a) because the dimensionality of the
            % input may be higher than the dimensionality of the distribution!
            hfdTransIdVals = hfdTransId.pdf(reshape(cat(ndims(a), a, b, c, d), [], 2*noiseDistribution.dim)');
            hfdTransSqrtVals = hfdTransSqrt.pdf(reshape(cat(ndims(a), a, b, c, d), [], 2*noiseDistribution.dim)');
            
            
            testCase.verifyEqual(fTransVals, hfdTransIdVals, 'AbsTol', 5E-7);
            testCase.verifyEqual(fTransVals, hfdTransSqrtVals, 'AbsTol', 2E-6);
            
            
            function [x1Propagated, x2Propagated] = a2D(x1, x2)
                % 2-D Igor's function
                param1 = 4;
                param2 = 3;
                x1Propagated = pi * (sin(sign(x1-pi)/2.*abs(x1-pi).^param1/pi^(param1 - 1)) + 1);
                x2Propagated = pi * (sin(sign(x2-pi)/2.*abs(x2-pi).^param2/pi^(param2 - 1)) + 1);
            end
            % Taken from HypertoroidalFourierFilter. Should be changed if
            % it is changed there.
            function p = fTrans(varargin)
                assert(all(diff(cellfun(@numel, varargin)) == 0), ...
                    'All input arguments need to be equally sized. Use ndgrid to generate appropriate input arguments.');
                fout = cell(1, this.hfd.dim);
                % Second d dimenisons are for x(k), propagate through f.
                [fout{:}] = f(varargin{this.hfd.dim+1:2*this.hfd.dim});
                % Concatenate on the last dimension to ensure that it is iterated
                % through LAST when using reshape.
                ws = cat(ndims(varargin{1}), varargin{1:this.hfd.dim}) - cat(ndims(varargin{1}), fout{:});
                % Flatten tensor to dim x (numel(varargin{1}) tensor for
                % evaluation in pdf of the noise. Because the dimension is
                % iterated through last, first reshape it into
                % (numel(varargin{1}) and then transpose to move the
                % concatenated input arguments to the first dimension.
                wsReshaped = reshape(ws, [], noiseDistribution.dim)';
                pdfvals = noiseDistribution.pdf(wsReshaped);
                p = reshape(pdfvals, size(varargin{1})); % Restore to input size
            end
        end
        
        function testJointDensity2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test that transition density is correct if .C is multiply
            % with 2*pi for identity and with sqrt(2*pi) for sqrt
            global enableExpensive
            if ~islogical(enableExpensive), enableExpensive = true; end
            noiseDistribution = ToroidalWNDistribution([1; 2], [1, -0.3; -0.3, 4]);
            coeffsPerDim = 81;
            coeffVector = coeffsPerDim * ones(1, noiseDistribution.dim);
            f = @a2D;
            % Initialize variable this to make sure that fTrans can be
            % used as in HypertoroidalFourierFilter
            this = struct('hfd', []);
            this.hfd = HypertoroidalFourierDistribution(diag([zeros(1, 50), 1 / (2 * pi)^2, zeros(1, 50)]), 'identity');
            
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdTransId = HypertoroidalFourierDistribution.fromFunction( ...
                @fTrans, coeffsPerDim*ones(1, 2*noiseDistribution.dim), 'identity');
            hfdTransSqrt = HypertoroidalFourierDistribution.fromFunction( ...
                @fTrans, coeffsPerDim*ones(1, 2*noiseDistribution.dim), 'sqrt');
            fixture.teardown;
            
            densityInit = ToroidalWNDistribution([3 * 2 / pi; pi], [3, 0.3; 0.3, 1]);
            hfdInitId = HypertoroidalFourierDistribution.fromDistribution(densityInit, coeffVector, 'identity');
            hfdInitSqrt = HypertoroidalFourierDistribution.fromDistribution(densityInit, coeffVector, 'sqrt');
            
            % Initially first and overwrite C (this ensures no additional
            % normalization takes place)
            hfdJointId = ToroidalFourierDistribution.fromDistribution(HypertoroidalUniformDistribution(4), [coeffVector, coeffVector], 'identity');
            hfdJointSqrt = ToroidalFourierDistribution.fromDistribution(HypertoroidalUniformDistribution(4), [coeffVector, coeffVector], 'sqrt');
            
            % This intermediate step is not performed in HypertoroidalFourierFilter for the identity
            % transformation. Next test case shows it is equal for
            % identity transformation. Scale up fTrans after calculation of
            % joint.
            hfdJointId.C = (2 * pi)^noiseDistribution.dim * convnc(hfdTransId.C, reshape(hfdInitId.C, [ones(1, hfdInitId.dim), size(hfdInitId.C)]), 'same');
            % Keep calculations as close to the code as possible. C of
            % hfdTrans is scaled up before convolution.
            hfdJointSqrt.C = convnc((2 * pi)^(this.hfd.dim / 2)*hfdTransSqrt.C, reshape(hfdInitSqrt.C, [ones(1, hfdInitSqrt.dim), size(hfdInitSqrt.C)]), 'same');
            
            % Denser grid too expensive
            [a, b, c, d] = ndgrid([0, pi], [0, pi], [0, pi], [0, pi]);
            
            % Flatten tensor to dim x (numel(varargin{1}) tensor for
            % evaluation in pdf of the noise. Because the values for
            % the different dimensions are iterated through last, use
            % reshape into numel(a) x dim first and
            % transpose. Make sure to make the concatenation this
            % dependent on ndims(a) because the dimensionality of the
            % input may be higher than the dimensionality of the distribution!
            fJointVals = reshape(fTrans(a, b, c, d), 1, []) .* densityInit.pdf(reshape(cat(ndims(c), c, d), [], noiseDistribution.dim)');
            hfdJointIdVals = hfdJointId.pdf(reshape(cat(ndims(a), a, b, c, d), [], 2*noiseDistribution.dim)');
            hfdJointSqrtVals = hfdJointSqrt.pdf(reshape(cat(ndims(a), a, b, c, d), [], 2*noiseDistribution.dim)');
            
            testCase.verifyEqual(fJointVals, hfdJointIdVals, 'AbsTol', 5E-8);
            testCase.verifyEqual(fJointVals, hfdJointSqrtVals, 'AbsTol', 5E-8);
            
            
            function [x1Propagated, x2Propagated] = a2D(x1, x2)
                % 2-D Igor's function
                param1 = 4;
                param2 = 3;
                x1Propagated = pi * (sin(sign(x1-pi)/2.*abs(x1-pi).^param1/pi^(param1 - 1)) + 1);
                x2Propagated = pi * (sin(sign(x2-pi)/2.*abs(x2-pi).^param2/pi^(param2 - 1)) + 1);
            end
            % Taken from HypertoroidalFourierFilter. Should be changed if
            % it is changed there.
            function p = fTrans(varargin)
                assert(all(diff(cellfun(@numel, varargin)) == 0), ...
                    'All input arguments need to be equally sized. Use ndgrid to generate appropriate input arguments.');
                fout = cell(1, this.hfd.dim);
                % Second d dimenisons are for x(k), propagate through f.
                [fout{:}] = f(varargin{this.hfd.dim+1:2*this.hfd.dim});
                % Concatenate on the last dimension to ensure that it is iterated
                % through LAST when using reshape.
                ws = cat(ndims(varargin{1}), varargin{1:this.hfd.dim}) - cat(ndims(varargin{1}), fout{:});
                % Flatten tensor to dim x (numel(varargin{1}) tensor for
                % evaluation in pdf of the noise. Because the dimension is
                % iterated through last, first reshape it into
                % (numel(varargin{1}) and then transpose to move the
                % concatenated input arguments to the first dimension.
                wsReshaped = reshape(ws, [], noiseDistribution.dim)';
                pdfvals = noiseDistribution.pdf(wsReshaped);
                p = reshape(pdfvals, size(varargin{1})); % Restore to input size
            end
        end
        
        function testTruncatePredictionJointId2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Test that using truncation with .valid when using the
            % identity transformation is valid
            noiseDistribution = ToroidalWNDistribution([1; 2], [1, -0.3; -0.3, 4]);
            coeffsPerDim = 21;
            coeffVector = coeffsPerDim * ones(1, noiseDistribution.dim);
            f = @a2D;
            % Initialize variable this to make sure that fTrans can be
            % used as in HypertoroidalFourierFilter
            this = struct('hfd', []);
            this.hfd = HypertoroidalFourierDistribution(diag([zeros(1, 50), 1 / (2 * pi)^2, zeros(1, 50)]), 'identity');
            
            fixture = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            hfdTransId = HypertoroidalFourierDistribution.fromFunction( ...
                @fTrans, coeffsPerDim*ones(1, 2*noiseDistribution.dim), 'identity');
            fixture.teardown;
            
            densityInit = ToroidalWNDistribution([3 * 2 / pi; pi], [3, 0.3; 0.3, 1]);
            hfdInitId = HypertoroidalFourierDistribution.fromDistribution(densityInit, coeffVector, 'identity');
            
            % Initially first and overwrite C (this ensures no additional
            % normalization takes place)
            hfdJointIdSame = ToroidalFourierDistribution.fromDistribution(HypertoroidalUniformDistribution(4), [coeffVector, coeffVector], 'identity');
            hfdJointIdFull = ToroidalFourierDistribution.fromDistribution(HypertoroidalUniformDistribution(4), [coeffVector, coeffVector], 'identity');
            
            % This intermediate step is not performed in HypertoroidalFourierFilter for the identity
            % transformation. Next test case shows it is equal for
            % identity transformation. Scale up fTrans after calculation of
            % joint.
            hfdJointIdSame.C = (2 * pi)^noiseDistribution.dim * convnc( ...
                hfdTransId.C, reshape(hfdInitId.C, [ones(1, hfdInitId.dim), size(hfdInitId.C)]), 'same');
            hfdJointIdFull.C = (2 * pi)^noiseDistribution.dim * convnc( ...
                hfdTransId.C, reshape(hfdInitId.C, [ones(1, hfdInitId.dim), size(hfdInitId.C)]), 'same');
            
            indices = num2cell((size(hfdJointIdSame.C) + 1)/2);
            indices(1:noiseDistribution.dim) = {':'};
            cPredictedIdUsingSame = (2 * pi)^noiseDistribution.dim * hfdJointIdSame.C(indices{:});
            cPredictedIdUsingFull = (2 * pi)^noiseDistribution.dim * hfdJointIdFull.C(indices{:});
            
            % Using HypertoroidalFourierFilter class
            fourierFilterNl = HypertoroidalFourierFilter(coeffVector, 'identity');
            fourierFilterNl.setState(hfdInitId);
            fourierFilterNl.predictNonlinear(f, noiseDistribution);
            
            cPredictedIdUsingValid = fourierFilterNl.hfd.C;
            
            testCase.verifyEqual(cPredictedIdUsingSame, cPredictedIdUsingValid, 'AbsTol', 1E-16);
            testCase.verifyEqual(cPredictedIdUsingFull, cPredictedIdUsingValid, 'AbsTol', 1E-16);
            
            function [x1Propagated, x2Propagated] = a2D(x1, x2)
                % 2-D Igor's function
                param1 = 4;
                param2 = 3;
                x1Propagated = pi * (sin(sign(x1-pi)/2.*abs(x1-pi).^param1/pi^(param1 - 1)) + 1);
                x2Propagated = pi * (sin(sign(x2-pi)/2.*abs(x2-pi).^param2/pi^(param2 - 1)) + 1);
            end
            % Taken from HypertoroidalFourierFilter. Should be changed if
            % it is changed there.
            function p = fTrans(varargin)
                assert(all(diff(cellfun(@numel, varargin)) == 0), ...
                    'All input arguments need to be equally sized. Use ndgrid to generate appropriate input arguments.');
                fout = cell(1, this.hfd.dim);
                % Second d dimenisons are for x(k), propagate through f.
                [fout{:}] = f(varargin{this.hfd.dim+1:2*this.hfd.dim});
                % Concatenate on the last dimension to ensure that it is iterated
                % through LAST when using reshape.
                ws = cat(ndims(varargin{1}), varargin{1:this.hfd.dim}) - cat(ndims(varargin{1}), fout{:});
                % Flatten tensor to dim x (numel(varargin{1}) tensor for
                % evaluation in pdf of the noise. Because the dimension is
                % iterated through last, first reshape it into
                % (numel(varargin{1}) and then transpose to move the
                % concatenated input arguments to the first dimension.
                wsReshaped = reshape(ws, [], noiseDistribution.dim)';
                pdfvals = noiseDistribution.pdf(wsReshaped);
                p = reshape(pdfvals, size(varargin{1})); % Restore to input size
            end
        end
        
        function testPredictNonlinear2D(testCase)
            noiseDistribution = ToroidalWNDistribution([1; 2], [1, -0.3; -0.3, 4]);
            coeffsPerDim = 51;
            coeffVector = coeffsPerDim * ones(1, noiseDistribution.dim);
            f = @a2D;
            % Initialize variable this to make sure that fTrans can be
            % used as in HypertoroidalFourierFilter
            this = struct('hfd', []);
            this.hfd = HypertoroidalFourierDistribution(diag([0, 1 / (2 * pi)^2, 0]), 'identity');
            
            densityInit = ToroidalWNDistribution([3 * 2 / pi; pi], [3, 0.3; 0.3, 1]);
            
            fPredUsingInt = @(xkk)arrayfun(@(i)integral2(@(xk1s, xk2s) ...
                fTrans(xkk(1, i)*ones(size(xk1s)), xkk(2, i)*ones(size(xk1s)), xk1s, xk2s).* ...
                reshape(densityInit.pdf([xk1s(:)'; xk2s(:)']), size(xk1s)), ...
                0, 2*pi, 0, 2*pi), 1:size(xkk, 2));
            
            [xVals, yVals] = meshgrid(0:0.4:2*pi);
            gridReshaped = [xVals(:)'; yVals(:)'];
            fvalsPredUsingInt = fPredUsingInt(gridReshaped);
            for transformation = {'identity', 'sqrt'}
                fourierFilterNl = HypertoroidalFourierFilter(coeffVector, [transformation{:}]);
                fourierFilterNl.setState(HypertoroidalFourierDistribution.fromDistribution(densityInit, coeffVector, [transformation{:}]));
                fourierFilterNl.predictNonlinear(f, noiseDistribution);
                testCase.verifyEqual(fourierFilterNl.hfd.pdf(gridReshaped), fvalsPredUsingInt, 'AbsTol', 1E-7);
            end
            
            function [x1Propagated, x2Propagated] = a2D(x1, x2)
                % 2-D Igor's function
                param1 = 4;
                param2 = 3;
                x1Propagated = pi * (sin(sign(x1-pi)/2.*abs(x1-pi).^param1/pi^(param1 - 1)) + 1);
                x2Propagated = pi * (sin(sign(x2-pi)/2.*abs(x2-pi).^param2/pi^(param2 - 1)) + 1);
            end
            % Taken from HypertoroidalFourierFilter. Should be changed if
            % it is changed there.
            function p = fTrans(varargin)
                assert(all(diff(cellfun(@numel, varargin)) == 0), ...
                    'All input arguments need to be equally sized. Use ndgrid to generate appropriate input arguments.');
                fout = cell(1, this.hfd.dim);
                % Second d dimenisons are for x(k), propagate through f.
                [fout{:}] = f(varargin{this.hfd.dim+1:2*this.hfd.dim});
                % Concatenate on the last dimension to ensure that it is iterated
                % through LAST when using reshape.
                ws = cat(ndims(varargin{1}), varargin{1:this.hfd.dim}) - cat(ndims(varargin{1}), fout{:});
                % Flatten tensor to dim x (numel(varargin{1}) tensor for
                % evaluation in pdf of the noise. Because the dimension is
                % iterated through last, first reshape it into
                % (numel(varargin{1}) and then transpose to move the
                % concatenated input arguments to the first dimension.
                wsReshaped = reshape(ws, [], noiseDistribution.dim)';
                pdfvals = noiseDistribution.pdf(wsReshaped);
                p = reshape(pdfvals, size(varargin{1})); % Restore to input size
            end
        end
        
        function testGetEstimateMean(testCase)
            dist = ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5);
            for transformation = {'identity', 'sqrt'}
                filter = HypertoroidalFourierFilter([9, 9], [transformation{:}]);
                filter.setState(HypertoroidalFourierDistribution.fromDistribution(dist, length(filter.hfd.C), [transformation{:}]));
                testCase.verifyEqual(filter.getEstimateMean, dist.mu, 'AbsTol', 1E-8);
            end
        end
    end
end

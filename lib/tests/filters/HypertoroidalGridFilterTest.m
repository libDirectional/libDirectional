classdef HypertoroidalGridFilterTest < matlab.unittest.TestCase
    methods(Test)
        function testInitialization(testCase)
            noGridPoints=[23, 23];
            filter = HypertoroidalGridFilter([23, 23],2);
            testCase.verifyClass(filter.getEstimate, 'HypertoroidalGridDistribution');
            testCase.verifySize(filter.getEstimate.getGrid(),[2,prod(noGridPoints)]);
            testCase.verifySize(filter.getEstimate.gridValues,[prod(noGridPoints),1]);
        end
        
        function testSetState(testCase)
            filter = HypertoroidalGridFilter([21, 21, 5]);
            hgd = HypertoroidalGridDistribution.fromDistribution(HypertoroidalWNDistribution([1; 1; 1], eye(3)),[21,21,5]);
            testCase.verifyWarningFree(@()filter.setState(hgd));
            testCase.verifyEqual(filter.gd,hgd);
            % Test warnining
            fullGrid = hgd.getGrid();
            hgdtmp=hgd;
            hgdtmp.grid=fullGrid(:,end);
            hgdtmp.gridValues(end)=[];
            hgdtmp.gridValues = hgdtmp.gridValues(1:end-1);
            
            testCase.verifyClass(filter.gd,'HypertoroidalGridDistribution')
            filter = filter.copy;
            testCase.verifyWarning(@()filter.setState(hgdtmp),'setState:gridDiffers');
            testCase.verifyClass(filter.gd,'HypertoroidalGridDistribution')
        end
        
        function testPredictIdentity(testCase)
            % Only test filter interface, validity of operation is tested
            % in HypertoroidalGridDistribution
            filter = HypertoroidalGridFilter([23, 25]);
            gd1 = HypertoroidalGridDistribution.fromDistribution(ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5), [23, 25]);
            gd2 = HypertoroidalGridDistribution.fromDistribution(ToroidalVMSineDistribution([0; 0], [0.8; 1.5], 0.2), [23, 25]);
            gdRes = gd1.convolve(gd2);
            filter.setState(gd1);
            filter.predictIdentity(gd2);
            testCase.verifyClass(filter.gd, 'HypertoroidalGridDistribution');
            testCase.verifySize(filter.gd.gridValues, size(gd1.gridValues));
            testCase.verifyEqual(filter.gd.getGrid(), gdRes.getGrid());
            testCase.verifyEqual(filter.gd.gridValues, gdRes.gridValues, 'Rel', 1E-5);
            testCase.verifyWarning(@()filter.predictIdentity(HypertoroidalWNDistribution([1; 1], eye(2))), 'PredictIdentity:automaticConversion');
        end
        
        function testUpdateIdentity(testCase)
            % As above, only test filter interface
            z = [3; 4];
            filter = HypertoroidalGridFilter([23, 25]);
            gd1 = HypertoroidalGridDistribution.fromDistribution(ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5), [23, 25]);
            gd2 = HypertoroidalGridDistribution.fromDistribution(ToroidalVMSineDistribution(z, [0.8; 1.5], 0.2), [23, 25]);
            hfd2zeroMean = HypertoroidalGridDistribution.fromDistribution(ToroidalVMSineDistribution([0; 0], [0.8; 1.5], 0.2), [23, 25]);
            gdRes = gd1.multiply(gd2);
            filter.setState(gd1);
            filter.updateIdentity(hfd2zeroMean, z);
            testCase.verifyClass(filter.gd, 'HypertoroidalGridDistribution');
            testCase.verifySize(filter.gd.gridValues, size(gd1.gridValues));
            testCase.verifyEqual(filter.gd.gridValues, gdRes.gridValues, 'AbsTol', 1E-6);
        end
        
        function testUpdateNonlinear2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            z = [1; 3];
            filter = HypertoroidalGridFilter([23, 25]);
            hfd1 = HypertoroidalGridDistribution.fromDistribution(ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5), [23, 25]);
            likelihood = @(z, x)1 ./ (sum(abs(x-z)) + .1);
            filter.setState(hfd1);
            filter.updateNonlinear(likelihood, z);
            testCase.verifyEqual(filter.getEstimateMean, z, 'AbsTol', 0.2); % The mean should stay approximately equal, but not entirely
        end
        
        function testPredictNonlinearForLinear2D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            densityInit = HypertoroidalWNDistribution([3; 3], eye(2));
            fNoiseDist = HypertoroidalWNDistribution([0; 0], [1, 0.1; 0.1, 1]);
            noCoeffs = [15, 15];
            tol = 1E-6;
            gdFilterLin = HypertoroidalGridFilter(noCoeffs);
            gdFilterLin.setState(HypertoroidalGridDistribution.fromDistribution(densityInit, noCoeffs));
            gdFilterNl = HypertoroidalGridFilter(noCoeffs);
            gdFilterNl.setState(HypertoroidalGridDistribution.fromDistribution(densityInit, noCoeffs));

            testCase.applyFixture(SuppressedWarningsFixture('PredictIdentity:automaticConversion'));
            gdFilterLin.predictIdentity(fNoiseDist);
            gdFilterNl.predictNonlinear(@(varargin)varargin{:}, fNoiseDist);
            testCase.verifyEqual(gdFilterNl.gd.getGrid(),gdFilterLin.gd.getGrid());
            testCase.verifyEqual(gdFilterNl.getEstimate.gridValues, gdFilterLin.getEstimate.gridValues, 'AbsTol', tol);

            fNoiseDistShifted = fNoiseDist.shift([1; 1]);
            gdFilterLin.predictIdentity(fNoiseDistShifted);
            gdFilterNl.predictNonlinear(@(x)deal(x+1), fNoiseDist);
            testCase.verifyEqual(gdFilterNl.gd.getGrid(),gdFilterNl.gd.getGrid());
            testCase.verifyEqual(gdFilterNl.getEstimate.gridValues, gdFilterLin.getEstimate.gridValues, 'AbsTol', tol);
        end
        
        function testPredictNonlinear1D(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
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
            filterNl = HypertoroidalGridFilter(noCoeffs);
            filterNl.setState(HypertoroidalGridDistribution.fromDistribution(densityInit, noCoeffs));
            filterNl.predictNonlinear(a, fNoiseDist);
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(filterNl.getEstimate.pdf(xvals), fPredUsingInt(xvals), 'AbsTol', 2E-5);
        end
        
        function testWarningFreePredictNonlinearViaTransitionDensity(testCase)
            filterNl = HypertoroidalGridFilter([11, 11]);
            gd = filterNl.getfTransAsGd(@(x)x,HypertoroidalWNDistribution([1; 1], eye(2)));
            testCase.verifyWarningFree(@()filterNl.predictNonlinearViaTransitionDensity(gd));
        end
        
        function testPredictNonlinear2D(testCase)
            noiseDistribution = ToroidalWNDistribution([1; 2], [1, -0.3; -0.3, 4]);
            coeffsPerDim = 51;
            coeffVector = coeffsPerDim * ones(1, noiseDistribution.dim);
            f = @a2Dvec;
            % Initialize variable this to make sure that fTrans can be
            % used as in HypertoroidalGridFilter
            densityInit = ToroidalWNDistribution([3 * 2 / pi; pi], [3, 0.3; 0.3, 1]);
            
            % fTrans defined as in HypertoroidalGridFilter
            fTrans= @(xkk,xk) noiseDistribution.pdf(xkk-a2Dvec(xk));
            fJoint = @(xkk,xk) fTrans(xkk,xk).*densityInit.pdf(xk);
            fJointMatrixInput2nd = @(xkk,xk1s,xk2s)reshape(fJoint(xkk,[xk1s(:)'; xk2s(:)']),size(xk1s));
            fPredUsingInt = @(xkk)arrayfun(@(i)integral2(@(xk1s, xk2s) ...
                fJointMatrixInput2nd(xkk(:,i),xk1s,xk2s), ...
                0, 2*pi, 0, 2*pi),1:size(xkk,2));
            
            [xVals, yVals] = meshgrid(0:0.4:2*pi);
            gridReshaped = [xVals(:)'; yVals(:)'];
            fvalsPredUsingInt = fPredUsingInt(gridReshaped);

            fourierFilterNl = HypertoroidalGridFilter(coeffVector);
            fourierFilterNl.setState(HypertoroidalGridDistribution.fromDistribution(densityInit, coeffVector));
            fourierFilterNl.predictNonlinear(f, noiseDistribution);
            testCase.verifyEqual(fourierFilterNl.gd.pdf(gridReshaped), fvalsPredUsingInt, 'AbsTol', 1E-7);
            
            function [xPropagated] = a2Dvec(x)
                % 2-D Igor's function
                param1 = 4;
                param2 = 3;
                xPropagated = [pi * (sin(sign(x(1,:)-pi)/2.*abs(x(1,:)-pi).^param1/pi^(param1 - 1)) + 1);...
                    pi * (sin(sign(x(2,:)-pi)/2.*abs(x(2,:)-pi).^param2/pi^(param2 - 1)) + 1)];
            end
        end
        
        function testGetEstimateMean(testCase)
            dist = ToroidalVMSineDistribution([1; 3], [0.3; 0.5], 0.5);
            filter = HypertoroidalGridFilter([9, 9], 2);
            filter.setState(HypertoroidalGridDistribution.fromDistribution(dist, [9,9], 'CartesianProd'));
            testCase.verifyEqual(filter.getEstimateMean, dist.mu, 'AbsTol', 1E-8);
        end
    end
end

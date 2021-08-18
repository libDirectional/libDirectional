classdef SphericalGridFilterTest < matlab.unittest.TestCase
    methods(Test)
        function testSetState(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 1001;
            sgf = SphericalGridFilter(noGridPoints);
            testCase.verifySize(sgf.getEstimate.gridValues, [noGridPoints,1]);
            % Test is is uniform at the beginning
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0,'AbsTol',1E-13);
            fixture.teardown;
            
            M = eye(3);
            Z = [-2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            
            sgdState=SphericalGridDistribution.fromDistribution(bd,noGridPoints);
            testCase.verifyWarningFree(@()sgf.setState(sgdState));
            % Verify that it is no longer a uniform distribution
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyGreaterThan(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0.2);
            fixture.teardown;

            % Verify warnings
            fullGrid = sgdState.getGrid();
            fullGrid(:,end)=[];
            sgdState.grid=fullGrid;
            sgdState.gridValues(end)=[];
            sgfTmp=sgf.copy;
            testCase.verifyWarning(@()sgfTmp.setState(sgdState),'setState:gridDiffers');
            
            testCase.verifyWarning(@()sgf.setState(bd),'setState:nonGrid');
        end
        function testPredictConvergecesToUniform(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 501;
            sgf = SphericalGridFilter(noGridPoints);
            M = eye(3);
            Z = [-2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            sgf.setState(SphericalGridDistribution.fromDistribution(bd,noGridPoints));
            
            % Verify that it is not a uniform distribution
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyGreaterThan(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0.2);
            fixture.teardown;
            % Predict 10 times with VM-distributed noise
            trans = @(xkk,xk)cell2mat(arrayfun(@(i)VMFDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            s2s2 = S2CondS2GridDistribution.fromFunction(trans,noGridPoints,true);
            for i=1:10
                valuesAlternativeFormula = (4*pi/size(sgf.getEstimate.getGrid(),2))*sum(sgf.getEstimate.gridValues'.*s2s2.gridValues,2); % Alternative formula
                sgf.predictNonlinearViaTransitionDensity(s2s2);
                testCase.verifyEqual(sgf.getEstimate.gridValues,valuesAlternativeFormula,'AbsTol',1E-12);
            end
            % Verify that is it approximately uniform now
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0,'AbsTol',1E-3);
        end
        function testUpdateWithVMF(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 101;
            sgf = SphericalGridFilter(noGridPoints);
            
            vmf = VMFDistribution(1/sqrt(2)*[1;-1;0],0.1);
            sgf.updateNonlinear(@(z,x)vmf.pdf(x),NaN(3,1));
            fixture=testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(vmf),0,'AbsTol',5e-06);
            fixture.teardown;
            
            vmfPosterior = vmf.multiply(vmf);
            sgf.updateNonlinear(@(z,x)vmf.pdf(x),NaN(3,1));
            
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(vmfPosterior),0,'AbsTol',1e-12);
        end
        function testUpdateIdentity(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmff = VMFFilter;
            noGridPoints = 301;
            sgf = SphericalGridFilter(noGridPoints);
            
            vmfPrior = VMFDistribution(1/sqrt(2)*[-1;0;1],2);
            vmfMeas =  VMFDistribution([0;-1;0],2);
            
            vmff.setState(vmfPrior);
            sgf.setState(SphericalGridDistribution.fromDistribution(vmfPrior,noGridPoints));
            
            vmff.updateIdentity(VMFDistribution([0;0;1],vmfMeas.kappa),vmfMeas.mu);
            sgf.updateIdentity(VMFDistribution([0;0;1],vmfMeas.kappa),vmfMeas.mu);

            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(vmff.getEstimate),0,'AbsTol',1e-10);
        end
        % Compare prediction with SphericalHarmonicsFilter
        function testPredictNonlinearForPseudoadditive(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 401;
            degree=30;
            densityInit = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[-1;0;1],2),VMFDistribution([0;-1;0],2)},[0.5,0.5]);
            kappa = 1;
            sysNoise = @(mu)VMFDistribution(mu, kappa);
            
            shdInitId = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(densityInit, degree, 'identity');
            shdNoise=SphericalHarmonicsDistributionComplex.fromDistribution(sysNoise([0;0;1]), degree);
            shf = SphericalHarmonicsFilter(degree, 'identity');
            shf.setState(shdInitId);
            shf.predictIdentity(shdNoise);
            
            sgdInit=SphericalGridDistribution.fromDistribution(densityInit,noGridPoints);
            % Using predictNonlinear
            trans=@(xkk,xk)cell2mat(arrayfun(@(i)...
                VMFDistribution(xk(:,i),kappa).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            s2s2=S2CondS2GridDistribution.fromFunction(trans,noGridPoints,true);
            
            sgfNonlinear = SphericalGridFilter(noGridPoints);
            sgfNonlinear.setState(sgdInit);
            % Since s2s2 is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgfNonlinear.predictNonlinearViaTransitionDensity(s2s2));
            
            % Using predictIdentity
            sgfIdentity = SphericalGridFilter(noGridPoints);
            sgfIdentity.setState(sgdInit);
            testCase.verifyWarning(@()...
                sgfIdentity.predictIdentity(sysNoise([0;0;1])),'PredictIdentity:Inefficient');
            
            % Verify results. It is okay to use interpolation.
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            % Spherical harmonics may not be perectly normalized, ignoring
            % warning.
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgfNonlinear.getEstimate.totalVariationDistanceNumerical(shf.getEstimate),0,'AbsTol',2e-4);
            testCase.verifyEqual(sgfIdentity.getEstimate.totalVariationDistanceNumerical(shf.getEstimate),0,'AbsTol',2e-4);
        end
        function testPredictNonlinearViaTransitionDensityNextStateIndependent(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Result of prediction is independent from current state
            noGridpoints = 5001;
            vmfInit = VMFDistribution([0;0;1], 100);
            sgf = SphericalGridFilter(noGridpoints);
            sgf.setState(SphericalGridDistribution.fromDistribution(vmfInit,noGridpoints));
            
            % First prediction with Bingham
            Z = [-7.44;-2.84;0];
            M = [-0.18 -0.172 -0.969;-0.738 -0.627 0.248;-0.65 0.76 -0.0141];
            distPredictionResult = BinghamDistribution(Z,M);
            distPredictionResult.F = distPredictionResult.F*distPredictionResult.integralNumerical;
            
            fTrans = @(xkk,xk)repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            s2s2_fTrans = S2CondS2GridDistribution.fromFunction(fTrans, noGridpoints, true);
            sgf.predictNonlinearViaTransitionDensity(s2s2_fTrans);
                        
            fixture1 = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            fixture2 = testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(distPredictionResult),0,'AbsTol',5e-6);
            fixture1.teardown;
            fixture2.teardown;
            
            % Then prediction with VMF
            distPredictionResult = VMFDistribution([0;-1;0], 100);
            fTrans = @(xkk,xk)repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            s2s2_fTrans = S2CondS2GridDistribution.fromFunction(fTrans, noGridpoints, true);
            % Since s2s2 is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgf.predictNonlinearViaTransitionDensity(s2s2_fTrans));
            
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(distPredictionResult),0,'AbsTol',1e-11);
        end
    end
end
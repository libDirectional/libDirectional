classdef HypersphericalGridFilterTest < matlab.unittest.TestCase
    methods(Test)
        function testSetStateS2(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 1001;
            sgf = HypersphericalGridFilter(noGridPoints, 3);
            testCase.verifySize(sgf.getEstimate.gridValues, [noGridPoints,1]);
            % Test is is uniform at the beginning
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0,'AbsTol',1E-13);
            fixture.teardown;
            
            M = eye(3);
            Z = [-2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            
            sgdState=HypersphericalGridDistribution.fromDistribution(bd,noGridPoints);
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
            testCase.verifyWarningFree(@()sgf.setState(sgdState));
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
            % Verify that it is no longer a uniform distribution
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyGreaterThan(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0.2);
            fixture.teardown;

            % Verify warnings
            fullGrid = sgdState.getGrid();
            sgdState.grid=fullGrid(:,end);
            sgdState.gridValues(end)=[];
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
            sgfTmp = sgf.copy;
            testCase.verifyWarning(@()sgfTmp.setState(sgdState),'setState:gridDiffers');
            testCase.verifyClass(sgfTmp.gd,'HypersphericalGridDistribution')
            testCase.verifyWarning(@()sgf.setState(bd),'setState:nonGrid');
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
        end
        function testSetStateS3(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 1001;
            sgf = HypersphericalGridFilter(noGridPoints, 4);
            testCase.verifySize(sgf.getEstimate.gridValues, [noGridPoints,1]);
            % Test is is uniform at the beginning
            testCase.verifyEqual(sum(abs(diff(sgf.getEstimate.gridValues'))),0);
            
            M = eye(4);
            Z = [-2 -1 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            
            sgdState = HypersphericalGridDistribution.fromDistribution(bd,noGridPoints);
            testCase.verifyWarningFree(@()sgf.setState(sgdState));
            % Verify that it is no longer a uniform distribution
            testCase.verifyGreaterThan(sum(abs(diff(sgf.getEstimate.gridValues'))),5);

            % Verify warnings
            fullGrid = sgdState.getGrid();
            sgdState.grid=fullGrid(:,end);
            sgdState.gridValues(end)=[];
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
            sgfTmp = sgf.copy;
            testCase.verifyWarning(@()sgfTmp.setState(sgdState),'setState:gridDiffers');
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
            testCase.verifyWarning(@()sgf.setState(bd),'setState:nonGrid');
            testCase.verifyClass(sgf.gd,'HypersphericalGridDistribution')
        end
        function testPredictConvergecesToUniformS2(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 501;
            sgf = HypersphericalGridFilter(noGridPoints, 3);
            M = eye(3);
            Z = [-2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            sgf.setState(HypersphericalGridDistribution.fromDistribution(bd,noGridPoints));
            
            % Verify that it is not a uniform distribution
            fixture = testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyGreaterThan(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0.2);
            fixture.teardown;
            % Predict 10 times with VM-distributed noise
            trans = @(xkk,xk)cell2mat(arrayfun(@(i)VMFDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            sdsd = SdCondSdGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set',6);
            for i=1:10
                valuesAlternativeFormula = (4*pi/size(sgf.getEstimate.getGrid(),2))*sum(sgf.getEstimate.gridValues'.*sdsd.gridValues,2); % Alternative formula
                sgf.predictNonlinearViaTransitionDensity(sdsd);
                testCase.verifyEqual(sgf.getEstimate.gridValues,valuesAlternativeFormula,'AbsTol',1E-12);
            end
            % Verify that is it approximately uniform now
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(HypersphericalUniformDistribution(3)),0,'AbsTol',1E-3);
        end
        function testPredictConvergecesToUniformS3(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 1001;
            sgf = HypersphericalGridFilter(noGridPoints, 4);
            M = eye(4);
            Z = [-2 -1 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            sgf.setState(HypersphericalGridDistribution.fromDistribution(bd,noGridPoints));
            
            % Verify that it is not a uniform distribution
            testCase.verifyGreaterThan(sum(abs(diff(sgf.getEstimate.gridValues'))),5);
            % Predict 10 times with VM-distributed noise
            trans = @(xkk,xk)cell2mat(arrayfun(@(i)VMFDistribution(xk(:,i),1).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            sdsd = SdCondSdGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set',8);
            manifoldSize = AbstractHypersphericalDistribution.computeUnitSphereSurface(4);
            for i=1:10
                valuesAlternativeFormula = (manifoldSize/size(sgf.getEstimate.getGrid(),2))*sum(sgf.getEstimate.gridValues'.*sdsd.gridValues,2); % Alternative formula
                sgf.predictNonlinearViaTransitionDensity(sdsd);
                testCase.verifyEqual(sgf.getEstimate.gridValues,valuesAlternativeFormula,'AbsTol',1E-12);
            end
            % Verify that is it approximately uniform now
            testCase.verifyEqual(sum(abs(diff(sgf.getEstimate.gridValues'))),0,'AbsTol',1E-3);
        end
        % No test for S2! Use spherical grid filter instead for better
        % results!
        function testUpdateWithVMFS3(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 1001;
            sgf = HypersphericalGridFilter(noGridPoints, 4);
            
            vmf = VMFDistribution(1/sqrt(3)*[1;-1;0;1],0.1);
            sgf.updateNonlinear(@(z,x)vmf.pdf(x),NaN(4,1));
            fixture=testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.verifyEqual(sgf.gd.gridValues,vmf.pdf(sgf.gd.getGrid())','AbsTol',1e-05);
            fixture.teardown;
            
            vmfPosterior = vmf.multiply(vmf);
            sgf.updateNonlinear(@(z,x)vmf.pdf(x),NaN(4,1));
            
            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.gd.gridValues,vmfPosterior.pdf(sgf.gd.getGrid())','AbsTol',5e-6);
        end
        function testUpdateIdentityS2(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmff = VMFFilter;
            noGridPoints = 301;
            sgf = HypersphericalGridFilter(noGridPoints, 3);
            
            vmfPrior = VMFDistribution(1/sqrt(2)*[-1;0;1],2);
            vmfMeas =  VMFDistribution([0;0;1],2);
            z = [0;-1;0];
            
            vmff.setState(vmfPrior);
            sgf.setState(SphericalGridDistribution.fromDistribution(vmfPrior,noGridPoints));
            
            vmff.updateIdentity(vmfMeas, z);
            sgf.updateIdentity(vmfMeas, z);

            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.getEstimate.totalVariationDistanceNumerical(vmff.getEstimate),0,'AbsTol',1e-10);
        end
        function testUpdateIdentityS3(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            vmff = VMFFilter;
            noGridPoints = 150001;
            sgf = HypersphericalGridFilter(noGridPoints, 4);
            
            vmfPrior = VMFDistribution(1/sqrt(3)*[-1;0;1;1],2);
            vmfMeas =  VMFDistribution([0;-1;0;0],2);
            
            vmff.setState(vmfPrior);
            sgf.setState(HypersphericalGridDistribution.fromDistribution(vmfPrior,noGridPoints));
            
            vmff.updateIdentity(VMFDistribution([0;0;0;1],vmfMeas.kappa),vmfMeas.mu);
            sgf.updateIdentity(VMFDistribution([0;0;0;1],vmfMeas.kappa),vmfMeas.mu);

            testCase.applyFixture(SuppressedWarningsFixture('PDF:UseInterpolated'));
            testCase.applyFixture(SuppressedWarningsFixture('Normalization:notNormalized'));
            testCase.verifyEqual(sgf.gd.gridValues',vmff.getEstimate.pdf(sgf.gd.getGrid()),'AbsTol',1e-4);
        end
        % Compare prediction with SphericalHarmonicsFilter
        function testPredictNonlinearForPseudoadditiveS2(testCase)
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
            
            sgfNonlinear = HypersphericalGridFilter(noGridPoints, 3);
            sgfNonlinear.setState(sgdInit);
            % Since s2s2 is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgfNonlinear.predictNonlinearViaTransitionDensity(s2s2));
            
            % Using predictIdentity
            sgfIdentity = HypersphericalGridFilter(noGridPoints, 3);
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
        function testPredictNonlinearForPseudoadditiveS3(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            noGridPoints = 10000;
            particles = 50000;
            densityInit = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(3)*[-1;0;1;1],2),VMFDistribution([0;-1;0;0],2)},[0.5,0.5]);
            kappa = 1;
            sysNoise = @(mu)VMFDistribution(mu, kappa);
            
            % Ensure rng is always the same
            rngState = rng(10);
            spf = HypersphericalParticleFilter(particles, 4);
            spf.setState(densityInit);
            spf.predictIdentity(sysNoise([0;0;0;1]));
            rng(rngState); % Reset it to current settings
            
            sgdInit = HypersphericalGridDistribution.fromDistribution(densityInit,noGridPoints);
            % Using predictNonlinear
            trans = @(xkk,xk)cell2mat(arrayfun(@(i)...
                VMFDistribution(xk(:,i),kappa).pdf(xkk),1:size(xk,2),'UniformOutput',false));
            sdsd = SdCondSdGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set',8);
            
            sgfNonlinear = HypersphericalGridFilter(noGridPoints, 4);
            sgfNonlinear.setState(sgdInit);
            % Since sdsd is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgfNonlinear.predictNonlinearViaTransitionDensity(sdsd));
            
            % Using predictIdentity
            sgfIdentity = HypersphericalGridFilter(noGridPoints, 4);
            sgfIdentity.setState(sgdInit);
            testCase.verifyWarning(@()...
                sgfIdentity.predictIdentity(sysNoise([0;0;1])),'PredictIdentity:Inefficient');
            
            testCase.verifyEqual(sgfNonlinear.gd.gridValues, sgfIdentity.gd.gridValues);
            testCase.verifyEqual(sgfIdentity.getEstimateMean, spf.getEstimateMean, 'AbsTol', 4e-2);
        end
        function testPredictNonlinearViaTransitionDensityNextStateIndependentS3(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture
            % Result of prediction is independent of the initial density
            % because the transition density is independent of the state
            % xk!
            noGridpoints = 10000;
            vmfInit = VMFDistribution([0;0;0;1], 1);
            sgf = HypersphericalGridFilter(noGridpoints, 4);
            % Ignore warning there because we use a density that is almost
            % certain, which is unusual.
            sgf.setState(HypersphericalGridDistribution.fromDistribution(vmfInit,noGridpoints));
            % FORCE normalization for this test case
            sgf.gd.gridValues = sgf.gd.gridValues/(sgf.gd.getManifoldSize * mean(sgf.gd.gridValues));
            
            % First prediction with Bingham
            Z = [-7.44;-2.84;-1;0];
            M = [-0.1565 0.07999 -0.8749 -0.4513;0.7173 -0.5805 -0.3038 0.2373;-0.6305 -0.7749 0.02211 0.03839;0.2521 -0.2369 0.3766 -0.8594];
            distPredictionResult = BinghamDistribution(Z,M);
            distPredictionResult.F = distPredictionResult.F*distPredictionResult.integralNumerical;
            distPredictionResult = HypersphericalMixture({VMFDistribution([0;0;1;0],1),VMFDistribution([0;0;-1;0],1)},[0.5,0.5]);
            
            fTrans = @(xkk,xk)repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = SdCondSdGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set', 8);
            sgf.predictNonlinearViaTransitionDensity(sdsd_fTrans);
                        
            testCase.verifyEqual(sgf.gd.gridValues', distPredictionResult.pdf(sgf.gd.getGrid()), 'AbsTol', 1e-12);

            % Next transition density also does not depend on current
            % state.
            % FORCE normalization again for robust results
            sgf.gd.gridValues = sgf.gd.gridValues/(sgf.gd.getManifoldSize * mean(sgf.gd.gridValues));
            distPredictionResult = VMFDistribution([0;-1;0;0], 10);
            fTrans = @(xkk,xk)repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = SdCondSdGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set', 8);
            % Since sdsd is normalized, should be normalized as well.
            sgf.predictNonlinearViaTransitionDensity(sdsd_fTrans)
            testCase.verifyEqual(sgf.gd.gridValues', distPredictionResult.pdf(sgf.gd.getGrid()), 'AbsTol', 5e-13);
            
        end
    end
end
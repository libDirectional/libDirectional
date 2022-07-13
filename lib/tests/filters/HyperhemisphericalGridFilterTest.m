classdef HyperhemisphericalGridFilterTest < matlab.unittest.TestCase
    methods(Test)
        function testSetStateS2(testCase)
            rng default
            noGridPoints = 1001;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 3);
            testCase.verifySize(sgf.getEstimate.gridValues, [noGridPoints,1]);
            % Test is is uniform at the beginning
            testCase.verifyEqual(sgf.getEstimate.gridValues,...
                1/sgf.getEstimate.getManifoldSize*ones(size(sgf.getEstimate.gridValues)),'AbsTol',1E-13);
            M = eye(3);
            Z = [-2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            
            sgdState = HyperhemisphericalGridDistribution.fromDistribution(bd,noGridPoints);
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
            testCase.verifyWarningFree(@()sgf.setState(sgdState));
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
            % Verify that it is no longer a uniform distribution
            testCase.verifyGreaterThan(sum(abs(sgf.getEstimate.gridValues-(1/sgf.getEstimate.getManifoldSize))),60);
            % Verify estimate matches a mode of the Bingham
            testCase.verifyEqual(min(vecnorm(sgf.getPointEstimate() - [bd.mode(),-bd.mode()])), 0, 'AbsTol',1e-11);

            % Verify warnings
            fullGrid = sgdState.getGrid();
            sgdState.grid=fullGrid(:,end);
            sgdState.gridValues(end)=[];
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
            sgftmp = sgf.copy;
            testCase.verifyWarning(@()sgftmp.setState(sgdState),'setState:gridDiffers');
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
            testCase.verifyWarning(@()sgf.setState(bd),'setState:nonGrid');
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
        end
        function testSetStateS3(testCase)
            noGridPoints = 1001;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 4);
            testCase.verifySize(sgf.getEstimate.gridValues, [noGridPoints,1]);
            % Test is is uniform at the beginning
            testCase.verifyEqual(sum(abs(diff(sgf.getEstimate.gridValues'))),0);
            
            M = eye(4);
            Z = [-2 -1 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            
            sgdState = HyperhemisphericalGridDistribution.fromDistribution(bd,noGridPoints);
            testCase.verifyWarningFree(@()sgf.setState(sgdState));
            % Verify that it is no longer a uniform distribution
            testCase.verifyGreaterThan(sum(abs(diff(sgf.getEstimate.gridValues'))),5);

            % Verify warnings
            fullGrid = sgdState.getGrid();
            sgdState.grid=fullGrid(:,end);
            sgdState.gridValues(end)=[];
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
            sgfTmp = sgf.copy;
            testCase.verifyWarning(@()sgfTmp.setState(sgdState),'setState:gridDiffers');
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
            testCase.verifyWarning(@()sgf.setState(bd),'setState:nonGrid');
            testCase.verifyClass(sgf.gd,'HyperhemisphericalGridDistribution')
        end
        function testPredictConvergecesToUniformS2(testCase)
            noGridPoints = 501;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 3);
            M = eye(3);
            Z = [-2 -1 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            sgf.setState(HyperhemisphericalGridDistribution.fromDistribution(bd,noGridPoints));
            
            % Verify that it is not a uniform distribution
            testCase.verifyGreaterThan(sum(abs(diff(sgf.getEstimate.gridValues'))),3);
            % Predict 10 times with VM-distributed noise
            trans = @(xkk,xk)2*cell2mat(arrayfun(@(i)...
                HypersphericalMixture({VMFDistribution(xk(:,i),1),VMFDistribution(-xk(:,i),1)},[0.5,0.5]).pdf(xkk),...
                1:size(xk,2),'UniformOutput',false));
            sdsd = SdHalfCondSdHalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set_symm',6);
            manifoldSize = sgf.getEstimate.getManifoldSize;
            for i = 1:10
                valuesAlternativeFormula = (manifoldSize/size(sgf.getEstimate.getGrid(),2))*sum(sgf.getEstimate.gridValues'.*sdsd.gridValues,2); % Alternative formula
                sgf.predictNonlinearViaTransitionDensity(sdsd);
                testCase.verifyEqual(sgf.getEstimate.gridValues,valuesAlternativeFormula,'AbsTol',1E-12);
            end
            % Verify that is it approximately uniform now
            testCase.verifyEqual(sum(abs(diff(sgf.getEstimate.gridValues'))),0,'AbsTol',5e-5);
        end
        function testPredictConvergecesToUniformS3(testCase)
            noGridPoints = 1001;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 4);
            M = eye(4);
            Z = [-2 -1 -0.5 0]';
            bd = BinghamDistribution(Z,M);
            bd.F = bd.F*bd.integralNumerical; % Because normalization for Bingham is only an approximation
            sgf.setState(HyperhemisphericalGridDistribution.fromDistribution(bd,noGridPoints));
            
            % Verify that it is not a uniform distribution
            testCase.verifyGreaterThan(sum(abs(diff(sgf.getEstimate.gridValues'))),5);
            % Predict 10 times with VM-distributed noise
            trans = @(xkk,xk)2*cell2mat(arrayfun(@(i)...
                HypersphericalMixture({VMFDistribution(xk(:,i),1),VMFDistribution(-xk(:,i),1)},[0.5,0.5]).pdf(xkk),...
                1:size(xk,2),'UniformOutput',false));
            sdsd = SdHalfCondSdHalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set',8);
            manifoldSize = sgf.getEstimate.getManifoldSize;
            for i=1:10
                valuesAlternativeFormula = (manifoldSize/size(sgf.getEstimate.getGrid(),2))*sum(sgf.getEstimate.gridValues'.*sdsd.gridValues,2); % Alternative formula
                sgf.predictNonlinearViaTransitionDensity(sdsd);
                testCase.verifyEqual(sgf.getEstimate.gridValues,valuesAlternativeFormula,'AbsTol',1E-12);
            end
            % Verify that is it approximately uniform now
            testCase.verifyEqual(sum(abs(diff(sgf.getEstimate.gridValues'))),0,'AbsTol',1E-3);
        end
        function testUpdateWithVMFS2(testCase)
            noGridPoints = 101;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 3);
            
            % The density is symmetric because of the mean
            % No warning is also given because setState is
            % never called, only multiply (which ignores warnings regarding
            % normalization because a likelihood need not be normalized).
            % It only works because it's 0 along last dimension (so 
            dist = VMFDistribution(1/sqrt(2)*[1;-1;0],0.1);

            sgf.updateNonlinear(@(z,x)dist.pdf(x),NaN(3,1));
            testCase.verifyEqual(0.5*sgf.getEstimate.gridValues', dist.pdf(sgf.getEstimate.getGrid()),'AbsTol',1e-05);
            
            vmfPosterior = dist.multiply(dist);
            sgf.updateNonlinear(@(z,x)dist.pdf(x),NaN(3,1));
            
            testCase.verifyEqual(0.5*sgf.getEstimate.gridValues', vmfPosterior.pdf(sgf.getEstimate.getGrid()),'AbsTol',1e-05);
        end
        function testUpdateWithVMFS3(testCase)
            noGridPoints = 1001;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 4);
            
            % For this evaluation, it does not matter that the density is
            % asymmetric. No warning is also given because setState is
            % never called, only multiply (which ignores warnings regarding
            % normalization because a likelihood need not be normalized).
            dist = VMFDistribution(1/sqrt(3)*[1;-1;1;0],0.1);

            sgf.updateNonlinear(@(z,x)dist.pdf(x),NaN(4,1));
            testCase.verifyEqual(0.5*sgf.getEstimate.gridValues', dist.pdf(sgf.getEstimate.getGrid()),'AbsTol',1e-05);
            
            vmfPosterior = dist.multiply(dist);
            sgf.updateNonlinear(@(z,x)dist.pdf(x),NaN(4,1));
            
            testCase.verifyEqual(0.5*sgf.getEstimate.gridValues', vmfPosterior.pdf(sgf.getEstimate.getGrid()),'AbsTol',1e-05);
        end
        function testUpdateIdentityAndNonlinearS2(testCase)
            vmff = VMFFilter;
            noGridPoints = 3001;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 3);
            
            vmfPrior = VMFDistribution(1/sqrt(2)*[-1;1;0],2);
            vmfMeas =  VMFDistribution([0;0;1],2);
            z = [0;-1;0];
            
            vmff.setState(vmfPrior);
            sgf.setState(HyperhemisphericalGridDistribution.fromDistribution(vmfPrior,noGridPoints));
            
            sgfNonlin = sgf.copy;
            
            vmff.updateIdentity(vmfMeas, z);
            sgf.updateIdentity(vmfMeas, z);
            sgfNonlin.updateNonlinear(@(z,x)VMFDistribution(z, 2).pdf(x), z); % VMF is commutative for x and z.
            
            testCase.verifyEqual(sgfNonlin.getEstimate.gridValues, sgf.getEstimate.gridValues,'AbsTol',1e-15); % Same results for lin and nonlin
            testCase.verifyEqual(0.5*sgf.getEstimate.gridValues', vmff.getEstimate.pdf(sgf.getEstimate.getGrid()),'AbsTol',1e-5);
        end
        function testUpdateIdentityAndNonlinearS3(testCase)
            vmff = VMFFilter;
            noGridPoints = 150001;
            sgf = HyperhemisphericalGridFilter(noGridPoints, 4);
            
            vmfPrior = VMFDistribution(1/sqrt(3)*[-1;1;1;0],2);
            vmfMeas =  VMFDistribution([0;0;0;1],2);
            z = [0;-1;0;0];
            
            vmff.setState(vmfPrior);
            sgf.setState(HyperhemisphericalGridDistribution.fromDistribution(vmfPrior,noGridPoints));
            
            sgfNonlin = sgf.copy;
            
            vmff.updateIdentity(vmfMeas,z);
            sgf.updateIdentity(vmfMeas,z);
            
            sgfNonlin.updateNonlinear(@(z,x)VMFDistribution(z, 2).pdf(x), z); % VMF is commutative for x and z.
            
            testCase.verifyEqual(sgfNonlin.getEstimate.gridValues, sgf.getEstimate.gridValues,'AbsTol',1e-15); % Same results for lin and nonlin
            testCase.verifyEqual(0.5*sgf.getEstimate.gridValues', vmff.getEstimate.pdf(sgf.getEstimate.getGrid()),'AbsTol',1e-5);
        end
        % Compare prediction with SphericalHarmonicsFilter
        function testPredictNonlinearForPseudoadditiveS2(testCase)
            noGridPoints = 2001;
            degree = 50;
            densityInit = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(2)*[1;1;0],2),VMFDistribution(1/sqrt(2)*[-1;-1;0],2)},[0.5,0.5]);
            kappa = 1;
            sysNoise = @(mu)HypersphericalMixture({VMFDistribution(mu, kappa),VMFDistribution(-mu, kappa)},[0.5,0.5]);
            
            shdInitId = SphericalHarmonicsDistributionComplex.fromDistributionNumericalFast(densityInit, degree, 'identity');
            shdNoise = SphericalHarmonicsDistributionComplex.fromDistribution(sysNoise([0;0;1]), degree);
            shf = SphericalHarmonicsFilter(degree, 'identity');
            shf.setState(shdInitId);
            shf.predictIdentity(shdNoise);
            
            sgdInit = HyperhemisphericalGridDistribution.fromDistribution(densityInit,noGridPoints);
            % Using predictNonlinear
            trans=@(xkk,xk)cell2mat(arrayfun(@(i)... % * 0.5 and * 2 cancel out
                VMFDistribution(xk(:,i),kappa).pdf(xkk)+VMFDistribution(xk(:,i),kappa).pdf(-xkk),1:size(xk,2),'UniformOutput',false));
            sdhalfsdhalf=SdHalfCondSdHalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set_symm',6);
            
            sgfNonlinear = HyperhemisphericalGridFilter(noGridPoints, 3);
            sgfNonlinear.setState(sgdInit);
            % Since sdHalfsdHalf is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgfNonlinear.predictNonlinearViaTransitionDensity(sdhalfsdhalf));
            
            % Using predictIdentity
            sgfIdentity = HyperhemisphericalGridFilter(noGridPoints, 3);
            sgfIdentity.setState(sgdInit);
            testCase.verifyWarning(@()...
                sgfIdentity.predictIdentity(sysNoise([0;0;1])),'PredictIdentity:Inefficient');
            
            testCase.verifyEqual(0.5*sgfNonlinear.getEstimate.gridValues', shf.getEstimate.pdf(sgfNonlinear.getEstimate.getGrid()),'AbsTol',1e-5);
            testCase.verifyEqual(0.5*sgfIdentity.getEstimate.gridValues', shf.getEstimate.pdf(sgfIdentity.getEstimate.getGrid()),'AbsTol',1e-5);
        end
        function testPredictNonlinearForPseudoadditiveS3VersusHypersphericalGF(testCase)
            noGridPoints = 101;
            densityInit = HypersphericalMixture(...
                {VMFDistribution(1/sqrt(3)*[1;1;0;-1],2),VMFDistribution(1/sqrt(3)*[-1;-1;0;1],2)},[0.5,0.5]);
            kappa = 1;
            sysNoise = @(mu)HypersphericalMixture({VMFDistribution(mu, kappa),VMFDistribution(-mu, kappa)},[0.5,0.5]);

            sgdFullInit = HypersphericalGridDistribution.fromDistribution(densityInit, 2*noGridPoints, 'eq_point_set_symm');
            % Using predictNonlinear
            trans=@(xkk,xk)cell2mat(arrayfun(@(i)...
                0.5 * (VMFDistribution(xk(:,i),kappa).pdf(xkk)+VMFDistribution(xk(:,i),kappa).pdf(-xkk)),1:size(xk,2),'UniformOutput',false));
            sdsd=SdCondSdGridDistribution.fromFunction(trans,2*noGridPoints,true,'eq_point_set_symm',8);
            
            sgfFullNonlinear = HypersphericalGridFilter(2*noGridPoints, 4, 'eq_point_set_symm');
            sgfFullNonlinear.setState(sgdFullInit);
            % Since sdHalfsdHalf is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgfFullNonlinear.predictNonlinearViaTransitionDensity(sdsd));
            
            sgdInit = HyperhemisphericalGridDistribution.fromDistribution(densityInit,noGridPoints);
            % Using predictNonlinear
            trans=@(xkk,xk)cell2mat(arrayfun(@(i)... % * 0.5 and * 2 cancel out
                VMFDistribution(xk(:,i),kappa).pdf(xkk)+VMFDistribution(xk(:,i),kappa).pdf(-xkk),1:size(xk,2),'UniformOutput',false));
            sdhalfsdhalf=SdHalfCondSdHalfGridDistribution.fromFunction(trans,noGridPoints,true,'eq_point_set_symm',8);
            
            sgfNonlinear = HyperhemisphericalGridFilter(noGridPoints, 4);
            sgfNonlinear.setState(sgdInit);
            % Since sdHalfsdHalf is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgfNonlinear.predictNonlinearViaTransitionDensity(sdhalfsdhalf));
            
            % Using predictIdentity
            sgfIdentity = HyperhemisphericalGridFilter(noGridPoints, 4);
            sgfIdentity.setState(sgdInit);
            testCase.verifyWarning(@()...
                sgfIdentity.predictIdentity(sysNoise([0;0;0;1])),'PredictIdentity:Inefficient');
            
            testCase.verifyEqual(0.5*sgfNonlinear.gd.gridValues,sgfFullNonlinear.gd.gridValues(1:noGridPoints),'AbsTol',1e-16);
            testCase.verifyEqual(0.5*sgfIdentity.gd.gridValues,sgfFullNonlinear.gd.gridValues(1:noGridPoints),'AbsTol',1e-16);
        end
        function testPredictNonlinearViaTransitionDensityNextStateIndependentS2(testCase)
            % Result of prediction is independent of the initial density
            % because the transition density is independent of the state
            % xk!
            noGridpoints = 5001;
            vmfInit = VMFDistribution([0;1;0], 1);
            sgf = HyperhemisphericalGridFilter(noGridpoints, 3);
            sgf.setState(HyperhemisphericalGridDistribution.fromDistribution(vmfInit, noGridpoints));
            
            % First prediction with Bingham
            Z = [-7.44;-2.84;0];
            M = [-0.18 -0.172 -0.969;-0.738 -0.627 0.248;-0.65 0.76 -0.0141];
            distPredictionResult = BinghamDistribution(Z,M);
            distPredictionResult.F = distPredictionResult.F*distPredictionResult.integralNumerical;
            
            fTrans = @(xkk,xk)2*repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = SdHalfCondSdHalfGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set', 6);
            sgf.predictNonlinearViaTransitionDensity(sdsd_fTrans);
                        
            testCase.verifyEqual(0.5*sgf.gd.gridValues', distPredictionResult.pdf(sgf.gd.getGrid()),'AbsTol',1e-5);
            
            % Then prediction with a transition density that is independent
            % of the current estimate
            distPredictionResult = VMFDistribution([0;-1;0], 100);
            fTrans = @(xkk,xk)2*repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = SdHalfCondSdHalfGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set', 6);
            % Since sdsd is normalized, should be normalized as well.
            testCase.verifyWarningFree(@()...
                sgf.predictNonlinearViaTransitionDensity(sdsd_fTrans));
            
            testCase.verifyEqual(0.5*sgf.gd.gridValues', distPredictionResult.pdf(sgf.gd.getGrid()),'AbsTol',1e-3);
        end
        function testPredictNonlinearViaTransitionDensityNextStateIndependentS3(testCase)
            % Result of prediction is independent of the initial density
            % because the transition density is independent of the state
            % xk!
            noGridpoints = 10000;
            vmfInit = VMFDistribution([0;0;1;0], 1);
            sgf = HyperhemisphericalGridFilter(noGridpoints, 4);
            % Ignore warning there because we use a density that is almost
            % certain, which is unusual.
            sgf.setState(HyperhemisphericalGridDistribution.fromDistribution(vmfInit,noGridpoints));
            % FORCE normalization for this test case
            sgf.gd.gridValues = sgf.gd.gridValues/(sgf.gd.getManifoldSize * mean(sgf.gd.gridValues));
            
            % First prediction with Bingham
            Z = [-7.44;-2.84;-1;0];
            M = [-0.1565 0.07999 -0.8749 -0.4513;0.7173 -0.5805 -0.3038 0.2373;-0.6305 -0.7749 0.02211 0.03839;0.2521 -0.2369 0.3766 -0.8594];
            distPredictionResult = BinghamDistribution(Z,M);
            distPredictionResult.F = distPredictionResult.F*distPredictionResult.integralNumerical;
            distPredictionResult = HypersphericalMixture({VMFDistribution([0;0;1;0],1),VMFDistribution([0;0;-1;0],1)},[0.5,0.5]);
            
            fTrans = @(xkk,xk)2*repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = SdHalfCondSdHalfGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set', 8);
            sgf.predictNonlinearViaTransitionDensity(sdsd_fTrans);
                        
            testCase.verifyEqual(0.5*sgf.gd.gridValues', distPredictionResult.pdf(sgf.gd.getGrid()), 'AbsTol', 1e-12);

            % Next transition density also does not depend on current
            % state.
            % FORCE normalization again for robust results
            sgf.gd.gridValues = sgf.gd.gridValues/(sgf.gd.getManifoldSize * mean(sgf.gd.gridValues));
            distPredictionResult = VMFDistribution([0;-1;0;0], 10);
            fTrans = @(xkk,xk)2*repmat(distPredictionResult.pdf(xkk)',[1,size(xk,2)]);
            sdsd_fTrans = SdHalfCondSdHalfGridDistribution.fromFunction(fTrans, noGridpoints, true, 'eq_point_set', 8);
            % Since sdsd is normalized, should be normalized as well.
            sgf.predictNonlinearViaTransitionDensity(sdsd_fTrans)
            testCase.verifyEqual(0.5*sgf.gd.gridValues', distPredictionResult.pdf(sgf.gd.getGrid()), 'AbsTol', 5e-13);
        end
    end
end
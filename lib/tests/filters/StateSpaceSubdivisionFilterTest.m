classdef StateSpaceSubdivisionFilterTest < matlab.unittest.TestCase
    methods(Test)
        function testInitialization(testCase)
            testCase.verifyWarningFree(@()StateSpaceSubdivisionFilter());
        end
        function testSetState(testCase)
            n=10;
            apd = StateSpaceSubdivisionGaussianDistribution(...
                HyperhemisphericalGridDistribution.fromDistribution(HyperhemisphericalUniformDistribution(3),n),...
                repmat(GaussianDistribution([0;0], 100*eye(2)),[1,n]));
            
            lpf = StateSpaceSubdivisionFilter;
            testCase.verifyWarningFree(@()lpf.setState(apd));
        end
        function testPrediction(testCase)
            n = 10;
            gdInit = HyperhemisphericalGridDistribution.fromDistribution(HyperhemisphericalUniformDistribution(3),n);
            apd = StateSpaceSubdivisionGaussianDistribution(gdInit,repmat(GaussianDistribution([1;1], 100*eye(2)),[1,n]));
            
            kappaSys = 1;
            fTrans = @(xkk,xk)cell2mat(arrayfun(@(i)... % * 0.5 and * 2 cancel out
                    VMFDistribution(xk(:,i),kappaSys).pdf(xkk)...
                    +VMFDistribution(xk(:,i),kappaSys).pdf(-xkk),1:size(xk,2),'UniformOutput', false));
              
             fTrans_sdhalfxsdhalf = ...
                SdHalfCondSdHalfGridDistribution.fromFunction(...
                fTrans, n, true, 'eq_point_set_symm',2*apd.gd.dim);
            
            lpf = StateSpaceSubdivisionFilter;
            lpf.setState(apd)
            
            % First, we consider a model that does not change anything
            testCase.verifyWarning(@()lpf.predictLinear([]),'StateSpaceSubdivisionFilter:NoParamsForPred');
            % Verify nothing has changed
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),ones(2,1,n),'AbsTol',1E-15);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(100*eye(2),[1,1,n]),'AbsTol',1E-30);
            testCase.verifyEqual(lpf.getEstimate().gd.gridValues, apd.gd.gridValues);
            
            % Now, we consider a model that quadruples the
            % uncertainties of the linear part in one dimension and reduces
            % those along the other
            lpf.predictLinear([],[],repmat(diag([2,0.5]),[1,1,n,1]));
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat([2;0.5],[1,1,n]),'AbsTol',1E-15);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(diag([400,25]),[1,1,n]),'AbsTol',1E-30);
            testCase.verifyEqual(lpf.getEstimate().gd.gridValues, apd.gd.gridValues);
            
            % Now, we consider a model that increases the uncertainties
            % by 1 without changing anything else
            lpf.predictLinear([],repmat(eye(2),[1,1,n,1]),[]);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat([2;0.5],[1,1,n]),'AbsTol',1E-12);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(diag([401,26]),[1,1,n]),'AbsTol',1E-12);
            testCase.verifyEqual(lpf.getEstimate().gd.gridValues, apd.gd.gridValues);
            
            % Now, a grid distribution is used for the angular part. The Gaussians are mixed
            % from the previous Gaussians. Since all of the are identical,
            % we should get out the same ones
            lpf.predictLinear(fTrans_sdhalfxsdhalf,[],[]);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat([2;0.5],[1,1,n]),'AbsTol',1E-12);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(diag([401,26]),[1,1,n]),'AbsTol',1E-12);
            
            lpf.predictLinear([],[],[],reshape([1:n;-(1:n)],2,1,n));
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),reshape([1:n;-(1:n)],2,1,n)+[2;0.5],'AbsTol',1E-12);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(diag([401,26]),[1,1,n]),'AbsTol',1E-12);
            
            filterNormal = HyperhemisphericalGridFilter(n,3,'eq_point_set_symm');
            filterNormal.setState(gdInit);
            filterNormal.predictNonlinearViaTransitionDensity(fTrans_sdhalfxsdhalf);
            testCase.verifyEqual(lpf.getEstimate().gd.gridValues, filterNormal.gd.gridValues,'AbsTol',1E-15);
        end
        function testPredictionConvergeToUniform(testCase)
            % If we repreatedly predict, we should converge to a uniform
            % distribution on the circle
            n = 11;
            hgd = HypertoroidalGridDistribution.fromDistribution(VMDistribution(1,1),n);
            apd = HypercylindricalStateSpaceSubdivisionGaussianDistribution(hgd,repmat(GaussianDistribution(1,1),[1,n]));

            kappa = 1;
            noiseAng = VMDistribution(0,kappa);
            alphakks = hgd.getGrid();
            alphaks = alphakks;

            gridGen = @(xk)alphakks'-alphaks-xk;
            transValsGen = @(xk)reshape(noiseAng.pdf(reshape(gridGen(xk),1,[])),[size(hgd.getGrid(),2),size(hgd.getGrid(),2)]);
            tdcondtd = TdCondTdGridDistribution(hgd.getGrid(),transValsGen(1));

            lppf = StateSpaceSubdivisionFilter;
            lppf.setState(apd);
            
            for i=1:20
                lppf.predictLinear(tdcondtd,ones(1,1,n,n),[])
            end
            % Verify we have a uniform density on the circle
            testCase.verifyEqual(lppf.getEstimate().gd.gridValues,1/(2*pi)*ones(n,1),'AbsTol',1e-7);
            % Since all are the same, the reduced Gaussians should have mu = 1 and the uncertainty should be 1*20*1
            testCase.verifyEqual([lppf.getEstimate().linearDistributions.mu], ones(1,n), 'AbsTol', 1e-14);
            testCase.verifyEqual([lppf.getEstimate().linearDistributions.C], 21*ones(1,n), 'AbsTol', 1e-14);
        end
        
        function testUpdate(testCase)
            n = 10;
            gdInit = HemisphericalGridDistribution.fromDistribution(HemisphericalUniformDistribution(),n);
            apd = StateSpaceSubdivisionGaussianDistribution(gdInit,repmat(GaussianDistribution([1;1], 100*eye(2)),[1,n]));
            
            lpf = StateSpaceSubdivisionFilter;
            lpf.setState(apd)
            
            % We update the distributions with covariance 50 on the diagonal with
            % the same one
            lpf.update([],apd.linearDistributions);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat([1;1],[1,1,n]))
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(50*eye(2),[1,1,n]))
            
            % Now update with a nondiagonal one. Show that we get the same
            % result as when multiplying the Gaussians
            gaussianLikelihood = GaussianDistribution([3;5],[25,5;5,25]);
            posterior = lpf.apd.linearDistributions(1).multiply(gaussianLikelihood);
            lpf.update([],gaussianLikelihood);
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat(posterior.mu,[1,1,n]),'RelTol',1e-15)
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(posterior.C,[1,1,n]),'RelTol',1e-15)
            
            % Now update with a VMF and no linear part. The linear part
            % should stay the same.
            % fromDistribution:asymmetricOnHypersphere
            vmfFullSphere = VMFDistribution([0;1;1]/sqrt(2),1);
            noise = CustomHemisphericalDistribution(@(x)vmfFullSphere.pdf(x)+vmfFullSphere.pdf(-x));
            lpf.update(noise)
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat(posterior.mu,[1,1,n]),'RelTol',1e-15)
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(posterior.C,[1,1,n]),'RelTol',1e-15)
            
            % Should be the same as directly approximating the density
            % (enforce normalization to be equal to the normalized result
            % stemming from the update step)
            gdForcedNorm = HemisphericalGridDistribution.fromDistribution(noise,n).normalize(tol=0);
            testCase.verifyEqual(lpf.apd.gd.gridValues,gdForcedNorm.gridValues,'RelTol',5e-16);
            
            % Now update with a VMF and Gaussian with all Gaussians equal.
            % This should be the same as updating them separately because
            % the influence of the Gaussians is the same on all.
            likelihood = GaussianDistribution([1;-1],[15,-0.5;-0.5,15]);
            posterior = arrayfun(@(likelihood)lpf.apd.linearDistributions(1).multiply(likelihood),likelihood);
            lpf.update(noise,likelihood)
            
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),repmat(posterior.mu,[1,1,n]),'RelTol',5e-15)
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),repmat(posterior.C,[1,1,n]),'RelTol',5e-15)
            posteriorNumerical = CustomHemisphericalDistribution(@(x)noise.pdf(x).*noise.pdf(x));
            gdForcedNorm = HemisphericalGridDistribution.fromDistribution(posteriorNumerical,n).normalize(tol=0);
            
            testCase.verifyEqual(lpf.apd.gd.gridValues,gdForcedNorm.gridValues,'RelTol',5e-16);
            
            % Final test, which is hard to validate: with different linear
            % distributions
            likelihoods = arrayfun(@(i)GaussianDistribution([i;-i],[15,i/10;i/10,15]),1:n);
            posteriors = arrayfun(@(likelihood)lpf.apd.linearDistributions(1).multiply(likelihood),likelihoods);
            lpf.update(noise,likelihoods)
            
            % The linear distributions should not be influenced by the
            % nonlinear one
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.mu),cat(3,posteriors.mu),'RelTol',5e-15)
            testCase.verifyEqual(cat(3,lpf.apd.linearDistributions.C),cat(3,posteriors.C),'RelTol',5e-15)
            % Assert that is has gotten more certain now (is less like a
            % uniform distribution)
            testCase.verifyGreaterThan(sum(abs(lpf.apd.gd.gridValues-(1/lpf.apd.gd.getManifoldSize()))),...
                sum(abs(gdForcedNorm.gridValues-(1/lpf.apd.gd.getManifoldSize()))));
            
        end
    end
end

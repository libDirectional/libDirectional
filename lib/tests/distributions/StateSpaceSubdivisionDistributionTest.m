classdef StateSpaceSubdivisionDistributionTest < matlab.unittest.TestCase
    methods(Test)
        function testConstructor(testCase)
            n = 100;
            rbd = StateSpaceSubdivisionDistribution(HyperhemisphericalGridDistribution.fromDistribution(...
                HyperhemisphericalUniformDistribution(3),n),...
                repmat(GaussianDistribution([0;0;0], 1000*eye(3)),[n,1]));
            testCase.verifySize(rbd.linearDistributions,[100,1]);
            testCase.verifySize(rbd.gd.gridValues,[100,1]);
        end
        
        function testPlotSxR(testCase)
            close all
            n = 100;
            means = rand(n,1);
            vars = 1+rand(n,1);
            linearDistributions = arrayfun(@(mean,var)GaussianDistribution(mean,var),means,vars);
            gd = FIGDistribution.fromDistribution(...
                CircularUniformDistribution,n);
            sed = StateSpaceSubdivisionDistribution(gd, linearDistributions);
            testCase.verifyWarningFree(@()sed.plot);
            close all
            testCase.verifyWarningFree(@()sed.plot(true));
        end
        
        function testPlotH2xR(testCase)
            close all
            n = 100;
            means = rand(n,1);
            vars = 1+rand(n,1);
            linearDistributions = arrayfun(@(mean,var)GaussianDistribution(mean,var),means,vars);
            gd = HyperhemisphericalGridDistribution.fromDistribution(...
                HyperhemisphericalUniformDistribution(3),n);
            rbd = StateSpaceSubdivisionDistribution(gd, linearDistributions);
            testCase.verifyWarningFree(@()rbd.plot);
            testCase.verifyWarning(@()rbd.plot(true),'StateSpaceSubdivisionDistribution:invalidPlotArgument');
        end
        
        function testPlotH2xR2(testCase)
            close all
            n = 100;
            means = rand(2,n);
            vars = eye(2)+rand(1,1,n);
            linearDistributions = arrayfun(@(i)GaussianDistribution(means(:,i),vars(:,:,i)),1:n);
            gd = HyperhemisphericalGridDistribution.fromDistribution(...
                HyperhemisphericalUniformDistribution(3),n);
            rbd = StateSpaceSubdivisionDistribution(gd, linearDistributions);
            h = rbd.plot(); % Not sure why this sometimes throws a warning
            testCase.verifyInstanceOf(h,'handle'); 
            
            testCase.verifyWarning(@()rbd.plot(true),'StateSpaceSubdivisionDistribution:invalidPlotArgument');
        end

        function testModeSxR(testCase)
            n = 100;
            means = rand(n,1);
            vars = ones(n,1);
            trueInd = 3;
            vars(trueInd)=0.5;
            linearDistributions = arrayfun(@(mean,var)GaussianDistribution(mean,var),means,vars);
            gd = FIGDistribution.fromDistribution(...
                CircularUniformDistribution,n);
            sed = StateSpaceSubdivisionDistribution(gd, linearDistributions);
            testCase.verifyEqual(sed.mode, [2*pi/100*(trueInd-1);means(trueInd)]);
        end
        
        function testFromDistributionUncorrelated(testCase)
            n = 1000;
            vm = VMDistribution(0,1);
            gauss = GaussianDistribution(0,1);
            customUncorrelated = CustomHypercylindricalDistribution(...
                @(x)vm.pdf(x(1,:)).*gauss.pdf(x(2,:)),1,1);
            
            hcrbd = StateSpaceSubdivisionDistribution.fromDistribution(...
                customUncorrelated,n);
            
            [xMesh, yMesh] = ndgrid(-4:0.1:4, -4:0.1:4);

            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'nearestNeighbor','nearestNeighbor'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',0.003);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'nearestNeighbor','gridDefault'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',2e-15);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'mixture','nearestNeighbor'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',0.003);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'mixture','gridDefault'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',1e-14);
        end
        
        function testFromDistributionUncorrelatedToGaussians(testCase)
            n = 1000;
            vm = VMDistribution(0,1);
            gauss = GaussianDistribution(0,1);
            customUncorrelated = CustomHypercylindricalDistribution(...
                @(x)vm.pdf(x(1,:)).*gauss.pdf(x(2,:)),1,1);
            
            hcrbd = StateSpaceSubdivisionDistribution.fromDistribution(...
                customUncorrelated,n);
            dists = arrayfun(@(dist)GaussianDistribution(dist.mean(),dist.covariance()),hcrbd.linearDistributions);
            hcrbd.linearDistributions = dists;
            
            [xMesh, yMesh] = ndgrid(-4:0.1:4, -4:0.1:4);

            % The approximation with a Gaussian may not be optimal so
            % higher tolerances can be necessary (compared with the
            % previous test case)
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'nearestNeighbor','nearestNeighbor'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',0.003);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'nearestNeighbor','gridDefault'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',2e-10);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'mixture','nearestNeighbor'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',0.003);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'mixture','gridDefault'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',5e-12);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'convexCombinationGaussian','nearestNeighbor'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',0.003);
            testCase.verifyEqual(hcrbd.pdf([xMesh(:)';yMesh(:)'],'convexCombinationGaussian','gridDefault'),...
                customUncorrelated.pdf([xMesh(:)';yMesh(:)']),'RelTol',5e-12);
        end
        
        function testFromFunctionH2xR(testCase)
            rng default
            hud = HemisphericalUniformDistribution();
            gauss = GaussianDistribution(0, 1);
            cd = CustomHypercylindricalDistribution(@(x)hud.pdf(x(1:3,:)).*gauss.pdf(x(4,:)),3,1);

            apd = StateSpaceSubdivisionDistribution.fromFunction(@(x)cd.pdf(x),101,3,1,'hyperhemisphere');

            pointsBounded = randn(3,100);
            pointsCart = [pointsBounded./vecnorm(pointsBounded);randn(1,100)];

            testCase.verifyEqual(apd.pdf(pointsCart,'nearestNeighbor','nearestNeighbor'),cd.pdf(pointsCart),'AbsTol',5e-17);
            testCase.verifyEqual(apd.pdf(pointsCart,'nearestNeighbor','gridDefault'),cd.pdf(pointsCart),'AbsTol',5e-17);
        end
    end
        
end
